"""
水分提取模块
计算植物根系从土壤中提取水分
"""
module WaterExtraction

using ..SoilParameters

export extraction


const POTLEAF = -153.0  # 叶片水势 (m)
const POTWILT = -153.0  # 凋萎点水势 (m)  
const POTFC = -3.366    # 田间持水量水势 (m)

"""
从土壤层中提取水分用于蒸腾

# 参数
- `i::Int`, `j::Int`: 网格坐标
- `nzg::Int`: 土壤层数
- `slz::Vector{Float64}`: 土壤层深度 (m) [nzg+1]
- `dz::Vector{Float64}`: 土壤层厚度 (m) [nzg]
- `Δt::Float64`: 时间步长 (s)
- `soiltxt::Int`: 土壤类型
- `wtd::Float64`: 地下水位深度 (m)
- `smoi::Vector{Float64}`: 土壤含水量 [nzg]
- `smoiwtd::Float64`: 地下水含水量
- `delta::Float64`, `gamma::Float64`, `lambda::Float64`: 气象参数
- `lai::Float64`: 叶面积指数
- `ra_a::Float64`, `ra_c::Float64`: 空气动力学阻力
- `rs_c_factor::Float64`: 冠层阻力因子
- `R_a::Float64`, `R_s::Float64`: 组合阻力
- `petfactor_s::Float64`, `petfactor_c::Float64`: 蒸散发因子
- `inactivedays::Vector{Int}`: 非活跃天数 [nzg+2]
- `maxinactivedays::Int`: 最大非活跃天数
- `fieldcp::Matrix{Float64}`: 田间持水量 [nzg×nstyp]
- `hhveg::Float64`: 植被高度 (m)
- `fdepth::Float64`: 根系深度因子 (m)
- `icefac::Vector{Int8}`: 冰冻因子 [nzg]

# 返回
- `Tuple`: (pet_s, pet_c, watdef, dsmoi, dsmoideep)
"""
function extraction(i::Int, j::Int, nzg::Int, slz::Vector{Float64}, dz::Vector{Float64}, 
                   Δt::Float64, soiltxt::Int, wtd::Float64, smoi::Vector{Float64}, 
                   smoiwtd::Float64, δ::Float64, γ::Float64, λ::Float64, 
                   lai::Float64, ra_a::Float64, ra_c::Float64, rs_c_factor::Float64,
                   R_a::Float64, R_s::Float64, petfactor_s::Float64, petfactor_c::Float64,
                   inactivedays::Vector{Int}, maxinactivedays::Int, fieldcp::Matrix{Float64},
                   hhveg::Float64, fdepth::Float64, icefac::Vector{Int8})
    
    
    
    hveg = 2.0 * hhveg / 3.0
    
    # 初始化变量
    easy = zeros(Float64, nzg)
    easydeep = 0.0
    dzwtd = 0.0
    rootmask = zeros(Int, nzg)
    dz2 = copy(dz)
    dz3 = 0.0
    
    # 计算地下水位所在层
    iwtd = 1
    for k in 1:nzg
        if wtd < slz[k]
            iwtd = k
            break
        end
    end
    kwtd = iwtd - 1
    
    # 调整地下水位层厚度
    if kwtd >= 1 && kwtd < nzg
        dz2[kwtd] = slz[iwtd] - wtd
    end
    
    # 计算根系最低层
    kroot = 0
    for k in 1:nzg
        if inactivedays[k] <= maxinactivedays
            kroot = k - 1
            break
        end
    end
    
    # 计算各层的水分提取便利性
    for k in max(kwtd, kroot, 1):nzg
        if inactivedays[k] <= maxinactivedays
            rootmask[k] = 1
        end
        
        vctr4_k = 0.5 * (slz[k] + slz[k+1])
        
        # 获取土壤参数
        soil_params = get_soil_params(soiltxt)
        
        # 计算饱和含水量和饱和基质势（考虑深度变化）
        smoisat = soil_params.slmsts * max(min(exp((vctr4_k + 1.5) / fdepth), 1.0), 0.1)
        psisat = soil_params.slpots * min(max(exp(-(vctr4_k + 1.5) / fdepth), 1.0), 10.0)
        
        # 计算基质势
        pot = psisat * (smoisat / smoi[k])^soil_params.slbs
        
        # 考虑冰冻因子
        soilfactor = icefac[k] == 0 ? 1.0 : 0.0
        
        # 计算水分提取便利性
        easy[k] = max(-(POTLEAF - pot) * soilfactor / (hveg - vctr4_k), 0.0)
    end
    
    # 初始化输出变量
    dsmoi = zeros(Float64, nzg)
    dsmoideep = 0.0
    watdef = 0.0
    
    # 消除小的根系活性
    maxeasy = maximum(easy[rootmask .== 1])
    
    # 去除小于最大值0.1%的活性
    for k in 1:nzg
        if easy[k] < 0.001 * maxeasy
            easy[k] = 0.0
        end
    end
    
    # 对于非活跃层，进一步限制根系活性
    for k in max(kroot, 1):nzg
        if inactivedays[k] > maxinactivedays && easy[k] < maxeasy
            easy[k] = 0.0
        end
    end
    
    # 计算总的便利性
    toteasy = sum(easy .* dz2)
    
    if toteasy == 0.0
        rootactivity = zeros(Float64, nzg)
    else
        rootactivity = min.(max.((easy .* dz2) ./ toteasy, 0.0), 1.0)
    end
    
    # 更新非活跃天数
    for k in 1:nzg
        if easy[k] == 0.0
            inactivedays[k] += 1
        else
            inactivedays[k] = 0
        end
    end
    
    # 限制最大非活跃天数
    inactivedays .= min.(inactivedays, maxinactivedays + 1)
    
    # 计算根区土壤湿度和田间持水量
    rootsmoi = 0.0
    rootfc = 0.0
    maxwat = zeros(Float64, nzg)
    
    for k in max(kwtd, 1):nzg
        vctr4_k = 0.5 * (slz[k] + slz[k+1])
        soil_params = get_soil_params(soiltxt)
        
        smoisat = soil_params.slmsts * max(min(exp((vctr4_k + 1.5) / fdepth), 1.0), 0.1)
        psisat = soil_params.slpots * min(max(exp(-(vctr4_k + 1.5) / fdepth), 1.0), 10.0)
        
        smoimin = smoisat * (psisat / POTWILT)^(1.0 / soil_params.slbs)
        smoifc = smoisat * (psisat / POTFC)^(1.0 / soil_params.slbs)
        
        maxwat[k] = max((smoi[k] - smoimin) * dz[k], 0.0)
        rootsmoi += max(rootactivity[k] * (smoi[k] - smoimin), 0.0)
        rootfc += max(rootactivity[k] * (smoifc - smoimin), 0.0)
    end
    
    # 计算土壤水分胁迫因子
    if rootsmoi <= 0.0
        fswp = 0.0
    elseif rootsmoi / rootfc <= 1.0
        fswp = rootsmoi / rootfc
    else
        fswp = 1.0
    end
    
    # 计算冠层和土壤阻力
    rs_c = fswp == 0.0 ? 5000.0 : min(rs_c_factor / fswp, 5000.0)
    
    soil_params = get_soil_params(soiltxt)
    rs_s = 33.5 + 3.5 * (soil_params.slmsts / smoi[nzg])^2.38
    
    R_c = (δ + γ) * ra_c + γ * rs_c
    R_s_final = R_s + γ * rs_s
    
    C_c = 1.0 / (1.0 + R_a * R_c / (R_s_final * (R_c + R_a)))
    C_s = 1.0 / (1.0 + R_a * R_s_final / (R_c * (R_s_final + R_a)))
    
    if lai < 0.001
        C_c = 0.0
    end
    
    # 计算蒸腾和土壤蒸发
    pet_c = C_c * petfactor_c / (δ + γ * (1.0 + rs_c / (ra_a + ra_c)))
    pet_c = max(Δt * pet_c / λ, 0.0)
    
    pet_s = C_s * petfactor_s / (δ + γ * (1.0 + rs_s / (ra_a + ra_c)))
    pet_s = max(Δt * pet_s / λ, 0.0)
    
    transpwater = pet_c * 1.0e-3  # 转换为m
    
    if toteasy == 0.0
        watdef = transpwater
        return pet_s, pet_c, watdef, dsmoi, dsmoideep
    end
    
    # 从各层提取水分
    for k in max(kwtd, 1):nzg
        extract = max(rootactivity[k] * transpwater, 0.0)
        
        if extract <= maxwat[k]
            dsmoi[k] = extract
        else
            dsmoi[k] = maxwat[k]
            watdef += (extract - maxwat[k])
        end
    end
    
    dsmoi .= max.(dsmoi, 0.0)
    
    return pet_s, pet_c, watdef, dsmoi, dsmoideep
end

end # module WaterExtraction
