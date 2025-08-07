"""
土壤初始化模块
初始化土壤层深度和厚度
"""
module SoilInitialization

export initializesoildepthclm, initializesoildepth

"""
CLM方法初始化土壤层深度

# 参数
- `nzg::Int`: 土壤层数

# 返回
- `Tuple{Vector{Float64}, Vector{Float64}}`: (土壤层深度slz, 土壤层厚度dz)
"""
function initializesoildepthclm(nzg::Int)
    slz = zeros(Float64, nzg + 1)
    slz2 = zeros(Float64, nzg + 1)
    dz = zeros(Float64, nzg)
    dz2 = zeros(Float64, nzg)
    vctr4 = zeros(Float64, nzg)
    
    # 计算节点深度（指数分布）
    for k in 1:nzg
        vctr4[k] = 0.025 * (exp(0.5 * (Float64(k) - 0.5)) - 1.0)
    end
    
    # 计算层厚度
    for k in 2:(nzg-1)
        dz2[k] = 0.5 * (vctr4[k+1] - vctr4[k-1])
    end
    dz2[1] = 0.5 * (vctr4[1] + vctr4[2])
    dz2[nzg] = vctr4[nzg] - vctr4[nzg-1]
    
    # 计算层中心深度
    for k in 1:nzg
        slz2[k] = 0.5 * (vctr4[k] + vctr4[k+1])
    end
    slz2[nzg] = vctr4[nzg] + 0.5 * dz2[nzg]
    
    # 反转深度（从地表向下为负值）
    for k in 1:nzg
        kk = nzg - k + 1
        slz[k] = -slz2[kk]
        dz[k] = dz2[kk]
    end
    
    slz[nzg+1] = 0.0  # 地表深度
    
    return slz, dz
end

"""
固定层厚方法初始化土壤层深度

# 参数
- `nzg::Int`: 土壤层数

# 返回  
- `Tuple{Vector{Float64}, Vector{Float64}}`: (土壤层深度slz, 土壤层厚度dz)
"""
function initializesoildepth(nzg::Int)
    slz = zeros(Float64, nzg + 1)
    dz = zeros(Float64, nzg)
    
    # 预定义的层厚度（40层）
    const DZ2 = [0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 
                  0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4, 0.5, 0.5, 0.6, 
                  0.7, 0.7, 0.8, 0.9, 1.0, 1.0, 1.2, 1.2, 1.5, 1.5, 
                  2.0, 2.0, 3.0, 6.0, 11.0, 20.0, 50.0, 100.0, 250.0, 540.0]
    
    if nzg > length(DZ2)
        error("土壤层数不能超过 $(length(DZ2))")
    end
    
    slz[nzg+1] = 0.0  # 地表深度
    
    # 从地表向下计算各层深度
    for k in nzg:-1:1
        dz[k] = DZ2[nzg - k + 1]
        slz[k] = slz[k+1] - dz[k]
    end
    
    return slz, dz
end

end # module SoilInitialization
