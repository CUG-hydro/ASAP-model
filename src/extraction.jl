export extraction


const POTLEAF = -153.0  # 叶片水势 (m)
const POTWILT = -153.0  # 凋萎点水势 (m)  
const POTFC = -3.366    # 田间持水量水势 (m)

"""
从土壤层中提取水分用于蒸腾

# 参数
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
function extraction(
  nzg::Int, z₋ₕ::Vector{Float64}, dz::Vector{Float64},
  Δt::Float64, soiltxt::Int,
  wtd::Float64, θ::Vector{Float64}, θ_wtd::Float64,
  Δ::Float64, γ::Float64, λ::Float64,
  lai::Float64, ra_a::Float64, ra_c::Float64, rs_c_factor::Float64,
  R_a::Float64, R_s::Float64, petfactor_s::Float64, petfactor_c::Float64,
  inactivedays::Vector{Int}, maxinactivedays::Int,
  hhveg::Float64, fdepth::Float64, icefac::Vector{Int8})

  hveg = 2 / 3 * hhveg # 叶片高度
  easy = zeros(Float64, nzg)
  rootmask = zeros(Int, nzg)
  dz2 = copy(dz)

  # 计算地下水位所在层
  jwt = find_jwt(wtd, z₋ₕ)  # 地下水所在的下一层
  kwtd = jwt - 1           # 地下水所在层

  # 调整地下水位层厚度
  if kwtd >= 1 && kwtd < nzg
    dz2[kwtd] = z₋ₕ[jwt] - wtd  # 非饱和层厚度, 
  end

  # 计算根系最低层
  kroot = 0
  for k in 1:nzg
    if inactivedays[k] <= maxinactivedays
      kroot = k - 1 # 为何返回的是k-1，而不是k? 从k应该更合理
      break
    end
  end
  # 现在的代码，根系没有直接从地下水中抽水

  # 计算各层的水分提取便利性
  for k in max(kwtd, kroot, 1):nzg
    inactivedays[k] <= maxinactivedays && (rootmask[k] = 1)

    _z = 0.5 * (z₋ₕ[k] + z₋ₕ[k+1])
    soil = get_soil_params(soiltxt)

    # 计算饱和含水量和饱和基质势（考虑深度变化）
    θ_sat = soil.θ_sat * clamp(exp((_z + 1.5) / fdepth), 0.1, 1.0)
    ψ_sat = soil.ψsat * clamp(exp(-(_z + 1.5) / fdepth), 1.0, 10.0)

    # 计算土壤水势
    ψ = ψ_sat * (θ_sat / θ[k])^soil.b        # 考虑冰冻因子
    soilfactor = icefac[k] == 0 ? 1.0 : 0.0

    # 计算水分提取便利性
    easy[k] = max(-(POTLEAF - ψ) * soilfactor / (hveg - _z), 0.0)
  end

  # 初始化输出变量
  dθ = zeros(Float64, nzg)
  dθ_deep = 0.0
  watdef = 0.0

  # 消除小的根系活性
  maxeasy = maximum(easy[rootmask.==1])

  # 去除小于最大值0.1%的活性
  for k in 1:nzg
    easy[k] < 0.001 * maxeasy && (easy[k] = 0.0)
  end

  # 对于非活跃层，进一步限制根系活性
  for k in max(kroot, 1):nzg
    if inactivedays[k] > maxinactivedays && easy[k] < maxeasy
      easy[k] = 0.0
    end
  end

  # 计算总的便利性
  toteasy = sum(easy .* dz2) # 不考虑地下水的部分？
  if toteasy == 0.0
    rootactivity = zeros(Float64, nzg)
  else
    rootactivity = clamp((easy .* dz2) ./ toteasy, 0.0, 1.0)
  end

  # 更新非活跃天数
  for k in 1:nzg
    inactivedays[k] = easy[k] == 0.0 ? inactivedays[k] + 1 : 0
  end
  inactivedays .= min.(inactivedays, maxinactivedays + 1)   # 限制最大非活跃天数

  # 计算根区土壤湿度和田间持水量
  θ_root = 0.0
  rootfc = 0.0
  maxwat = zeros(Float64, nzg)

  for k in max(kwtd, 1):nzg
    _z = 0.5 * (z₋ₕ[k] + z₋ₕ[k+1])
    soil = get_soil_params(soiltxt)

    θ_sat = soil.θ_sat * max(min(exp((_z + 1.5) / fdepth), 1.0), 0.1)
    ψ_sat = soil.ψsat * min(max(exp(-(_z + 1.5) / fdepth), 1.0), 10.0)

    θ_min = θ_sat * (ψ_sat / POTWILT)^(1.0 / soil.b)
    θ_fc = θ_sat * (ψ_sat / POTFC)^(1.0 / soil.b)

    maxwat[k] = max((θ[k] - θ_min) * dz[k], 0.0)
    θ_root += max(rootactivity[k] * (θ[k] - θ_min), 0.0)
    rootfc += max(rootactivity[k] * (θ_fc - θ_min), 0.0)
  end

  # 计算土壤水分胁迫因子
  fswp = clamp(θ_root / rootfc, 0.0, 1.0)

  # 计算冠层和土壤阻力
  rs_c = fswp == 0.0 ? 5000.0 : min(rs_c_factor / fswp, 5000.0)

  soil = get_soil_params(soiltxt)
  rs_s = 33.5 + 3.5 * (soil.θ_sat / θ[nzg])^2.38

  R_c = (Δ + γ) * ra_c + γ * rs_c
  R_s_final = R_s + γ * rs_s

  C_c = 1.0 / (1.0 + R_a * R_c / (R_s_final * (R_c + R_a)))
  C_s = 1.0 / (1.0 + R_a * R_s_final / (R_c * (R_s_final + R_a)))

  lai < 0.001 && (C_c = 0.0)

  # 计算蒸腾和土壤蒸发
  pet_c = C_c * petfactor_c / (Δ + γ * (1.0 + rs_c / (ra_a + ra_c)))
  pet_c = max(Δt * pet_c / λ, 0.0)

  pet_s = C_s * petfactor_s / (Δ + γ * (1.0 + rs_s / (ra_a + ra_c)))
  pet_s = max(Δt * pet_s / λ, 0.0)

  transpwater = pet_c * 1.0e-3  # 转换为m

  if toteasy == 0.0
    watdef = transpwater
    return pet_s, pet_c, watdef, dθ, dθ_deep
  end

  # 从各层提取水分，每层都会抽取水分
  for k in max(kwtd, 1):nzg
    extract = max(rootactivity[k] * transpwater, 0.0)

    if extract <= maxwat[k]
      dθ[k] = extract
    else
      dθ[k] = maxwat[k]
      watdef += (extract - maxwat[k]) # 没有得到满足的部分
    end
  end
  dθ .= max.(dθ, 0.0)
  return pet_s, pet_c, watdef, dθ, dθ_deep
end
