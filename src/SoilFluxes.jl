export soilfluxes, tridag
export find_jwt

function find_jwt(wtd, z₋ₕ)
  nzg = length(z₋ₕ) - 1
  jwt = 1
  for k in 1:nzg
    if wtd < z₋ₕ[k]
      jwt = k # 地下水水位所在的上一层
      break
    end
  end
  return jwt
end

"""
计算土壤中的水分通量和运动

# 参数
- `i::Int`, `j::Int`: 网格坐标
- `nzg::Int`: 土壤层数
- `freedrain::Int`: 自由排水标志
- `dtll::Float64`: 时间步长 (s)
- `slz::Vector{Float64}`: 土壤层深度 (m) [nzg+1]
- `dz::Vector{Float64}`: 土壤层厚度 (m) [nzg]
- `soiltxt::Int`: 土壤类型
- `smoiwtd::Float64`: 地下水含水量
- `transp::Vector{Float64}`: 蒸腾提取量 (m) [nzg]
- `transpdeep::Float64`: 深层蒸腾 (m)
- `smoi::Vector{Float64}`: 土壤含水量 [nzg]
- `wtd::Float64`: 地下水位深度 (m)
- `precip::Float64`: 降水 (mm)
- `pet_s::Float64`: 土壤蒸发潜力 (mm)
- `fdepth::Float64`: 根系深度因子 (m)
- `qlat::Float64`: 侧向流 (m)
- `qrf::Float64`: 河流反馈流 (m)
- `flood::Float64`: 洪水深度 (m)
- `icefactor::Vector{Int8}`: 冰冻因子 [nzg]
- `smoieq::Vector{Float64}`: 平衡含水量 [nzg]
- `o18::Vector{Float64}`: 氧18同位素 [nzg]
- `precipo18::Float64`: 降水氧18
- `tempsfc::Float64`: 地表温度 (K)
- `qlato18::Float64`: 侧向流氧18
- `transpo18::Float64`: 蒸腾氧18

# 返回
- `Tuple`: (et_s, runoff, rech, flux, qrfcorrect, transpo18, updated_smoi, updated_o18)
"""
function soilfluxes(
  nzg::Int, dt::Float64,
  z₋ₕ::Vector{Float64}, Δz::Vector{Float64}, soiltxt::Int,
  θ_wtd::Float64, transp::Vector{Float64}, transpdeep::Float64,
  θ::Vector{Float64}, wtd::Float64, precip::Float64, pet_s::Float64,
  fdepth::Float64, qlat::Float64, qrf::Float64, flood::Float64,
  icefactor::Vector{Int8}; freedrain::Bool=true)
  # θ_eq::Vector{Float64}
  # o18::Vector{Float64}, precipo18::Float64, tempsfc::Float64, qlato18::Float64, transpo18::Float64

  # 计算辅助变量
  z = 0.5 .* (z₋ₕ[1:nzg] .+ z₋ₕ[2:nzg+1]) # 层中心深度, z 
  Δz₊ₕ = zeros(Float64, nzg)
  for k in 2:nzg
    Δz₊ₕ[k] = z[k] - z[k-1] # 两层之间的中心位置
  end

  # 初始化变量
  K_mid = zeros(Float64, nzg)
  D_mid = zeros(Float64, nzg)
  Q = zeros(Float64, nzg + 1)
  gravflux = zeros(Float64, nzg + 1)
  capflux = zeros(Float64, nzg + 1)
  flux = zeros(Float64, nzg + 1)
  qlatflux = zeros(Float64, nzg + 2)

  rech = 0.0
  runoff = 0.0

  # 复制初始值
  θ_old = copy(θ)
  qgw = qlat - qrf

  # 顶部边界条件：入渗 + 土壤蒸发
  soil = get_soil_params(soiltxt)
  θ_cp = soil.θ_cp


  pet_s_actual = θ[nzg] <= θ_cp ? 0.0 : pet_s
  Q[nzg+1] = (-precip + pet_s_actual) * 1.0e-3 - flood # [m], 

  Imax = soil.Ksat * dt
  # 检查入渗能力
  if -Q[nzg+1] > Imax
    runoff = -Q[nzg+1] - Imax
    Q[nzg+1] = -Imax
  end

  # 确定地下水位位置
  if !freedrain
    jwt = 1
    for k in 1:nzg
      if wtd < z₋ₕ[k]
        jwt = k # 地下水水位所在的上一层
        break
      end
    end
  else
    jwt = 0
  end

  # 计算中间层的导水率和扩散率
  for k in max(jwt - 1, 2):nzg
    θmid = θ[k] + (θ[k] - θ[k-1]) * (z₋ₕ[k] - z[k]) / Δz₊ₕ[k] # 边界处θ

    # 获取土壤参数
    soil = get_soil_params(soiltxt)
    (; Ksat, θ_sat, ψsat, b) = soil

    _Ksat = soil.Ksat * clamp(exp((z₋ₕ[k] + 1.5) / fdepth), 0.1, 1.0)
    _θsat = soil.θ_sat * clamp(exp((z₋ₕ[k] + 1.5) / fdepth), 0.1, 1.0)
    _ψsat = soil.ψsat * clamp(exp(-(z₋ₕ[k] + 1.5) / fdepth), 1.0, 10.0)

    θmid = min(θmid, _θsat)

    # 考虑冰冻因子
    f_ice = icefactor[k] == 0 ? 1.0 : 0.0

    K_mid[k] = f_ice * _Ksat * (θmid / _θsat)^(2.0 * b + 3.0)
    D_mid[k] = -f_ice * (_Ksat * _ψsat * b / _θsat) *
               (θmid / _θsat)^(b + 2.0) # D = K(θ) pdv(ϕ, θ), 经过严格推导
  end

  # 计算三对角矩阵元素
  aa = zeros(Float64, nzg)
  bb = zeros(Float64, nzg)
  cc = zeros(Float64, nzg)
  rr = zeros(Float64, nzg)

  for k in max(jwt, 3):(nzg-1)
    aa[k] = D_mid[k] / Δz₊ₕ[k]
    cc[k] = D_mid[k+1] / Δz₊ₕ[k+1]
    bb[k] = -(aa[k] + cc[k] + Δz[k] / dt)
    rr[k] = -θ[k] * Δz[k] / dt - K_mid[k+1] + K_mid[k] + transp[k] / dt
  end

  # 顶部边界条件
  if jwt - 1 == nzg # 若地下水没过地表, Q_k-1近似重力排水
    aa[nzg] = 0.0
    cc[nzg] = 0.0
    bb[nzg] = -Δz[nzg] / dt
    rr[nzg] = Q[nzg+1] / dt - θ[nzg] * Δz[nzg] / dt + transp[nzg] / dt +
              min(K_mid[nzg] + D_mid[nzg] / Δz₊ₕ[nzg] * (θ[nzg] - θ[nzg-1]), 0.0)
  else
    aa[nzg] = D_mid[nzg] / Δz₊ₕ[nzg]
    cc[nzg] = 0.0
    bb[nzg] = -aa[nzg] - Δz[nzg] / dt
    rr[nzg] = Q[nzg+1] / dt - θ[nzg] * Δz[nzg] / dt + K_mid[nzg] + transp[nzg] / dt
  end

  # 底部边界条件
  if !freedrain # 非自由排水
    if jwt <= 2
      aa[1] = 0.0
      cc[1] = D_mid[2] / Δz₊ₕ[2]
      bb[1] = -(cc[1] + Δz[1] / dt)
      rr[1] = -θ[1] * Δz[1] / dt - K_mid[2] + transp[1] / dt

      k = 2
      aa[k] = D_mid[k] / Δz₊ₕ[k]
      if k < nzg
        cc[k] = D_mid[k+1] / Δz₊ₕ[k+1]
        bb[k] = -(aa[k] + cc[k] + Δz[k] / dt)
        rr[k] = -θ[k] * Δz[k] / dt - K_mid[k+1] + K_mid[k] + transp[k] / dt
      else
        cc[k] = 0.0
        bb[k] = -(aa[k] + Δz[k] / dt)
        rr[k] = -θ[k] * Δz[k] / dt + K_mid[k] + transp[k] / dt
      end
    else
      # 其他层保持不变
      for k in 1:(jwt-3)
        aa[k] = 0.0
        cc[k] = 0.0
        bb[k] = 1.0
        rr[k] = θ[k]
      end

      k = jwt - 1  # 地下水位层
      aa[k] = 0.0
      if k < nzg
        cc[k] = D_mid[k+1] / Δz₊ₕ[k+1]
        bb[k] = -(cc[k] + Δz[k] / dt)
        rr[k] = -θ[k] * Δz[k] / dt - K_mid[k+1] + transp[k] / dt +
                min(K_mid[k] + D_mid[k] / Δz₊ₕ[k] * (θ[k] - θ[k-1]), 0.0)
      else
        cc[k] = 0.0
        bb[k] = -Δz[k] / dt
        rr[k] = -θ[k] * Δz[k] / dt + transp[k] / dt +
                min(K_mid[k] + D_mid[k] / Δz₊ₕ[k] * (θ[k] - θ[k-1]), 0.0)
      end

      k = jwt - 2
      aa[k] = 0.0
      cc[k] = 0.0
      bb[k] = -Δz[k] / dt
      rr[k] = -θ[k] * Δz[k] / dt +
              max(-K_mid[k+1] - D_mid[k+1] / Δz₊ₕ[k+1] * (θ[k+1] - θ[k]), 0.0)
    end
  else
    # 重力排水边界
    soil = get_soil_params(soiltxt)
    K = soil.Ksat * max(min(exp((z₋ₕ[1] + 1.5) / fdepth), 1.0), 0.1)
    smoisat = soil.θ_sat * max(min(exp((z₋ₕ[1] + 1.5) / fdepth), 1.0), 0.1)

    K_mid[1] = K * (θ[1] / smoisat)^(2.0 * soil.b + 3.0)
    aa[1] = 0.0
    cc[1] = D_mid[2] / Δz₊ₕ[2]
    bb[1] = -(cc[1] + Δz[1] / dt)
    rr[1] = -θ[1] * Δz[1] / dt - K_mid[2] + K_mid[1] + transp[1] / dt
  end

  # 求解三对角系统
  tridag!(aa, bb, cc, rr, θ)

  # 计算通量
  for k in max(jwt, 3):nzg
    gravflux[k] = -K_mid[k] * dt
    capflux[k] = -aa[k] * (θ[k] - θ[k-1]) * dt
    Q[k] = capflux[k] + gravflux[k]
  end

  if jwt <= 2
    capflux[1] = 0.0
    gravflux[1] = 0.0
    Q[1] = 0.0
    gravflux[2] = -K_mid[2] * dt
    capflux[2] = -aa[2] * (θ[2] - θ[1]) * dt
    Q[2] = capflux[2] + gravflux[2]
  else
    for k in 1:(jwt-2)
      capflux[k] = 0.0
      gravflux[k] = 0.0
      Q[k] = 0.0
    end

    k = jwt - 1
    gravflux[k] = -K_mid[k] * dt
    capflux[k] = -D_mid[k] / Δz₊ₕ[k] * (θ_old[k] - θ_old[k-1]) * dt

    if capflux[k] > -gravflux[k]
      Q[k] = capflux[k] + gravflux[k]
    elseif capflux[k] > 0.0
      gravflux[k] = -capflux[k]
      Q[k] = 0.0
    else
      capflux[k] = 0.0
      gravflux[k] = 0.0
      Q[k] = 0.0
    end
  end

  Q[1] = freedrain ? -K_mid[1] * dt : 0.0 # 为何乘dt ? 

  # 重新计算土壤含水量
  θ .= θ_old
  for k in 1:nzg
    θ_old[k] = θ_old[k] + (Q[k] - Q[k+1] - transp[k]) / Δz[k]
  end

  # 检查并修正土壤含水量边界
  for k in 1:nzg
    soil = get_soil_params(soiltxt)
    smoisat = soil.θ_sat * max(min(exp((z[k] + 1.5) / fdepth), 1.0), 0.1)

    if θ_old[k] > smoisat
      dsmoi = max((θ_old[k] - smoisat) * Δz[k], 0.0)
      if k < nzg
        θ_old[k+1] = θ_old[k+1] + dsmoi / Δz[k+1]
        Q[k+1] = Q[k+1] + dsmoi
      else
        Q[k+1] = Q[k+1] + dsmoi
        runoff = runoff + dsmoi
      end
      θ_old[k] = smoisat

      if capflux[k+1] < 0.0
        gravflux[k+1] = gravflux[k+1] + capflux[k+1]
        capflux[k+1] = 0.0
      end
      gravflux[k+1] = gravflux[k+1] + dsmoi
      if gravflux[k+1] > 0.0
        capflux[k+1] = capflux[k+1] + gravflux[k+1]
        gravflux[k+1] = 0.0
      end
    end
  end

  # 处理最上层的干燥限制
  k = nzg
  soil = get_soil_params(soiltxt)
  θ_cp = soil.ρb * max(min(exp((z[k] + 1.5) / fdepth), 1.0), 0.1)

  et_s = pet_s
  if θ_old[k] < θ_cp
    dsmoi = max((θ_cp - θ_old[k]) * Δz[k], 0.0)
    if Q[k+1] > dsmoi
      et_s = max(0.0, pet_s - dsmoi * 1.0e3)
      θ_old[k] = θ_cp
      Q[k+1] = Q[k+1] - dsmoi
    else
      et_s = max(0.0, pet_s - max(Q[k+1], 0.0) * 1.0e3)
      Q[k+1] = min(Q[k+1], 0.0)
      θ_old[k] = θ_old[k] + max(Q[k+1], 0.0) / Δz[k]

      dsmoi = max((θ_cp - θ_old[k]) * Δz[k], 0.0)
      θ_old[k-1] = θ_old[k-1] - dsmoi / Δz[k-1]
      Q[k] = Q[k] + dsmoi
      θ_old[k] = θ_cp
    end
  end

  # 继续处理下层
  for k in (nzg-1):-1:1
    soil = get_soil_params(soiltxt)
    θ_cp = soil.ρb * max(min(exp((z[k] + 1.5) / fdepth), 1.0), 0.1)

    if θ_old[k] < θ_cp
      dsmoi = max((θ_cp - θ_old[k]) * Δz[k], 0.0)
      if k > 1
        θ_old[k-1] = θ_old[k-1] - dsmoi / Δz[k-1]
      end
      Q[k] = Q[k] + dsmoi
      θ_old[k] = θ_cp
    end
  end

  qrfcorrect = 0.0
  if Q[1] > 0.0
    qrfcorrect = -min(Q[1], max(qrf, 0.0))
  end

  # 累积通量
  flux .+= Q

  # 处理自由排水
  if freedrain == 1
    rech = Q[1]
    # smoiwtd = smoiwtd - Q[1]  # 这个变量需要在调用函数中处理
  end

  # o18ratio = o18 ./ θ
  # 氧18同位素计算（简化版本）
  # transpo18_out = 0.0
  # for k in 1:nzg
  #     transpo18_out += 0.5 * (o18[k] / θ_old[k] + o18ratio[k]) * transp[k]
  #     if o18[k] < 0.0
  #         println("警告: O18小于零! i=$i, j=$j, k=$k")
  #     end
  #     if o18[k] > θ_old[k]
  #         println("警告: O18大于土壤含水量! i=$i, j=$j, k=$k")
  #     end
  # end
  # o18 .= max.(o18, 0.0)
  θ .= θ_old
  return et_s, runoff, rech, flux, qrfcorrect, θ
end


"""
三对角矩阵求解器（原地修改版本）

# 参数
- `a::Vector{Float64}`: 下对角线元素
- `b::Vector{Float64}`: 主对角线元素  
- `c::Vector{Float64}`: 上对角线元素
- `r::Vector{Float64}`: 右端向量
- `u::Vector{Float64}`: 解向量（输出）
"""
function tridag!(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64},
  r::Vector{Float64}, u::Vector{Float64})
  n = length(b)
  gam = zeros(Float64, n)

  if b[1] == 0.0
    error("tridag: 需要重写方程")
  end

  bet = b[1]
  u[1] = r[1] / bet

  for j in 2:n
    gam[j] = c[j-1] / bet
    bet = b[j] - a[j] * gam[j]
    if bet == 0.0
      error("tridag 失败")
    end
    u[j] = (r[j] - a[j] * u[j-1]) / bet
  end

  for j in (n-1):-1:1
    u[j] = u[j] - gam[j+1] * u[j+1]
  end
end
