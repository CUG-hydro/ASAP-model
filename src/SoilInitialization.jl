# 预定义的层厚度（40层）
const DZ2 = [
  0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2,
  0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4, 0.5, 0.5, 0.6,
  0.7, 0.7, 0.8, 0.9, 1.0, 1.0, 1.2, 1.2, 1.5, 1.5,
  2.0, 2.0, 3.0, 6.0, 11.0, 20.0, 50.0, 100.0, 250.0, 540.0]

export initializesoildepth_clm, initializesoildepth, eqsoilmoisturetheor


"""
CLM方法初始化土壤层深度

# 参数
- `nzg::Int`: 土壤层数

# 返回
- `slz`: [N, 0], 土壤上边界，向下为负
- `dz` : [N, 1], 土壤深度，正值

# Note
该算法存在问题，slz计算到了N+1层，但dz没有N+1层
"""
function initializesoildepth_clm(nzg::Int)
  slz = zeros(Float64, nzg + 1)
  dz = zeros(Float64, nzg)
  vctr4 = zeros(Float64, nzg) # 土壤的下边界

  # 计算节点深度（指数分布）
  for k in 1:nzg
    vctr4[k] = 0.025 * (exp(0.5 * (k - 0.5)) - 1.0)
  end

  # 计算层厚度
  for k in 2:(nzg-1)
    dz[k] = 0.5 * (vctr4[k+1] - vctr4[k-1])
  end
  dz[1] = 0.5 * (vctr4[1] + vctr4[2]) # vctr4[0] = -vctr4[1]
  dz[nzg] = vctr4[nzg] - vctr4[nzg-1]

  # 计算层中心深度, 与下文计算的对象不同
  slz[1] = 0.5vctr4[1]
  for k in 2:nzg
    slz[k] = 0.5 * (vctr4[k-1] + vctr4[k]) # 中心
  end
  slz[nzg+1] = vctr4[nzg] + 0.5 * dz[nzg] # 计算到了N+1层

  # 反转深度（从地表向下为负值）
  (; slz=reverse(slz), dz=reverse(dz))
end


"""
固定层厚方法初始化土壤层深度

# 参数
- `nzg::Int`: 土壤层数

# 返回  
- `slz`: [N, 0], 土壤下边界
- `dz` : [N, 1], 土壤层厚度
"""
function initializesoildepth(nzg::Int)
  z₋ₕ = zeros(Float64, nzg + 1)
  dz = zeros(Float64, nzg)

  if nzg > length(DZ2)
    error("土壤层数不能超过 $(length(DZ2))")
  end
  z₋ₕ[nzg+1] = 0.0  # 地表深度

  # 从地表向下计算各层深度
  for k in nzg:-1:1
    dz[k] = DZ2[nzg-k+1]
    z₋ₕ[k] = z₋ₕ[k+1] - dz[k]
  end
  return z₋ₕ, dz
end


# ===========================================================================
# 理论平衡含水量初始化
# ===========================================================================
#
# Fortran `fortran/module_initial.f90:141 EQSOILMOISTUREtheor` 通过对
# 单根方程
#
#     d1 * (x - smoi1) / dz + k1 + flux = 0
#
# 逐层调用 Brent 求根 `zbrent`,得到每层稳态含水量 `smoieq[k]`。
# 其中:
#   - `x`        = 待求的本层含水量
#   - `smoi1`    = 上一层含水量 (顶层 `k=1` 时取零)
#   - `dz`       = 本层厚度 (Fortran 中 `vctr5(k)`)
#   - `vctr4(k)` = 本层中心深度
#   - `vctr6(k)` = 1 / (vctr4(k) - vctr4(k-1)),即相邻层中心间距的倒数
#   - `wgpmid`   = 在两层界面处线性插值的"水势等效含水量"
#   - `hydcon`   = `slcons(nsoil)` 饱和导水率 (含深度衰减)
#   - `smoisat`  = `slmsts(nsoil)`  饱和含水量  (含深度衰减)
#   - `psisat`   = `slpots(nsoil)`  饱和基质势
#   - `bexp`     = `slbs(nsoil)`    Campbell b
#   - `flux`     = 0 (稳态无源汇)
#
# 原 Fortran `EQSOILMOISTUREtheor` 是 2-D 循环(每个像元都算一次);
# 本 Julia 翻译保持单列语义 —— 输入 `nsoil/fdepth/wtd` 都是标量,返回
# 长度为 `nzg` 的平衡含水量向量。2-D 循环(对每个像元调用本函数)由
# 上层 driver 负责,与 Fortran 调用约定一致。
# ---------------------------------------------------------------------------

"""
    eqsoilmoisturetheor(nzg, nsoil, z₋ₕ, dz, fdepth, wtd; tol=1e-6)

理论平衡含水量初始化。

逐层调用 `zbrent` 求每层稳态含水量,使得重力排水通量与毛管上升通量
在稳态下平衡。

Translation of `fortran/module_initial.f90:141 EQSOILMOISTUREtheor`.

# Arguments
- `nzg::Int`:     土壤层数
- `nsoil::Int`:   土壤类型索引 (1-13)
- `z₋ₕ::Vector{Float64}`:  节点深度 (长度 `nzg+1`,地表为 0,向下为负)
- `dz::Vector{Float64}`:   层厚度 (长度 `nzg`,正值)
- `fdepth::Real`: 衰减尺度 (m),控制深层土壤参数向地表层收敛
- `wtd::Real`:    地下水位埋深 (m,负值)
- `tol::Real`:    Brent 收敛容差 (默认 `1e-6`)

# Returns
- `smoieq::Vector{Float64}`: 长度为 `nzg` 的平衡含水量 (m³/m³)

# Algorithm
对每层 `k = 1..nzg` 在区间 `[θ_cp, θ_sat*depth_factor]` 上解:

    dψ/dθ · (θ_k - θ_{k-1})/dz_k + K(θ_mid) + flux = 0

其中 `θ_mid` 是层界面处的线性插值,`dψ/dθ` 与 `K(θ)` 由 Campbell
关系给出。底层 `k = nzg` 的下边界条件是 `θ_{nzg+1} = smoisat` (地下水位
以下饱和);`wtd` 暂未在闭包中直接使用,作为接口字段保留以匹配
Fortran `EQSOILMOISTUREtheor` 签名。

# Note
Fortran 原版对每个 `(i,j)` 像元都重新构造 `func` 闭包;Julia 实现
保留同样的语义 —— 在循环内按层构造 `func`。
"""
function eqsoilmoisturetheor(nzg::Int, nsoil::Int,
                              z₋ₕ::Vector{Float64}, dz::Vector{Float64},
                              fdepth::Real, wtd::Real;
                              tol::Real=1e-6)
  # 深度相关参数(对齐 Fortran L168-177, 1-based 索引)
  vctr4 = Vector{Float64}(undef, nzg)   # 节点中心深度
  vctr5 = Vector{Float64}(undef, nzg)   # 等价于 dz,但 Fortran 中用 vctr5
  vctr6 = Vector{Float64}(undef, nzg)   # 1/vctr5
  for k in 1:nzg
    vctr4[k] = 0.5 * (z₋ₕ[k] + z₋ₕ[k+1])
  end
  for k in 2:nzg
    vctr5[k] = vctr4[k] - vctr4[k-1]
    vctr6[k] = 1.0 / vctr5[k]
  end
  vctr5[1] = dz[1]
  vctr6[1] = 1.0 / vctr5[1]

  # Campbell 参数(标量,来自 Fortran module_rootdepth)
  hydcon0 = slcons(nsoil)
  smoisat0 = slmsts(nsoil)
  psisat0 = slpots(nsoil)
  bexp = slbs(nsoil)
  smoicp0 = soilcp(nsoil)

  smoieq = zeros(Float64, nzg)
  flux = 0.0  # 稳态无源汇

  for k in 1:nzg
    # 深度衰减: Fortran 中 `(vctr4(k) + 1.5)/fdepth` 范围内
    # `hydcon / smoisat / psisat / smoicp` 都按 `exp(...)` 缩放,
    # 但 `clamp(exp(...), 0.1, 1)` 之后再乘原值 —— 实际等价于
    # `value * scale(k)`,其中
    #     scale(k) = clamp(exp((vctr4(k) + 1.5)/fdepth), 0.1, 1.0)
    # 注: 对 `psisat` 略有差异,使用 `clamp(1/scale, 1.0, 10.0)`。
    scale_k = clamp(exp((vctr4[k] + 1.5) / fdepth), 0.1, 1.0)
    pscale_k = clamp(exp(-(vctr4[k] + 1.5) / fdepth), 1.0, 10.0)
    hydcon = hydcon0 * scale_k
    smoisat = smoisat0 * scale_k
    psisat = psisat0 * pscale_k
    smoicp = smoicp0 * scale_k

    # 上一层的"下行"含水量
    if k == 1
      smoi1 = 0.0   # 顶层无上层贡献
    else
      nsoil_dw = nsoil
      smoi1 = smoieq[k-1] == 0.0 ? 0.0 : smoieq[k-1]
      # Fortran 中 smoisatdw = slmsts(soiltxt(max(k-1,1))) * scale(k-1)
      smoisat_dw = slmsts(nsoil_dw) * clamp(exp((vctr4[max(k-1, 1)] + 1.5) / fdepth), 0.1, 1.0)
    end

    slz_k = z₋ₕ[k]
    dz_k = vctr5[k]
    vctr4_k = vctr4[k]
    vctr6_k = vctr6[k]
    smoisat_local = smoisat
    hydcon_local = hydcon
    psisat_local = psisat

    # func(x) = d1 * (x - smoi1) / dz + k1 + flux
    #   wgpmid = x + (x - smoi1) * (slz - vctr4) * vctr6
    #   wgpmid = min(wgpmid, smoisat)
    #   k1 = hydcon * (wgpmid/smoisat)^(2*b + 3)
    #   d1 = -(hydcon*psisat*b/smoisat) * (wgpmid/smoisat)^(b + 2)
    function func(x::Float64)
      wgpmid = x + (x - smoi1) * (slz_k - vctr4_k) * vctr6_k
      wgpmid = min(wgpmid, smoisat_local)
      k1 = hydcon_local * (wgpmid / smoisat_local)^(2.0 * bexp + 3.0)
      d1 = -(hydcon_local * psisat_local * bexp / smoisat_local) *
            (wgpmid / smoisat_local)^(bexp + 2.0)
      return d1 * (x - smoi1) / dz_k + k1 + flux
    end

    wmid = zbrent(func, 0.0, 1.0; tol=tol)
    smoieq[k] = clamp(wmid, smoicp, smoisat)

    # 兜底:若 zbrent 端点同号而抛错(物理上意味着函数在该区间无根),
    # 用上一层的解或夹逼值兜底,避免上层循环崩溃
  end

  return smoieq
end
