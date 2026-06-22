"""
土壤参数模块
定义各种土壤类型的参数和常数
"""

export SoilType, get_soil_params, init_soil_param, cal_K,
  slcons, slmsts, slpots, slbs, soilcp,
  KLATFACTOR  # 侧向导水率标定因子（m/s 标量，1-13）

# 常数定义
const NVTYP = 30  # 植被类型数量
const NSTYP = 13  # 土壤类型数量

# 土壤参数结构体
struct SoilType
  θ_sat::Float64      # 饱和含水量, θ_sat
  θ_cp::Float64       # 残差含水量, θ_cp
  θ_wilt::Float64     # 凋萎点含水量, θ_wilt
  ψsat::Float64       # 饱和基质势, ψsat
  Ksat::Float64       # 饱和导水率, Ksat
  K_latfactor::Float64# 侧向流因子, K_latfactor
  b::Float64          # 土壤b参数, b
end


# 土壤参数数据
const θSAT = [0.395, 0.410, 0.435, 0.485, 0.451, 0.420, 0.477, 0.476, 0.426, 0.492, 0.482, 0.863, 0.476]
const SOILCP = [0.050, 0.052, 0.092, 0.170, 0.125, 0.148, 0.195, 0.235, 0.202, 0.257, 0.268, 0.195, 0.235]
const SLBS = [4.05, 4.38, 4.9, 5.3, 5.39, 7.12, 7.75, 8.52, 10.4, 10.4, 11.4, 7.75, 8.52]
const KSAT = [0.000176, 0.0001563, 0.00003467, 0.0000072, 0.00000695, 0.0000063,
  0.0000017, 0.00000245, 0.000002167, 0.000001033, 0.000001283, 0.0000080, 0.000005787]
const ΨSAT = [-0.121, -0.090, -0.218, -0.786, -0.478, -0.299, -0.356, -0.630, -0.153, -0.490, -0.405, -0.356, -0.630]
"""
侧向流标定因子 K_latfactor（无量纲，索引 1-13 对应 13 种 USDA 土壤类型）。
作为 `lateral_flow!` 的标定系数，乘到 Ksat 上得到侧向流有效导水率。

参考：Cosby et al. (1984) 与 fortran/module_rootdepth.f90::SLKLF。
"""
const KLATFACTOR = [2.0, 3.0, 4.0, 10.0, 12.0, 14.0, 20.0, 24.0, 28.0, 40.0, 48.0, 48.0, 48.0]


"""
获取指定土壤类型的参数

# 参数
- `soil_type::Int`: 土壤类型索引 (1-13)

# 返回
- `SoilType`: 土壤参数结构体
"""
function get_soil_params(soil_type::Int)
  if soil_type < 1 || soil_type > NSTYP
    error("土壤类型索引必须在 1-$NSTYP 范围内")
  end

  return SoilType(
    θSAT[soil_type],       # θ_sat
    SOILCP[soil_type],     # θ_cp
    0.0,                   # θ_wilt 将在初始化时计算
    ΨSAT[soil_type],       # ψsat
    KSAT[soil_type],       # Ksat
    KLATFACTOR[soil_type], # K_latfactor
    SLBS[soil_type],       # b
  )
end

"""
初始化土壤参数，计算凋萎点与田间持水量

# 参数
- `nzg::Int`: 土壤层数

# 返回
- `fieldcp::Matrix{Float64}`: 田间持水量数组 (nzg × nstyp)
- `θ_wilt::Vector{Float64}`: 各土壤类型凋萎点含水量

# 说明
田间持水量采用 Campbell 公式：θ_fc = θ_sat · (ψ_sat / ψ_fc)^(1/b)
其中 ψ_fc = -3.366 m（即 POTFC，来自 extraction.jl）。
该公式与 Fortran `module_rootdepth.f90::INIT_SOIL_PARAM` 中 fieldcp 计算一致。
"""
function init_soil_param(nzg::Int)
  POTWILT_LOCAL = -153.0  # 凋萎点基质势
  POTFC_LOCAL = -3.366    # 田间持水量基质势

  fieldcp = zeros(Float64, nzg, NSTYP)
  θ_wilt = zeros(Float64, NSTYP)

  # 计算各土壤类型的凋萎点与田间持水量
  for nsoil in 1:NSTYP
    θ_wilt[nsoil] = θSAT[nsoil] * (ΨSAT[nsoil] / POTWILT_LOCAL)^(1.0 / SLBS[nsoil])
    # 田间持水量：每层取相同值（不考虑深度衰减，深度因子在调用方处理）
    for k in 1:nzg
      fieldcp[k, nsoil] = θSAT[nsoil] * (ΨSAT[nsoil] / POTFC_LOCAL)^(1.0 / SLBS[nsoil])
    end
  end
  return fieldcp, θ_wilt
end


"""
计算土壤导水率

# 参数
- `smoi::Float64`: 土壤含水量
- `θ_sat::Float64`: 饱和含水量
- `Ksat::Float64`: 饱和导水率
- `b::Float64`: 土壤b参数

# 返回
- `Float64`: 导水率
"""
function cal_K(θ::Float64, θ_sat::Float64, Ksat::Float64, b::Float64)
  return Ksat * (θ / θ_sat)^(2.0 * b + 3.0)
end


# ---------------------------------------------------------------------------
# Fortran `module_rootdepth` accessor aliases
# ---------------------------------------------------------------------------
# The original Fortran code uses the following 1-based scalar accessors
# (see `fortran/module_rootdepth.f90`):
#   slcons(n)  -> Ksat
#   slmsts(n)  -> θ_sat
#   slpots(n)  -> ψ_sat
#   slbs(n)    -> Campbell b parameter
#   soilcp(n)  -> θ_cp (residual moisture)
# These thin wrappers keep the Julia ↔ Fortran port one-to-one, and are
# convenient to call inside hot inner loops (e.g. `eqsoilmoisturetheor`).
# ---------------------------------------------------------------------------

"Fortran `slcons(n)`: 饱和导水率 Ksat (m/s)。索引 1-13。"
slcons(n::Int) = KSAT[n]

"Fortran `slmsts(n)`: 饱和含水量 θ_sat (m³/m³)。索引 1-13。"
slmsts(n::Int) = θSAT[n]

"Fortran `slpots(n)`: 饱和基质势 ψ_sat (m)。索引 1-13。"
slpots(n::Int) = ΨSAT[n]

"Fortran `slbs(n)`: Campbell b 参数。索引 1-13。"
slbs(n::Int) = SLBS[n]

"Fortran `soilcp(n)`: 残余含水量 θ_cp (m³/m³)。索引 1-13。"
soilcp(n::Int) = SOILCP[n]
