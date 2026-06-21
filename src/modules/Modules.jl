# include("helper.jl")
include("flooding.jl")
include("lateral_flow.jl")
include("gw2river.jl")

include("rivers_dw_flood.jl")
include("rivers_kw_flood.jl")

include("Tracing/IsotopeTracing.jl")

# 注意：以下函数未在此模块中实现（实现位于 src/backup/，已被新模块替代）：
#   - wtable! / updatewtd! → 由 updatewtd_shallow + lateral_flow! + gw2river! 组合替代
#   - gw2river! / moveqrf! 在 gw2river.jl 中定义并自导出
export lateral_flow!,
  rivers_kw_flood!, rivers_dw_flood!, flooding!,         # 河流-地下水相互作用
  lateral_isotope!, updatedeepwtable!              # 同位素追踪
