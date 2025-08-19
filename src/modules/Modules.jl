# include("helper.jl")
include("flooding.jl")
include("lateral_flow.jl")
include("gw2river.jl")

include("rivers_dw_flood.jl")
include("rivers_kw_flood.jl")

include("Tracing/IsotopeTracing.jl")

export wtable!, updatewtd!, lateral_flow!,
  gw2river!,
  rivers_kw_flood!, rivers_dw_flood!, flooding!, moveqrf!,         # 河流-地下水相互作用
  lateral_isotope!, updatedeepwtable!              # 同位素追踪
