# 水位模块主接口
module WaterTable

const π2r = 0.0174532925199  # π/180 度转弧度

using ASAP

# 导入子模块
include("helper.jl")
include("flooding.jl")
include("lateral_flow.jl")
include("GWRiverInteraction.jl")

include("rivers_dw_flood.jl")
include("rivers_kw_flood.jl")

include("IsotopeTracing.jl")

export flowdir
export wtable!, updatewtd!, lateral_flow!,
  gw2river!,
  rivers_kw_flood!, rivers_dw_flood!, flooding!, moveqrf!,         # 河流-地下水相互作用
  lateral_isotope!, updatedeepwtable!              # 同位素追踪

end # module
