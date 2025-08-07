# 水位模块主接口
module WaterTable

# 常数定义
const π₄ = 3.1415927 * 4.0
const π2r = 0.0174532925199  # π/180 度转弧度
const g0 = 9.81  # 重力加速度 (m/s²)


# 导入子模块
include("helper.jl")
include("WaterTableCalculations.jl")
include("GroundwaterRiverInteraction.jl")
include("IsotopeTracing.jl")

# 导出主要函数
export wtable!, updatewtd!, lateral_flow!,                                      # 水位计算
    gw2river!, rivers_kw_flood!, rivers_dw_flood!, flooding!, moveqrf!,         # 河流-地下水相互作用
    lateral_isotope!, updatedeepwtable!                                         # 同位素追踪

export flowdir

end # module
