using Test

# 基础功能测试
include("test_watertable.jl")

# 河流路由测试
include("test_river_routing.jl")

# 同位素追踪测试
include("test_isotope_tracing.jl")
