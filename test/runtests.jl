using ASAP, Test


# 运行各个模块的测试
include("test_soil_initialization.jl")
include("test_evapotranspiration.jl")
include("test_interception.jl")
include("test_soil_parameters.jl")
include("test_soil_fluxes.jl")
include("test_extraction.jl")
include("test_updatewtd_shallow.jl")

include("wtable/runtests.jl")
