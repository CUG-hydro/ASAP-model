"""
ASAP模型主模块
整合所有子模块，提供完整的水文模型功能
"""
module ASAP

using Reexport

# 导入所有子模块
include("SoilParameters.jl")
include("Evapotranspiration.jl")
include("Interception.jl")
include("WaterExtraction.jl")
include("SoilFluxes.jl")
include("WaterTableDynamics.jl")
include("SoilInitialization.jl")
include("RootDepth.jl")

# 导入新的水位模块
include("wtable/WaterTable.jl")

@reexport using .SoilParameters
@reexport using .Evapotranspiration
@reexport using .WaterExtraction
@reexport using .WaterTable

export rootdepth_main

end # module ASAP
