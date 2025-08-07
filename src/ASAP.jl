"""
根系深度主模块
整合所有子模块，提供主要的根系深度计算功能
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

@reexport using .SoilParameters
@reexport using .Evapotranspiration
@reexport using .WaterExtraction
@reexport using .SoilFluxes
@reexport using .WaterTableDynamics


export rootdepth_main


end # module ASAP
