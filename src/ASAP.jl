"""
ASAP模型主模块
整合所有子模块，提供完整的水文模型功能
"""
module ASAP

using Reexport

include("helper.jl")
include("SoilParameters.jl")
include("Evapotranspiration.jl")
include("Interception.jl")
include("extraction.jl")
include("SoilInitialization.jl")
include("SoilFluxes.jl")
include("updatewtd_shallow.jl")
include("updatewtd_qlat.jl")
include("RootDepth.jl")

# 导入新的水位模块
include("wtable/WaterTable.jl")

@reexport using .Evapotranspiration
@reexport using .WaterTable


export updatewtd_shallow, updatewtd_qlat
export rootdepth_main

end # module ASAP
