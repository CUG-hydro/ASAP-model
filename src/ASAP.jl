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
# include("updatewtd_qlat.jl")
include("RootDepth.jl")

# 导入新的水位模块
include("modules/Modules.jl")

@reexport using .Evapotranspiration


export updatewtd_shallow
export rootdepth_main
export zbrent
export eqsoilmoisturetheor

# NetCDF I/O (read-only subset of fortran/module_io.f90)
include("io/NetCDF.jl")
using .NetCDFIO: read_initial, read_wtdnc
export read_initial, read_wtdnc

# ERA5 forcing readers (read-only subset of fortran/module_forcings.f90)
include("Forcings/ERA5.jl")
using .ERA5Forcings
export read_hourly_forcings, read_accumulated_forcings,
       read_snow, read_snow_hour, read_soil_temps,
       read_topo_era5, read_lai_climatology, compute_qair

end # module ASAP
