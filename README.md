# ASAP-model
The ASAP model is composed of the following modules:
  
-module_driver.f90  Driver

-module_forcings.f90 Reads the forcings and interpolates to model grid in each processor

-module_rootdepth.f90 Contains the routines for vertical water fluxes and the soil moisture calculation

-module_wtable.f90 Contains the hydrology routines (groundwater lateral flow and rivers/flooding)

-module_initial.f90 Initialization routines for soil moisture

-module_io.f90 Read in and write out routines

-module_nrtype.f90 Some parameter definitions

-module_parallel.f90 Paralellization routines

-interp_lib.f90 Interpolation routines

The input files and parameters in the model code in the example are set up for South America. In
global simulations, the code runs for each continent separately.

Input parameters needed are top soil textures, topography, efolding depth for transmisivity, initial water table
depth, land cover, vegetation height and several river related parameters (length, width, depth, flow direction, river slope, flood stage).

Forcing variables needed are wind, temperature, pressure, relative humidity,
rainfall, snowmelt, downward solar radiation, net radiation at the surface and
soil temperatures. In the example, they come from ERA5 (hourly), excepth for snowmelt,
which is from ERA5-Land. In addition, LAI at 8 day intervals is also read in
from data from Yuan et al (2011).

Links to download these data are provided below:


-ECMWF-Reanalysis-ERA5

https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5

-Leaf Area Index (LAI)

MODIS Mod15A3-based reprocessed LAI (Yuan et al 2011) (http://globalchange.bnu.edu.cn/research/lai)

-Vegetation Types

MODIS land use (http://www2.mmm.ucar.edu/wrf/users/download/get_sources_wps_geog.html)

-Max Leaf Conductance

NASA Land Data Assimilation System (http://ldas.gsfc.nasa.gov/nldas/web/web.veg.table.html)

-Vegetation Height

ORNL DAAC Global Forest Canopy Height (http://webmap.ornl.gov/wcsdown/dataset.jsp?ds_id=10023)

-Soil Texture

FAO Global Harmonized Soil Data Base (http://webarchive.iiasa.ac.at/Research/LUC/External-World-soil-database/HTML/)

-Topography

Merit DEM (Yamazaki et al 2017) (http://hydro.iis.u-tokyo.ac.jp/~yamadai/MERIT_DEM/)

 

 
