"""
ERA5 forcing readers (subset of `fortran/module_forcings.f90`).

Translates 6 Fortran subroutines (`READFORCINGS`, `READFORCINGSACC`,
`READFORCINGSSNOW`, `READFORCINGSSOILT`, `READTOPOERA5`, `READLAICLIM`)
into a serial, MPI-free Julia interface, plus a `compute_qair` helper.

## Region portability

Every reader takes a `grid::Grid` keyword that bundles the four pieces of
metadata needed to map between a source data grid and a model grid:

- `grid.origin` — `(xswlon, xswlat)` SW-corner longitude/latitude (deg)
- `grid.res`    — `(dx, dy)` cell size (deg)
- `grid.shape`  — `(nx, ny)` global grid dimensions

Two production grids ship as module constants: `ERA5_GRID` (0.25°,
1440×721) and `ERA5LAND_GRID` (0.1°, 3600×1801). Tests supply a custom
`Grid` matching their mock files; production callers can omit the
keyword and get the same behaviour as the original Fortran port.

The regional crop range is passed as `xrange` / `yrange` keywords
(defaults preserve the Fortran `[1067:1316, 298:587]` window).
"""
module ERA5Forcings

using NCDatasets

export Grid,
       read_hourly_forcings, read_accumulated_forcings,
       read_snow, read_snow_hour, read_soil_temps,
       read_topo_era5, read_lai_climatology,
       compute_qair

# ===========================================================================
# Grid metadata
# ===========================================================================

"""
    Grid(name, origin, res, shape)

Source-data grid descriptor. `origin = (xswlon, xswlat)`,
`res = (dx, dy)` in degrees, `shape = (nx_global, ny_global)`.

Two production grids are pre-defined:
- `ERA5_GRID`     — 0.25° ERA5 (1440×721), origin (266.5, -56.5)
- `ERA5LAND_GRID` — 0.1°  ERA5-Land (3600×1801), origin (266.5, -56.5)
"""
struct Grid
    name::String
    origin::Tuple{Float64,Float64}
    res::Tuple{Float64,Float64}
    shape::Tuple{Int,Int}
end

"""Production ERA5 grid: 0.25°, 1440×721, SW corner at (266.5°E, -56.5°N)."""
const ERA5_GRID = Grid("ERA5", (266.5, -56.5), (0.25, 0.25), (1440, 721))

"""Production ERA5-Land grid: 0.1°, 3600×1801, SW corner at (266.5°E, -56.5°N)."""
const ERA5LAND_GRID = Grid("ERA5-Land", (266.5, -56.5), (0.10, 0.10), (3600, 1801))

# Fortran `varread(1067:1316, 298:587, :)` (1-based inclusive). Defaults
# preserved for backward compatibility with the original port.
const ERA5_REGION_XRANGE = 1067:1316
const ERA5_REGION_YRANGE = 298:587

# ERA5-Land snow: 2666:3285 × 747:1466 (Fortran L628).
const ERA5LAND_SNOW_XRANGE = 2666:3285
const ERA5LAND_SNOW_YRANGE = 747:1466

const HOURS_PER_DAY = 24
const N_SOIL_LEVELS = 4
const T_FREEZE_K    = 273.15

# ===========================================================================
# Path template (matches fortran/module_forcings.f90 L69-L211)
# ===========================================================================

"""
    era5_hourly_path(root, varname, date) -> String

Build the canonical NetCDF path `\$root/\$varname/ERA5_\$varname_\$date.nc`.
"""
era5_hourly_path(root::String, varname::String, date::String) =
    joinpath(root, varname, "ERA5_$(varname)_$(date).nc")

# ===========================================================================
# Internal helpers
# ===========================================================================

"""
    _crop3d(path, varname; xrange, yrange) -> Array{Float64,3}

Read the full 3-D `varname` from `path` and crop to the regional
subdomain `(length(xrange), length(yrange), :)`.
"""
function _crop3d(path::String, varname::String;
                 xrange::AbstractUnitRange, yrange::AbstractUnitRange)
    full = NCDataset(path) do ds
        Array{Float64}(ds[varname][:, :, :])
    end
    return full[xrange, yrange, :]
end

"""
    _grid_index(lats, lons, src_grid) -> (grx, gry)

Map model-grid lat/lon (degrees) to fractional source-grid indices.
Longitude is wrapped to `[0, 360)` if negative. Mirrors Fortran
`gry = ny + (xswlat - lat)/dy` (north-up indexing).
"""
@inline function _grid_index(lats::AbstractMatrix, lons::AbstractMatrix,
                             src_grid::Grid)
    (xswlon, xswlat) = src_grid.origin
    (dx, dy)         = src_grid.res
    ny_src            = src_grid.shape[2]
    nx, ny = size(lats)
    grx = Array{Float64}(undef, nx, ny)
    gry = Array{Float64}(undef, nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        xlon = lons[i, j] < 0.0 ? lons[i, j] + 360.0 : lons[i, j]
        grx[i, j] = (xlon - xswlon) / dx + 1.0
        gry[i, j] = ny_src + (xswlat - lats[i, j]) / dy
    end
    return grx, gry
end

"""
    bilinear(src, src_grid, lats, lons; landmask, flip_lat=false) -> Matrix

Bilinearly interpolate a 2-D source field `src` (indexed as the regional
crop — i.e. `(length(xrange), length(yrange))`) to the model grid
defined by `lats`, `lons`. Ocean cells (`landmask == 0`) are left at zero.

When `flip_lat = true`, the source is reversed along the y-axis before
interpolation (matches Fortran `topoera(:, ydim-j+1) = ...`). This is the
behaviour `read_topo_era5` needs for some legacy ERA5 static files.
"""
function bilinear(src::AbstractMatrix, src_grid::Grid,
                  lats::AbstractMatrix, lons::AbstractMatrix;
                  landmask::AbstractMatrix,
                  flip_lat::Bool=false)
    src_used = flip_lat ? reverse(src, dims=2) : src
    rx, ry   = size(src_used)
    nx, ny   = size(lats)
    out      = zeros(Float64, nx, ny)

    grx, gry = _grid_index(lats, lons, src_grid)

    @inbounds for j in 1:ny, i in 1:nx
        landmask[i, j] == 0 && continue
        gx = grx[i, j]; gy = gry[i, j]
        ii = clamp(floor(Int, gx), 1, rx - 1)
        jj = clamp(floor(Int, gy), 1, ry - 1)
        fx = gx - ii; fy = gy - jj
        out[i, j] = (1.0 - fx) * (1.0 - fy) * src_used[ii,     jj]     +
                    fx       * (1.0 - fy) * src_used[ii + 1, jj]     +
                    (1.0 - fx) * fy       * src_used[ii,     jj + 1] +
                    fx       * fy       * src_used[ii + 1, jj + 1]
    end
    return out
end

# ===========================================================================
# 1. READFORCINGS  →  read_hourly_forcings
# ===========================================================================

"""
    read_hourly_forcings(date::String, root::String;
                         grid::Grid = ERA5_GRID,
                         xrange = ERA5_REGION_XRANGE,
                         yrange = ERA5_REGION_YRANGE) -> NamedTuple

Reads 24 h of ERA5 hourly fields and packs them into a regional
`varpack` of shape `(rx, ry, 24, 11)` where `(rx, ry) = (length(xrange),
length(yrange))`. The 11 slots are documented in the module header.
"""
function read_hourly_forcings(date::String, root::String;
                              grid::Grid = ERA5_GRID,
                              xrange::AbstractUnitRange = ERA5_REGION_XRANGE,
                              yrange::AbstractUnitRange = ERA5_REGION_YRANGE)
    rx, ry = length(xrange), length(yrange)
    path = (v) -> era5_hourly_path(root, v, date)

    wind    = _crop3d(path("WIND"),     "ws10"; xrange=xrange, yrange=yrange)
    tempera = _crop3d(path("TEMP"),     "t2m";  xrange=xrange, yrange=yrange)
    press   = _crop3d(path("SFCPRESS"), "sp";   xrange=xrange, yrange=yrange)
    dtemp   = _crop3d(path("DEWPOINT"), "d2m";  xrange=xrange, yrange=yrange)

    # RH from saturation vapour pressure (Fortran ports the no-`*611.2` form).
    pressesat  = exp.(17.67 .* (tempera .- 273.15) ./
                      (243.5 .+ (tempera .- 273.15)))
    pressesattd = exp.(17.67 .* (dtemp .- 273.15) ./
                       (243.5 .+ (dtemp .- 273.15)))
    rh = clamp.(pressesattd ./ pressesat, 0.0, 1.0)

    # Precipitation: tp - sf (m → mm, clamped ≥ 0).
    tp = _crop3d(path("PRECIP"), "tp"; xrange=xrange, yrange=yrange)
    sf = _crop3d(path("PRECIP"), "sf"; xrange=xrange, yrange=yrange)
    rain = max.(tp .- sf, 0.0) .* 1000.0

    swdown  = _crop3d(path("RADIATION"), "ssrd"; xrange=xrange, yrange=yrange)
    rnetera = _crop3d(path("RADIATION"), "sr";   xrange=xrange, yrange=yrange)

    # 4 soil-temperature levels in one file.
    st_arrs = ntuple(Val(N_SOIL_LEVELS)) do k
        _crop3d(path("SOILTEMP"), "stl$k"; xrange=xrange, yrange=yrange)
    end

    varpack = Array{Float64}(undef, rx, ry, HOURS_PER_DAY, 11)
    varpack[:, :, :,  1] = wind
    varpack[:, :, :,  2] = tempera
    varpack[:, :, :,  3] = press
    varpack[:, :, :,  4] = rh
    varpack[:, :, :,  5] = rain
    varpack[:, :, :,  6] = swdown
    varpack[:, :, :,  7] = rnetera
    varpack[:, :, :,  8] = st_arrs[1]
    varpack[:, :, :,  9] = st_arrs[2]
    varpack[:, :, :, 10] = st_arrs[3]
    varpack[:, :, :, 11] = st_arrs[4]

    return (varpack=varpack, wind=wind, tempera=tempera, press=press,
            rh=rh, rain=rain, swdown=swdown, rnetera=rnetera,
            st1=st_arrs[1], st2=st_arrs[2], st3=st_arrs[3], st4=st_arrs[4])
end

# ===========================================================================
# 2. READFORCINGSACC  →  read_accumulated_forcings
# ===========================================================================

"""
    read_accumulated_forcings(hour::Int, varpack::Array{Float64,4};
                              lats, lons, landmask,
                              grid::Grid = ERA5_GRID) -> NamedTuple

Bilinearly interpolate the 3 accumulated fields (rain, shortwave, net
radiation) from `varpack` to the model grid. Radiation fluxes are
converted from J/m² to W/m² (÷ 3600).
"""
function read_accumulated_forcings(hour::Int, varpack::Array{Float64,4};
                                   lats::AbstractMatrix, lons::AbstractMatrix,
                                   landmask::AbstractMatrix,
                                   grid::Grid = ERA5_GRID)
    nx, ny = size(lats)
    rain_h    = view(varpack, :, :, hour, 5)
    swdown_h  = view(varpack, :, :, hour, 6)
    rnetera_h = view(varpack, :, :, hour, 7)

    pcpgl  = bilinear(rain_h,    grid, lats, lons; landmask=landmask)
    rshort = bilinear(swdown_h,  grid, lats, lons; landmask=landmask) ./ 3600.0
    netrad = bilinear(rnetera_h, grid, lats, lons; landmask=landmask) ./ 3600.0
    return (pcpgl=pcpgl, rshort=rshort, netrad=netrad)
end

# ===========================================================================
# 3. READFORCINGSSNOW  →  read_snow_hour, read_snow
# ===========================================================================

"""
    read_snow_hour(date::String, hour::Int, root::String;
                   lats, lons, landmask,
                   grid::Grid = ERA5LAND_GRID,
                   xrange = ERA5LAND_SNOW_XRANGE,
                   yrange = ERA5LAND_SNOW_YRANGE) -> Matrix{Float64}

Read a single hour of ERA5-Land `smlt` (snowmelt), convert m → mm,
clamp ≥ 0, and bilinearly interpolate to the model grid.
"""
function read_snow_hour(date::String, hour::Int, root::String;
                        lats::AbstractMatrix, lons::AbstractMatrix,
                        landmask::AbstractMatrix,
                        grid::Grid = ERA5LAND_GRID,
                        xrange::AbstractUnitRange = ERA5LAND_SNOW_XRANGE,
                        yrange::AbstractUnitRange = ERA5LAND_SNOW_YRANGE)
    smlt = _crop3d(joinpath(root, "SNOW", "ERA5L_snow_$(date).nc"),
                   "smlt"; xrange=xrange, yrange=yrange)
    smlt = max.(smlt, 0.0) .* 1000.0
    return bilinear(smlt[:, :, hour], grid, lats, lons; landmask=landmask)
end

"""
    read_snow(date::String, root::String; kwargs...) -> Matrix{Float64}

Convenience alias for `read_snow_hour(date, 1, root; kwargs...)`.
"""
read_snow(date::String, root::String; kwargs...) =
    read_snow_hour(date, 1, root; kwargs...)

# ===========================================================================
# 4. READFORCINGSSOILT  →  read_soil_temps
# ===========================================================================

"""
    read_soil_temps(hour::Int, varpack::Array{Float64,4};
                    nzg::Int, lats, lons, landmask,
                    topo, topoera,
                    grid::Grid = ERA5_GRID) -> Array{Int8,3}

Build the Fortran `icefactor` array from the 4 ERA5 soil-temperature
levels in `varpack`. Each level is bilinearly interpolated, then
lapse-rate-corrected (`topofactor = -0.0065*(topo - topoera)`), and
flagged frozen (`1`) when ≤ 273.15 K. The third axis has length 15
(corresponds to Fortran `icefactor(:, :, nzg-14:nzg)`).
"""
function read_soil_temps(hour::Int, varpack::Array{Float64,4};
                         nzg::Int,
                         lats::AbstractMatrix, lons::AbstractMatrix,
                         landmask::AbstractMatrix,
                         topo::AbstractMatrix, topoera::AbstractMatrix,
                         grid::Grid = ERA5_GRID)
    nx, ny = size(lats)
    nlay   = 15  # Fortran: indices nzg-14 .. nzg
    icefactor = zeros(Int8, nx, ny, nlay)

    grx, gry = _grid_index(lats, lons, grid)
    rx, ry   = size(varpack, 1), size(varpack, 2)

    @inbounds for j in 1:ny, i in 1:nx
        landmask[i, j] == 0 && continue
        gx = grx[i, j]; gy = gry[i, j]
        ii = clamp(floor(Int, gx), 1, rx - 1)
        jj = clamp(floor(Int, gy), 1, ry - 1)
        fx = gx - ii; fy = gy - jj
        w11 = (1.0 - fx) * (1.0 - fy)
        w21 = fx * (1.0 - fy)
        w12 = (1.0 - fx) * fy
        w22 = fx * fy

        topofactor = -0.0065 * (topo[i, j] - topoera[i, j])

        for nlev in 1:N_SOIL_LEVELS
            slice = view(varpack, :, :, hour, nlev + 7)  # stl1..stl4
            tempint = w11 * slice[ii,     jj]     +
                      w21 * slice[ii + 1, jj]     +
                      w12 * slice[ii,     jj + 1] +
                      w22 * slice[ii + 1, jj + 1]
            soilt = tempint + topofactor
            soilt <= T_FREEZE_K || continue
            if nlev == 1
                icefactor[i, j, nlay] = 1
            elseif nlev == 2
                icefactor[i, j, nlay - 2 : nlay - 1] .= 1
            elseif nlev == 3
                icefactor[i, j, nlay - 7 : nlay - 3] .= 1
            else  # nlev == 4
                icefactor[i, j, nlay - 14 : nlay - 8] .= 1
            end
        end
    end
    return icefactor
end

# ===========================================================================
# 5. READTOPOERA5  →  read_topo_era5
# ===========================================================================

"""
    read_topo_era5(path::String;
                   lats, lons, landmask,
                   grid::Grid = ERA5_GRID,
                   scale_factor::Real = 1.0,
                   add_offset::Real   = 0.0) -> Matrix{Float64}

Reads the static ERA5 `z` (geopotential, m²/s²), converts to geometric
height (m, ÷ 9.81), reverses along the y-axis to match Fortran's
north-up convention, and bilinearly interpolates to the model grid.
"""
function read_topo_era5(path::String;
                        lats::AbstractMatrix, lons::AbstractMatrix,
                        landmask::AbstractMatrix,
                        grid::Grid = ERA5_GRID,
                        scale_factor::Real = 1.0,
                        add_offset::Real   = 0.0)
    nx_g, ny_g = grid.shape
    zraw = NCDataset(path) do ds
        Array{Float64}(ds["z"][:, :])
    end
    # Fortran L839-L843 north-up reversal; preserve for legacy compatibility.
    topoera = zeros(Float64, nx_g, ny_g)
    @inbounds for j in 1:ny_g
        topoera[:, ny_g - j + 1] .= zraw[:, j] .* scale_factor .+ add_offset
    end
    topoera ./= 9.81
    return bilinear(topoera, grid, lats, lons;
                    landmask=landmask, flip_lat=false)
end

# ===========================================================================
# 6. READLAICLIM  →  read_lai_climatology
# ===========================================================================

"""
    read_lai_climatology(path::String, month::Int;
                         target_shape::Tuple{Int,Int}) -> Matrix{Float64}

Reads the LAI monthly-climatology variable (3-D `LAI(x, y, month)`) and
returns the `month`-th slice. The source file is assumed to already be
on the model grid (matches the Fortran `READLAICLIM` crop convention).

# Errors
Throws `DimensionMismatch` if the source `LAI` slice does not match
`target_shape`.
"""
function read_lai_climatology(path::String, month::Int;
                              target_shape::Tuple{Int,Int})
    n2, n3 = target_shape
    NCDataset(path) do ds
        full = Array{Float64}(ds["LAI"][:, :, month])
        size(full, 1) == n2 && size(full, 2) == n3 ||
            throw(DimensionMismatch(
                "LAI slice size $(size(full)) ≠ target_shape=$target_shape"))
        return full
    end
end

# ===========================================================================
# 7. qair from dewpoint and surface pressure
# ===========================================================================

"""
    compute_qair(d2m_K, sp_Pa) -> Array{Float64}

Specific humidity (kg/kg) from 2 m dew-point (`d2m`, K) and surface
pressure (`sp`, Pa) via the Tetens formula:
`e = 6.112·exp(17.67·(Td − 273.15)/(Td − 29.65))·100 Pa`,
`q = 0.622·e/(sp − 0.378·e)`.
"""
compute_qair(d2m_K::AbstractArray, sp_Pa::AbstractArray) =
    let e = 6.112 .* exp.(17.67 .* (d2m_K .- 273.15) ./ (d2m_K .- 29.65)) .* 100.0
        0.622 .* e ./ (sp_Pa .- 0.378 .* e)
    end

end # module ERA5Forcings