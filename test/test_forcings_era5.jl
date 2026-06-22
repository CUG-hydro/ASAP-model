# Tests for `ASAP.ERA5Forcings` — translates `fortran/module_forcings.f90`
# subroutines READFORCINGS, READFORCINGSACC, READFORCINGSSNOW,
# READFORCINGSSOILT, READTOPOERA5 and READLAICLIM into Julia.
#
# The tests do not depend on real ERA5 files: they write small mock
# NetCDF files in a temporary directory using `NCDatasets`, then call
# each reader and verify shapes, units and key transformations.
using Test
using NCDatasets

# Functions are re-exported by the `ASAP` module.
using ASAP: read_hourly_forcings, read_accumulated_forcings,
            read_snow, read_soil_temps, read_topo_era5,
            read_lai_climatology, compute_qair

# ---------------------------------------------------------------------------
# Mock data layout — small synthetic grids (no real ERA5 data needed)
# ---------------------------------------------------------------------------
#
# Production path:  global 1440×721 ERA5 grid, regional crop 1067:1316 × 298:587
#                   (250 × 290); ERA5-Land snow uses 3600×1801, crop 2666:3285 × 747:1466.
# Test path:        tiny 30×30 global grids, crop 5:14 × 5:14 (10×10).
#                   This makes each mock file < 1 MB instead of > 150 MB.
#
# The reader functions accept `xrange`/`yrange` kwargs (added to ERA5.jl)
# that default to the production constants, so default behaviour is unchanged.
# Here we pass the small test ranges explicitly.

const GX = 30          # small mock global x dimension
const GY = 30          # small mock global y dimension
const XR = 5:14        # regional crop x range (10 elements)
const YR = 5:14        # regional crop y range (10 elements)
const RX = 10          # length(XR)
const RY = 10          # length(YR)
const NHOURS = 24

# ERA5-Land snow uses a different (0.1°) global grid; same small approach.
const GX_SNOW = 30
const GY_SNOW = 30
const XR_SNOW = 5:14
const YR_SNOW = 5:14
const RX_SNOW = 10     # length(XR_SNOW)
const RY_SNOW = 10     # length(YR_SNOW)

# Model grid for accumulation/interpolation
const NX = 12
const NY = 8

const TMPDIR = mktempdir()

# ---------------------------------------------------------------------------
# Mock-file writers  (all use small GX×GY global grids)
# ---------------------------------------------------------------------------

"""
    write_mock_era5_3d(dir, dirname, ncvar, value) -> String

Write a 3-D mock ERA5 file `ERA5_<dirname>_<date>.nc` of shape
`(GX, GY, 24)`, containing the NetCDF variable `ncvar`. Every voxel
(including outside the crop) is filled with `value`.
"""
function write_mock_era5_3d(dir::String, dirname::String, ncvar::String, value::Real)
    subdir = joinpath(dir, dirname)
    mkpath(subdir)
    path = joinpath(subdir, "ERA5_$(dirname)_20200101.nc")
    NCDataset(path, "c") do ds
        defDim(ds, "x", GX)
        defDim(ds, "y", GY)
        defDim(ds, "time", NHOURS)
        v = defVar(ds, ncvar, Float64, ("x", "y", "time"))
        v[:, :, :] = fill(Float64(value), GX, GY, NHOURS)
    end
    return path
end

"""
    write_mock_era5_3d_in_region(dir, dirname, ncvar, value) -> String

Same as `write_mock_era5_3d` but writes `value` only inside the
regional crop `XR × YR`, and zero outside.
"""
function write_mock_era5_3d_in_region(dir::String, dirname::String, ncvar::String, value::Real)
    subdir = joinpath(dir, dirname)
    mkpath(subdir)
    path = joinpath(subdir, "ERA5_$(dirname)_20200101.nc")
    NCDataset(path, "c") do ds
        defDim(ds, "x", GX)
        defDim(ds, "y", GY)
        defDim(ds, "time", NHOURS)
        v = defVar(ds, ncvar, Float64, ("x", "y", "time"))
        v[:, :, :] = fill(0.0, GX, GY, NHOURS)
        for k in 1:NHOURS
            v[XR, YR, k] = fill(Float64(value), RX, RY)
        end
    end
    return path
end

"""
    write_mock_era5_precip(dir) -> String

Write a mock PRECIP file with `tp` (total precipitation, m) and `sf`
(snowfall, m). Inside the regional subdomain, `tp = 0.005`, `sf = 0.001`,
so the resulting rain is `(0.005 - 0.001) * 1000 = 4.0` mm.
"""
function write_mock_era5_precip(dir::String)
    mkpath(joinpath(dir, "PRECIP"))
    path = joinpath(dir, "PRECIP", "ERA5_PRECIP_20200101.nc")
    NCDataset(path, "c") do ds
        defDim(ds, "x", GX)
        defDim(ds, "y", GY)
        defDim(ds, "time", NHOURS)
        tp = defVar(ds, "tp", Float64, ("x", "y", "time"))
        sf = defVar(ds, "sf", Float64, ("x", "y", "time"))
        tp[:, :, :] = fill(0.0, GX, GY, NHOURS)
        sf[:, :, :] = fill(0.0, GX, GY, NHOURS)
        for k in 1:NHOURS
            tp[XR, YR, k] = fill(0.005, RX, RY)   # m
            sf[XR, YR, k] = fill(0.001, RX, RY)   # m
        end
    end
    return path
end

"""
    write_mock_era5_radiation(dir) -> String

Write a mock RADIATION file holding both `ssrd` (J/m²) and `sr` (J/m²)
in the regional subdomain. Both filled with the same value (3600 J/m²)
so the post-conversion flux is 1.0 W/m².
"""
function write_mock_era5_radiation(dir::String, value::Real=3600.0)
    mkpath(joinpath(dir, "RADIATION"))
    path = joinpath(dir, "RADIATION", "ERA5_RADIATION_20200101.nc")
    NCDataset(path, "c") do ds
        defDim(ds, "x", GX)
        defDim(ds, "y", GY)
        defDim(ds, "time", NHOURS)
        ssrd = defVar(ds, "ssrd", Float64, ("x", "y", "time"))
        sr   = defVar(ds, "sr",   Float64, ("x", "y", "time"))
        ssrd[:, :, :] = fill(0.0, GX, GY, NHOURS)
        sr[:, :, :] = fill(0.0, GX, GY, NHOURS)
        for k in 1:NHOURS
            ssrd[XR, YR, k] = fill(Float64(value), RX, RY)
            sr[XR, YR, k]   = fill(Float64(value), RX, RY)
        end
    end
    return path
end

"""
    write_mock_era5_soilt(dir) -> String

Write a mock SOILTEMP file holding stl1..stl4 inside the regional
subdomain. Layers 1 and 2 are filled with `T_cold` (frozen), layers
3 and 4 with `T_warm` (liquid).
"""
function write_mock_era5_soilt(dir::String; T_cold::Real=270.0, T_warm::Real=280.0)
    mkpath(joinpath(dir, "SOILTEMP"))
    path = joinpath(dir, "SOILTEMP", "ERA5_SOILTEMP_20200101.nc")
    NCDataset(path, "c") do ds
        defDim(ds, "x", GX)
        defDim(ds, "y", GY)
        defDim(ds, "time", NHOURS)
        for (k, name) in enumerate(("stl1", "stl2", "stl3", "stl4"))
            v = defVar(ds, name, Float64, ("x", "y", "time"))
            v[:, :, :] = fill(0.0, GX, GY, NHOURS)
            target = k <= 2 ? T_cold : T_warm
            for h in 1:NHOURS
                v[XR, YR, h] = fill(Float64(target), RX, RY)
            end
        end
    end
    return path
end

"""
    write_mock_era5_snow(dir) -> String

Write a mock ERA5-Land snowmelt file of size `GX_SNOW × GY_SNOW × 24`.
Inside the regional subdomain `XR_SNOW × YR_SNOW`, fill `smlt` (m)
with `smlt` value, so after clamp + ×1000 the value is `smlt * 1000` mm.
"""
function write_mock_era5_snow(dir::String; smlt::Real=0.002)
    mkpath(joinpath(dir, "SNOW"))
    path = joinpath(dir, "SNOW", "ERA5L_snow_20200101.nc")
    NCDataset(path, "c") do ds
        defDim(ds, "x", GX_SNOW)
        defDim(ds, "y", GY_SNOW)
        defDim(ds, "time", NHOURS)
        v = defVar(ds, "smlt", Float64, ("x", "y", "time"))
        v[:, :, :] = fill(0.0, GX_SNOW, GY_SNOW, NHOURS)
        for h in 1:NHOURS
            v[XR_SNOW, YR_SNOW, h] = fill(Float64(smlt), RX_SNOW, RY_SNOW)
        end
    end
    return path
end

"""
    write_mock_era5_static(path, z_value) -> path

Write a static ERA5 file containing the `z` (geopotential) variable
of shape `(1440, 721)`. The test value is `z_value` (m²/s²) inside
the global grid; the resulting `topoera5` after `÷ 9.81` is
`z_value / 9.81`.
"""
function write_mock_era5_static(path::String, z_value::Real)
    NCDataset(path, "c") do ds
        defDim(ds, "x", 1440)
        defDim(ds, "y", 721)
        v = defVar(ds, "z", Float64, ("x", "y"))
        v[:, :] = fill(Float64(z_value), 1440, 721)
    end
    return path
end

"""
    write_mock_lai_clim(path) -> path

Write a 3-D LAI climatology with dimensions `(nx=10, ny=8, time=12)`.
Each month has a constant value: month `m` holds `0.1 * m`.
"""
function write_mock_lai_clim(path::String; nx::Int=10, ny::Int=8)
    NCDataset(path, "c") do ds
        defDim(ds, "x", nx)
        defDim(ds, "y", ny)
        defDim(ds, "month", 12)
        v = defVar(ds, "LAI", Float64, ("x", "y", "month"))
        for m in 1:12
            v[:, :, m] = fill(0.1 * m, nx, ny)
        end
    end
    return path
end

# ---------------------------------------------------------------------------
# Set up a model grid centred on the regional subdomain so that the
# bilinear interpolation is well-defined.
# ---------------------------------------------------------------------------
#
# The reader functions use a fixed SW corner (xswlat=-56.5, xswlon=266.5)
# and spacing (dxv=dyv=0.25) inherited from the Fortran code.
# With the small test crop (rx=ry=10), model-grid indices are clamped
# to [1, rx-1] / [1, ry-1].  Because the entire cropped array is filled
# with a constant value, bilinear interpolation of a constant field always
# returns that constant regardless of where the clamped index falls.
# Therefore any lat/lon that lands on a land cell will return the expected
# analytic value.
function make_model_grid()
    lats = zeros(Float64, NX, NY)
    lons = zeros(Float64, NX, NY)
    # Use coordinates near the centre of the regional subdomain.
    # Subdomain SW: lat = -56.5, lon = 266.5. Centre ≈ (lat + RY*0.25/2)
    for j in 1:NY, i in 1:NX
        lats[i, j] = -56.5 + (j - 1) * 2.0   # ~every 2° lat
        lons[i, j] = 266.5 + (i - 1) * 2.0   # ~every 2° lon
    end
    landmask = ones(Int, NX, NY)
    return lats, lons, landmask
end

# ---------------------------------------------------------------------------
# Test sets
# ---------------------------------------------------------------------------

@testset "ERA5Forcings — compute_qair analytic values" begin
    # Tetens at 10°C, 1013.25 hPa → saturation mixing ratio ≈ 7.6 g/kg.
    # At 100% RH (Td = T) we expect qair ≈ 0.0075.
    d2m = fill(283.15, 2, 3)
    sp  = fill(101325.0, 2, 3)
    qair = compute_qair(d2m, sp)
    @test all(abs.(qair .- 0.0075) .< 5e-4)
    @test size(qair) == (2, 3)
    @test eltype(qair) <: AbstractFloat

    # d2m = 0 → exp(0) = 1, e = 6.112 hPa = 611.2 Pa; sp - 0.378e > 0
    # for normal sp, so qair is well-defined and > 0.
    qair_zero = compute_qair(fill(273.15, 1, 1), fill(101325.0, 1, 1))
    @test qair_zero[1, 1] > 0
end

@testset "ERA5Forcings — read_hourly_forcings" begin
    # Make a fresh tempdir so this test does not depend on others.
    dir = mktempdir()
    write_mock_era5_3d_in_region(dir, "WIND",     "ws10", 5.0)       # 5 m/s
    write_mock_era5_3d_in_region(dir, "TEMP",     "t2m",  290.0)     # 17 °C
    write_mock_era5_3d_in_region(dir, "SFCPRESS", "sp",   101325.0)  # Pa
    # d2m such that T-Td = 5 K, so RH ≈ exp((17.67*5)/(243.5+17)) / 1 ≈ 1.27
    # which clamps to 1.0.
    write_mock_era5_3d_in_region(dir, "DEWPOINT", "d2m",  285.0)
    write_mock_era5_precip(dir)                             # tp=0.005, sf=0.001
    write_mock_era5_radiation(dir, 3600.0)                  # ⇒ 1 W/m² after /3600
    write_mock_era5_soilt(dir; T_cold=270.0, T_warm=280.0)

    out = read_hourly_forcings("20200101", dir; xrange=XR, yrange=YR)

    @test size(out.varpack) == (RX, RY, NHOURS, 11)
    @test all(out.wind    .== 5.0)
    @test all(out.tempera .== 290.0)
    @test all(out.press   .== 101325.0)
    # RH = exp((Td-273.15)*17.67/(243.5+(Td-273.15))) /
    #      exp((T -273.15)*17.67/(243.5+(T -273.15)))
    # With T=290, Td=285: numerator = exp(17.67*12/255.5) ≈ 2.287
    # denominator          = exp(17.67*17/260.5)  ≈ 3.165
    # RH ≈ 0.7224
    @test all(abs.(out.rh .- 0.7224) .< 0.01)
    @test all(out.rain .== 4.0)                # (0.005 - 0.001) * 1000 mm
    @test all(out.swdown  .== 3600.0)          # J/m²
    @test all(out.rnetera .== 3600.0)
    @test all(out.st1 .== 270.0)
    @test all(out.st2 .== 270.0)
    @test all(out.st3 .== 280.0)
    @test all(out.st4 .== 280.0)
    # varpack slice order
    @test out.varpack[:, :, 1, 1] == out.wind[:, :, 1]
    @test out.varpack[:, :, 1, 8] == out.st1[:, :, 1]
end

@testset "ERA5Forcings — read_accumulated_forcings" begin
    # Use the same data as the previous test (the unit test in this file
    # only needs varpack and a model grid).
    dir = mktempdir()
    write_mock_era5_3d_in_region(dir, "WIND",     "ws10", 5.0)
    write_mock_era5_3d_in_region(dir, "TEMP",     "t2m",  290.0)
    write_mock_era5_3d_in_region(dir, "SFCPRESS", "sp",   101325.0)
    write_mock_era5_3d_in_region(dir, "DEWPOINT", "d2m",  285.0)
    write_mock_era5_precip(dir)
    write_mock_era5_radiation(dir, 3600.0)
    write_mock_era5_soilt(dir; T_cold=270.0, T_warm=280.0)
    out = read_hourly_forcings("20200101", dir; xrange=XR, yrange=YR)

    lats, lons, landmask = make_model_grid()
    acc = read_accumulated_forcings(1, out.varpack;
                                    lats=lats, lons=lons, landmask=landmask)

    @test size(acc.pcpgl)  == (NX, NY)
    @test size(acc.rshort) == (NX, NY)
    @test size(acc.netrad) == (NX, NY)
    @test all(acc.pcpgl  .== 4.0)        # mm
    @test all(acc.rshort .== 1.0)        # 3600 J/m² ÷ 3600 = 1 W/m²
    @test all(acc.netrad .== 1.0)
end

@testset "ERA5Forcings — read_soil_temps (ice factor)" begin
    # stl1 = stl2 = 270 K (frozen), stl3 = stl4 = 280 K (liquid).
    dir = mktempdir()
    write_mock_era5_3d_in_region(dir, "WIND",     "ws10", 5.0)
    write_mock_era5_3d_in_region(dir, "TEMP",     "t2m",  290.0)
    write_mock_era5_3d_in_region(dir, "SFCPRESS", "sp",   101325.0)
    write_mock_era5_3d_in_region(dir, "DEWPOINT", "d2m",  285.0)
    write_mock_era5_precip(dir)
    write_mock_era5_radiation(dir, 3600.0)
    write_mock_era5_soilt(dir; T_cold=270.0, T_warm=280.0)
    out = read_hourly_forcings("20200101", dir; xrange=XR, yrange=YR)

    nzg = 30  # arbitrary, just needs to be ≥ 14 + 1
    lats, lons, landmask = make_model_grid()
    topo    = fill(500.0, NX, NY)
    topoera = fill(450.0, NX, NY)
    ice = read_soil_temps(1, out.varpack; nzg=nzg, lats=lats, lons=lons,
                          landmask=landmask, topo=topo, topoera=topoera)
    @test size(ice) == (NX, NY, 15)
    @test eltype(ice) == Int8
    # The Fortran code sets icefactor = 1 (frozen) when soilt ≤ 273.15 K.
    # Here stl1 = stl2 = 270 K, so for nzg=30:
    #   nlev 1 → icefactor[i, j, 15] = 1
    #   nlev 2 → icefactor[i, j, 13:14] = 1
    # stl3 = stl4 = 280 K → no change.
    @test all(ice[:, :, 15] .== 1)
    @test all(ice[:, :, 13:14] .== 1)
    @test all(ice[:, :, 1:12] .== 0)
end

@testset "ERA5Forcings — read_snow (mm + non-negativity)" begin
    dir = mktempdir()
    write_mock_era5_snow(dir; smlt=0.002)  # ⇒ 2.0 mm after × 1000

    lats, lons, landmask = make_model_grid()
    snow = read_snow("20200101", dir;
                     lats=lats, lons=lons, landmask=landmask,
                     xrange=XR_SNOW, yrange=YR_SNOW)
    @test size(snow) == (NX, NY)
    @test all(snow .== 2.0)
    @test all(snow .>= 0)
end

@testset "ERA5Forcings — read_topo_era5 (static)" begin
    path = joinpath(TMPDIR, "static_test.nc")
    z_value = 9810.0  # ⇒ 1000 m after ÷ 9.81
    write_mock_era5_static(path, z_value)

    lats, lons, landmask = make_model_grid()
    # We need lat/lon that fall inside the global ERA5 grid
    # (lat ∈ [-90, 90], lon ∈ [0, 360]).
    for j in 1:NY, i in 1:NX
        lats[i, j] = -45.0 + 2 * (j - 1)   # -45..-31
        lons[i, j] = 280.0 - 2 * (i - 1)   # 280..258 (wrapped > 0 OK)
    end
    topoera5 = read_topo_era5(path; lats=lats, lons=lons, landmask=landmask)
    @test size(topoera5) == (NX, NY)
    # After the Fortran L839-L843 north-up flip, the value remains
    # z_value / 9.81 across the whole grid (since the input is constant).
    @test all(abs.(topoera5 .- 1000.0) .< 1.0)
end

@testset "ERA5Forcings — read_lai_climatology" begin
    path = joinpath(TMPDIR, "lai_clim_test.nc")
    write_mock_lai_clim(path; nx=10, ny=8)
    lai = read_lai_climatology(path, 7; target_shape=(10, 8))
    @test size(lai) == (10, 8)
    @test all(abs.(lai .- 0.7) .< 1e-10)   # 0.1 * month (FP: 0.1*7 ≈ 0.7)

    # Wrong target_shape should fail with DimensionMismatch.
    @test_throws DimensionMismatch read_lai_climatology(path, 1;
                                                       target_shape=(5, 8))
end
