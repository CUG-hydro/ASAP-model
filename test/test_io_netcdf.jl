using Test
using NCDatasets

# `read_initial` and `read_wtdnc` are re-exported by the `ASAP` module
# (see src/ASAP.jl). They are also reachable as `ASAP.NetCDFIO.read_initial`
# / `ASAP.NetCDFIO.read_wtdnc`.
using ASAP: read_initial, read_wtdnc

const TMPDIR = mktempdir()

@testset "NetCDF I/O — read_wtdnc roundtrip" begin
    nx, ny = 5, 5
    # Positive raw values get negated; negatives stay at 0; min(-raw, 0) yields ≤0.
    wtd_in_raw = fill(1.5, nx, ny)               # ⇒ -1.5 after the min(-·, 0) post-process
    path = joinpath(TMPDIR, "wtd_test.nc")
    NCDataset(path, "c") do ds
        defDim(ds, "x", nx)
        defDim(ds, "y", ny)
        v = defVar(ds, "WTD", Float64, ("x", "y"))
        v[:] = wtd_in_raw
    end

    wtd_out = read_wtdnc(path)
    @test wtd_out == fill(-1.5, nx, ny)
    @test size(wtd_out) == (nx, ny)
    @test eltype(wtd_out) <: AbstractFloat
    @test all(wtd_out .<= 0)  # sign convention: depth below surface is non-positive
end

@testset "NetCDF I/O — read_wtdnc edge cases" begin
    # Negative raw values get clamped to 0 (min(-raw, 0) = 0 when raw > 0? actually
    # for negative raw: -raw > 0, min(>0, 0) = 0; for positive raw: -raw < 0).
    nx, ny = 3, 3
    path = joinpath(TMPDIR, "wtd_edge.nc")
    NCDataset(path, "c") do ds
        defDim(ds, "x", nx)
        defDim(ds, "y", ny)
        v = defVar(ds, "WTD", Float64, ("x", "y"))
        v[:] = [-1.0  0.0  2.0;
                 0.5 -0.5  3.0;
                 0.0  0.0 -1.0]
    end

    wtd_out = read_wtdnc(path)
    @test wtd_out == [0.0  0.0  -2.0;
                      -0.5  0.0  -3.0;
                       0.0  0.0   0.0]
end

@testset "NetCDF I/O — read_initial roundtrip" begin
    nx, ny = 4, 4
    path = joinpath(TMPDIR, "initial_test.nc")
    NCDataset(path, "c") do ds
        defDim(ds, "x", nx)
        defDim(ds, "y", ny)
        # STXT (int), topo (float), F (float)
        v_stxt = defVar(ds, "STXT", Int32, ("x", "y"))
        v_stxt[:] = fill(Int32(5), nx, ny)
        v_topo = defVar(ds, "topo", Float64, ("x", "y"))
        v_topo[:] = fill(100.0, nx, ny)
        v_f = defVar(ds, "F", Float64, ("x", "y"))
        v_f[:] = fill(2.0, nx, ny)
    end

    out = read_initial(path)
    @test out.soiltxt  == fill(5, nx, ny)
    @test out.topo     == fill(100.0, nx, ny)
    @test out.fdepth   == fill(2.0, nx, ny)
    @test out.landmask == ones(Int, nx, ny)
end

@testset "NetCDF I/O — read_initial fdepth clamping" begin
    # Fortran:  where (fdepth .lt. 1.e-6) fdepth = 100.
    nx, ny = 2, 2
    path = joinpath(TMPDIR, "fdepth_clamp.nc")
    NCDataset(path, "c") do ds
        defDim(ds, "x", nx)
        defDim(ds, "y", ny)
        defVar(ds, "STXT", Int32, ("x", "y"))[:] = fill(Int32(3), nx, ny)
        defVar(ds, "topo", Float64, ("x", "y"))[:] = fill(50.0, nx, ny)
        defVar(ds, "F",    Float64, ("x", "y"))[:] = [0.0  1.0e-7;
                                                       0.5  0.0]
    end

    out = read_initial(path)
    # Values < 1e-6 should be replaced with 100.0
    @test out.fdepth == [100.0  100.0;
                          0.5  100.0]
end

@testset "NetCDF I/O — read_initial landmask/topo masking" begin
    # Fortran:  landmask = 1; where (topo < -1e5) landmask = 0, topo = 0.
    nx, ny = 2, 2
    path = joinpath(TMPDIR, "landmask_test.nc")
    NCDataset(path, "c") do ds
        defDim(ds, "x", nx)
        defDim(ds, "y", ny)
        defVar(ds, "STXT", Int32, ("x", "y"))[:] = fill(Int32(1), nx, ny)
        defVar(ds, "topo", Float64, ("x", "y"))[:] = [ 100.0  -1.0e6;
                                                       -5.0e4  200.0]
        defVar(ds, "F",    Float64, ("x", "y"))[:] = fill(1.0, nx, ny)
    end

    out = read_initial(path)
    @test out.landmask == [1  0;
                            1  1]
    @test out.topo     == [100.0  0.0;
                            -5.0e4 200.0]
end
