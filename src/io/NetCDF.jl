"""
NetCDF I/O for ASAP model (read-only subset).

Translation of `fortran/module_io.f90` reading subset:
- `READINITIAL` (line 10) → [`read_initial`](@ref)
- `READWTDNC`   (line 309) → [`read_wtdnc`](@ref)

P1/P2 subroutines (`READLATLON`, `READVEG`, `READHISTORYNC`, …) are intentionally
out of scope for this PR and will be ported in a follow-up.
"""
module NetCDFIO

using NCDatasets

export read_initial, read_wtdnc

# ---------------------------------------------------------------------------
# Constants mirrored from Fortran
# ---------------------------------------------------------------------------

"""Topography threshold below which a grid cell is treated as ocean (m).

`fortran/module_io.f90` (READINITIAL) sets:
```fortran
landmask = 1
where (topo .lt. -1.e5) landmask = 0
where (topo .lt. -1.e5) topo = 0.
```
"""
const TOPO_MISSING_THRESHOLD = -1.0e5

"""Minimum physically-meaningful root depth factor (m). Smaller values
are clamped to `FDEPTH_MISSING_REPLACE` by READINITIAL.

`fortran/module_io.f90` (READINITIAL) sets:
```fortran
where (fdepth .lt. 1.e-6) fdepth = 100.
```
"""
const FDEPTH_MISSING_THRESHOLD = 1.0e-6
const FDEPTH_MISSING_REPLACE   = 100.0

# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

"""
    read_initial(path::String) -> NamedTuple

Read the ASAP static fields from a NetCDF file:
- `soiltxt`  (Int,      nx × ny): USDA soil texture class 1..13
- `topo`     (Float64,  nx × ny): surface elevation (m)
- `fdepth`   (Float64,  nx × ny): root-depth factor (m)
- `landmask` (Int,      nx × ny): 0 = ocean, 1 = land

Mirrors `fortran/module_io.f90::READINITIAL` (line 10). NetCDF variable
names follow the convention used by the Fortran `nf90_inq_varid` calls:
`STXT` (int16/int32), `F` (fdepth), and a topography variable named `topo`
(the Fortran source reads topography from a binary direct-access file rather
than NetCDF, so the variable name `topo` is the convention adopted here for
the Julia port; rename if the production dataset uses a different name).

# Arguments
- `path::String`: NetCDF file path (e.g. `"/data/static_fields.nc"`)

# Returns
A `NamedTuple` with fields `(soiltxt, topo, fdepth, landmask)`.
"""
function read_initial(path::String)
    soiltxt  = Matrix{Int}(undef, 0, 0)
    topo     = Matrix{Float64}(undef, 0, 0)
    fdepth   = Matrix{Float64}(undef, 0, 0)
    landmask = Matrix{Int}(undef, 0, 0)

    NCDataset(path) do ds
        # soiltxt: STXT (Int16/Int32)
        soiltxt = Array{Int}(ds["STXT"][:, :])
        # topo: topography (Float64)
        topo = Array{Float64}(ds["topo"][:, :])
        # fdepth: F (Float64)
        fdepth_raw = Array{Float64}(ds["F"][:, :])

        # Reproduce Fortran post-processing:
        #   where (fdepth .lt. 1.e-6) fdepth = 100.
        fdepth = ifelse.(fdepth_raw .< FDEPTH_MISSING_THRESHOLD,
                         FDEPTH_MISSING_REPLACE, fdepth_raw)

        #   landmask = 1
        #   where (topo .lt. -1.e5) landmask = 0
        #   where (topo .lt. -1.e5) topo = 0.
        landmask = ifelse.(topo .< TOPO_MISSING_THRESHOLD, 0, 1)
        topo     = ifelse.(topo .< TOPO_MISSING_THRESHOLD, 0.0, topo)
    end

    return (soiltxt=soiltxt, topo=topo, fdepth=fdepth, landmask=landmask)
end

"""
    read_wtdnc(path::String) -> Matrix{Float64}

Read the initial water-table depth (WTD) field from a NetCDF file.

Mirrors `fortran/module_io.f90::READWTDNC` (line 309). The Fortran routine
applies `wtd = min(-raw, 0.)` (clamp to non-positive); this implementation
applies the same transformation so that returned depths are in metres below
the surface, with `0.0` meaning "at the surface" and negative values
for cells where the water table lies below the surface.

# Arguments
- `path::String`: NetCDF file path (e.g. `"/data/wtd.nc"`)

# Returns
- `wtd::Matrix{Float64}` of size `(nx, ny)` containing the water-table depth (m).
"""
function read_wtdnc(path::String)
    wtd = Matrix{Float64}(undef, 0, 0)
    NCDataset(path) do ds
        raw = Array{Float64}(ds["WTD"][:, :])
        # Fortran:  varread(1:n2, 1:n3) = min(-varreadbig(nw:ne, ns:nn), 0.)
        wtd = min.(-raw, 0.0)
    end
    return wtd
end

end # module NetCDFIO
