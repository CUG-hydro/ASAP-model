# Unit tests for `example/regional_example.jl::generate_mock_dataset`.
#
# The example file itself is a manual driver and is not `include`-d by
# runtests.jl (examples are designed to be invoked from the CLI). We only
# unit-test the dataset generator here, verifying that:
#
# 1. The function returns a NamedTuple with all required fields.
# 2. The generated NetCDF files exist on disk and have the correct shapes.
# 3. Round-tripping `read_initial` and `read_wtdnc` over the mock files
#    yields the expected post-processed values (mirrors the Fortran
#    READINITIAL / READWTDNC post-processing).
using ASAP, Test, NCDatasets

# Include the example file in a sandboxed namespace so we can reach the
# internal `generate_mock_dataset` without launching `main(...)`.
const REGIONAL_EXAMPLE_PATH = joinpath(@__DIR__, "..", "example", "regional_example.jl")
isfile(REGIONAL_EXAMPLE_PATH) ||
    error("example/regional_example.jl not found at $REGIONAL_EXAMPLE_PATH")

# Evaluate the example in a module so the `if abspath(PROGRAM_FILE) == @__FILE__`
# entry guard at the bottom does not trigger.
const REGIONAL_EXAMPLE = Module(:RegionalExampleForTest)
Core.eval(REGIONAL_EXAMPLE, :(using NCDatasets, Random, Dates, Printf, Statistics))
Core.eval(REGIONAL_EXAMPLE, Meta.parseall(read(REGIONAL_EXAMPLE_PATH, String)))
generate_mock_dataset::Function = getfield(REGIONAL_EXAMPLE, :generate_mock_dataset)

@testset "regional_example: generate_mock_dataset 结构与字段" begin
    date, nx, ny, nzg = "20200101", 3, 3, 40
    paths = generate_mock_dataset(date, nx, ny, nzg)

    @test isdir(paths.root)
    for k in (:static, :wtd, :era5_wind, :era5_temp, :era5_dewpoint,
              :era5_press, :era5_strd, :era5_ssrd, :era5_soilt, :era5_tp,
              :lai_clim)
        @test k in keys(paths)
        @test isfile(getfield(paths, k))
    end
end

@testset "regional_example: static.nc 内容与 read_initial 往返" begin
    nx, ny = 4, 4
    paths = generate_mock_dataset("20200202", nx, ny, 40)

    initial = read_initial(paths.static)
    @test size(initial.soiltxt)  == (nx, ny)
    @test size(initial.topo)     == (nx, ny)
    @test size(initial.fdepth)   == (nx, ny)
    @test size(initial.landmask) == (nx, ny)

    # 新 mock：soiltxt 横向梯度 1..13；topo 含山脊 0..500 m；
    # landmask 由 topo ∈ (0, 480) 派生；fdepth 与 soiltxt 负相关
    @test all(1 .<= initial.soiltxt .<= 13)
    @test all(initial.topo       .>= -50)    # 含噪声
    @test all(initial.topo       .<= 500)
    @test all(initial.fdepth     .>= 0.5)
    @test all(initial.fdepth     .<= 3.0)
    @test all((initial.landmask .== 0) .| (initial.landmask .== 1))
end

@testset "regional_example: wtd.nc 符号约定（Fortran: wtd = min(-raw, 0)）" begin
    nx, ny = 3, 3
    paths = generate_mock_dataset("20200303", nx, ny, 40)

    wtd = read_wtdnc(paths.wtd)
    @test size(wtd) == (nx, ny)
    @test all(wtd .<= 0)                   # Fortran min(-raw, 0)
    @test all(-5.0 .<= wtd .<= -0.1)       # 新 mock: 原始 0.2..4.5 m → 读出负值
end

@testset "regional_example: ERA5 强迫变量维度" begin
    nx, ny, nzg = 3, 3, 40
    paths = generate_mock_dataset("20200404", nx, ny, nzg)

    NCDataset(paths.era5_wind) do ds
        @test haskey(ds, "wind")
        @test size(ds["wind"]) == (nx, ny, 24)
    end
    NCDataset(paths.era5_temp) do ds
        @test haskey(ds, "t2m")
        @test size(ds["t2m"]) == (nx, ny, 24)
    end
    NCDataset(paths.era5_press) do ds
        @test haskey(ds, "sp")
        @test size(ds["sp"]) == (nx, ny, 24)
    end
    NCDataset(paths.era5_tp) do ds
        @test haskey(ds, "tp")
        @test size(ds["tp"]) == (nx, ny, 24)
    end
    NCDataset(paths.era5_soilt) do ds
        @test haskey(ds, "STL1")
        @test size(ds["STL1"]) == (nx, ny, 24, 4)
    end
    NCDataset(paths.lai_clim) do ds
        @test haskey(ds, "lai")
        @test size(ds["lai"]) == (nx, ny, 12)
    end
end

@testset "regional_example: 端到端 smoke test（短时步 + mock 数据）" begin
    # 调用主流程 --duration 1；预期无 UndefVarError，输出包含完成字符串
    # 用 `--mock 1 --nx 3 --ny 3 --nzg 40 --duration 1` 跑最短路径，
    # 并把输出隔离到 base_dir 避免污染 /tmp。
    base_dir = mktempdir()
    out_path = joinpath(base_dir, "out.nc")

    cmd = `$(Base.julia_cmd()) --project=$(joinpath(@__DIR__, "..")) $(REGIONAL_EXAMPLE_PATH) --mock 1 --duration 1 --nx 3 --ny 3 --nzg 40 --out $(out_path)`
    try
        output = read(cmd, String)
        @test occursin("区域应用示例运行完成", output)
        @test occursin("Mock 模式", output)
        @test occursin("时步循环", output)
        @test isfile(out_path)  # 输出 NetCDF 应已落盘
    catch e
        rethrow(e)  # @test 框架会把异常记为 error 并打印 stacktrace
    end
end
