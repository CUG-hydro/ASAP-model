"""
ASAP 模型 Regional 区域应用示例
================================

演示从 NetCDF 真实数据驱动的端到端流程：

1. 读静态场 (`read_initial`)              ← `fortran/module_io.f90:READINITIAL`
2. 读初始水位 (`read_wtdnc`)               ← `fortran/module_io.f90:READWTDNC`
3. 计算平衡含水量 (`eqsoilmoisturetheor`)   ← `fortran/module_initial.f90:EQSOILMOISTUREtheor`
4. 时步循环（小时）
   - 读 ERA5 强迫（`read_hourly_forcings`；本示例内嵌最小实现）
   - 调 `rootdepth_main`                    ← `fortran/module_rootdepth.f90:ROOTDEPTH`
5. 输出最终状态

# 用法

```bash
# Mock 模式（默认，无外部数据；自动在 mktempdir() 生成数据集）
julia --project example/regional_example.jl --mock 1 --duration 3

# 真实数据模式（需提供 NetCDF 文件）
julia --project example/regional_example.jl \\
      --date 20200101 --duration 24 \\
      --static /data/ASAP/static.nc \\
      --wtd    /data/ASAP/wtd.nc \\
      --era5   /data/ERA5 \\
      --out    /tmp/regional_out.nc
```

# 命令行参数

- `--date YYYYMMDD`     起始日期（默认 20200101）
- `--duration N`        模拟小时数（默认 3）。**支持多日**：每跨过 24h 自动
                        推进一天，真实数据模式按滚动日期读取对应强迫文件
                        （需为每一天准备 ERA5 文件）；mock 模式复用单日合成数据
- `--nx N --ny N`       网格大小（默认 3，mock 模式）
- `--nzg N`             土壤层数（默认 40）
- `--mock 0|1`          1 = 合成数据（默认），0 = 读真实 NetCDF
- `--static PATH`       静态场 NetCDF（真实模式必填）
- `--wtd PATH`          初始水位 NetCDF（真实模式必填）
- `--era5 DIR`          ERA5 数据根目录（真实模式必填）
- `--out PATH`          输出 NetCDF（默认 `/tmp/asap_regional_<date>.nc`）

# 风格约定

- 与 `example/complete_example.jl` 一致：直接 `using ASAP`、所有状态原地修改
- `rootdepth_main` 签名与 `test/test_rootdepth.jl::make_rootdepth_inputs` 对齐
- mock 模式 100% 自给自足，不依赖任何外部文件
"""

using ASAP
using NCDatasets
using Statistics
using Dates
using Printf
using Random

# ===========================================================================
# 命令行参数解析（轻量、不引新依赖）
# ===========================================================================

"""
    parse_cli_args(args) -> NamedTuple

极简的 key=value / --flag 解析器。返回 `(; date, duration, nx, ny, nzg,
mock, static, wtd, era5, out)`。
"""
function parse_cli_args(args::Vector{String})
    kw = Dict{String,Any}(
        "date"     => "20200101",
        "duration" => 3,
        "nx"       => 3,
        "ny"       => 3,
        "nzg"      => 40,
        "mock"     => 1,
        "static"   => "",
        "wtd"      => "",
        "era5"     => "",
        "out"      => "",
    )
    i = 1
    while i <= length(args)
        a = args[i]
        if startswith(a, "--")
            key = String(chop(a; head=2, tail=0))
            if haskey(kw, key)
                # --flag 1 / --flag value
                if i < length(args) && !startswith(args[i+1], "--")
                    raw = args[i+1]
                    i += 2
                else
                    # 裸 flag（默认 1）
                    raw = "1"
                    i += 1
                end
                T = typeof(kw[key])
                kw[key] = T === String ? raw : parse(T, raw)
            else
                @warn "Unknown argument --$key (ignored)"
                i += 1
            end
        else
            i += 1
        end
    end
    return (; (Symbol(k) => v for (k, v) in kw)...)
end

# ===========================================================================
# Mock 数据集生成
# ===========================================================================

"""
    generate_mock_dataset(date::String, nx::Int, ny::Int, nzg::Int) -> NamedTuple

在 `mktempdir()` 下生成完整的 mock 数据集，覆盖 ASAP 区域应用所需的所有输入。
所有合成场都含**空间异质 + 时间节律**，与真实区域行为量级一致：

- 静态场：`soiltxt` 1→13（砂→粘）横向梯度；`topo` 山脊形 0..500 m；
  `fdepth` 与 soiltxt 相关（砂土深根 2.5 m，粘土浅根 0.8 m）；
  `landmask = (topo > 0 && topo < 480)`。
- 初始水位：`wtd = -1.0 - 0.005 * topo + 0.5 * noise`，符号约定 `min(-raw, 0)`。
- ERA5 强迫（日循环，UTC）：
  - `t2m`  日正弦 280..300 K，海陆 ±2 K 噪声；
  - `d2m`  `t2m - ΔTd`，ΔTd ∈ 2..6 K（RH ≈ 60..85%）；
  - `wind` 2..6 m/s，地形加速 ×(1 + 0.005·topo)；
  - `sp`   按压高公式 101325·(1 - 0.0065·topo/288.15)^5.255；
  - `ssrd` 0..800 W/m²，仅昼间（hour 6..18）非零，云量衰减 0.7..1.0；
  - `strd` 200 + 2·(t2m − 273)，即 ~Stefan-Boltzmann 粗近似；
  - `tp`   间歇性降水（90% 小时为 0；事件小时指数分布 2..8 mm/h）；
  - `stl1..4`  深−0..4 层滞后于 t2m（lag 1..4 h，年均 285 K）。
- LAI：月度正弦（1 月 0.5，7 月 4.0），日循环微扰 ±0.1。

返回字段：`root / static / wtd / era5_wind / era5_temp / era5_dewpoint /
era5_press / era5_strd / era5_ssrd / era5_soilt / era5_tp / lai_clim`。
"""
function generate_mock_dataset(date::String, nx::Int, ny::Int, nzg::Int)
    root = mktempdir(prefix="asap_regional_mock_")
    Random.seed!(42)  # 可重复

    # ---- 空间场 -----------------------------------------------------------
    # 构造 nx×ny 网格 (i = 1..nx 列, j = 1..ny 行)
    x_norm = range(0, 1; length=nx) |> collect       # length-nx
    y_norm = range(0, 1; length=ny) |> collect       # length-ny

    # soiltxt: 1..13 横向梯度 + 纵向微扰 + ±1 噪声
    soiltxt_raw = [12.0 * (xi - 0.5) + 2.0 * (yj - 0.5) + rand() - 0.5
                   for xi in x_norm, yj in y_norm]   # nx × ny
    soiltxt = clamp.(round.(Int, soiltxt_raw .+ 7.0), 1, 13)

    # topo: 山脊 (cos-bell 沿 x, 缓坡沿 y) + 噪声，范围 0..500 m
    topo = [250.0 * (1.0 - (2xi - 1)^2) + 150.0 * yj + 30.0 * randn()
            for xi in x_norm, yj in y_norm]          # nx × ny

    # fdepth: 与 soiltxt 相关，砂土 2.5 m，粘土 0.8 m
    fdepth = [2.8 - 0.15 * soiltxt[i, j] + 0.05 * randn()
              for i in 1:nx, j in 1:ny]

    # ---- static.nc --------------------------------------------------------
    static_path = joinpath(root, "static.nc")
    NCDataset(static_path, "c") do ds
        defDim(ds, "x", nx); defDim(ds, "y", ny)
        defVar(ds, "STXT", Int32,   ("x", "y"))[:] = Int32.(soiltxt)
        defVar(ds, "topo", Float64, ("x", "y"))[:] = topo
        defVar(ds, "F",    Float64, ("x", "y"))[:] = fdepth
    end

    # ---- wtd.nc (写入正值,符号约定 min(-raw, 0)) ------------------------
    # wtd_raw = -1.0 - 0.005·topo + 0.5·noise；正值 1.0..4.5 m
    wtd_raw = @. max(0.2, 1.0 + 0.005 * topo + 0.3 * randn())
    wtd_path = joinpath(root, "wtd.nc")
    NCDataset(wtd_path, "c") do ds
        defDim(ds, "x", nx); defDim(ds, "y", ny)
        defVar(ds, "WTD", Float64, ("x", "y"))[:] = wtd_raw
    end

    # ---- ERA5 目录 --------------------------------------------------------
    era5_dir = joinpath(root, "ERA5")
    for sub in ("WIND","TEMP","DEWPOINT","SFCPRESS","STRD","SSRD","SOILT","TP")
        mkpath(joinpath(era5_dir, sub))
    end
    mkpath(joinpath(root, "LAI"))

    hour_dim = 24
    n_soil_level = 4
    month_dim = 12

    # ---- 时间场 (3-D) -----------------------------------------------------
    # t2m(i,j,h): 日正弦 中心 290 K, 振幅 8 K, 海陆 ±2 K
    T_base = 290.0 .+ 2.0 .* randn(nx, ny)            # 海陆差异 (nx × ny)
    T_amp  = 8.0
    hours  = collect(1:hour_dim)
    # phase = (h - 9) / 12 → h=9 (UTC 09) 峰值, h=21 谷值
    diurnal = T_amp .* sin.(π .* (hours' .- 9) ./ 12)         # 1 × hour_dim
    t2m = repeat(T_base, 1, 1, hour_dim) .+ reshape(diurnal, 1, 1, hour_dim)
    # d2m = t2m - ΔTd，ΔTd ∈ 2..6 K (RH ≈ 60..85%)
    delta_td = 3.0 .+ 2.0 .* rand(nx, ny)
    d2m = t2m .- reshape(delta_td, nx, ny, 1)

    # wind: 2..6 m/s, 地形加速
    topo_factor = [1.0 + 0.003 * topo[i, j] for i in 1:nx, j in 1:ny]
    wind_base   = clamp.(3.0 .+ 0.8 .* randn(nx, ny), 1.0, 6.0)
    wind_diag   = 1.0 .+ 0.4 .* sin.(2π .* hours' ./ 24)
    wind        = repeat(wind_base .* topo_factor, 1, 1, 1) .* reshape(wind_diag, 1, 1, hour_dim)

    # sp: 压高公式
    sp_static = @. 101325.0 * (1.0 - 0.0065 * topo / 288.15)^5.255
    sp = repeat(sp_static, 1, 1, hour_dim)

    # ssrd: 0..800 W/m², 仅昼间 (h=6..18) 非零, 云量 0.7..1.0
    cloud = 0.7 .+ 0.3 .* rand(nx, ny)
    ssrd_solar = @. max(0.0, 800.0 * sin(π * (hours' - 6) / 12))
    ssrd = repeat(cloud, 1, 1, hour_dim) .* reshape(ssrd_solar, 1, 1, hour_dim)

    # strd: 200 + 2·(t2m - 273), 全天
    strd = @. 200.0 + 2.0 * (t2m - 273.15)
    strd = max.(strd, 100.0)  # 兜底

    # tp: 间歇性降水 (90% 小时 = 0, 10% 小时 = 1..5 mm/h → m/h)
    rng_tp = Random.MersenneTwister(43)
    tp_hourly_m = zeros(nx, ny, hour_dim)
    for h in 1:hour_dim, j in 1:ny, i in 1:nx
        if rand(rng_tp) < 0.10  # 10% 小时有降水
            # 1..5 mm/h (m/h = mm/h × 0.001)
            tp_hourly_m[i, j, h] = (1.0 + 4.0 * rand(rng_tp)) * 0.001
        end
    end

    # soilt 4 层: 深−0..4 层滞后于 t2m
    # lag [0, 2, 6, 24] 小时, deep 层年均 285 K (深−常数)
    stl = zeros(nx, ny, hour_dim, n_soil_level)
    for k in 1:n_soil_level
        lag = [0, 2, 6, 24][k]
        # 平移 t2m: 层 k = t2m 平移 lag 小时 (深层滞后)
        for h in 1:hour_dim
            src_h = mod1(h + lag, hour_dim)
            stl[:, :, h, k] = t2m[:, :, src_h] .- 0.5 .* (k - 1)
        end
    end
    # 强制 stl4 ≈ 285 K (深层年平均)
    stl[:, :, :, 4] .= 285.0

    # ---- 写文件 ----------------------------------------------------------
    era5_wind     = joinpath(era5_dir, "WIND",     "ERA5_wind_speed_$(date).nc")
    era5_temp     = joinpath(era5_dir, "TEMP",     "ERA5_2m_temperature_$(date).nc")
    era5_dewpoint = joinpath(era5_dir, "DEWPOINT", "ERA5_2m_dewpoint_$(date).nc")
    era5_press    = joinpath(era5_dir, "SFCPRESS", "ERA5_surface_pressure_$(date).nc")
    era5_strd     = joinpath(era5_dir, "STRD",     "ERA5_strd_$(date).nc")
    era5_ssrd     = joinpath(era5_dir, "SSRD",     "ERA5_ssrd_$(date).nc")
    era5_soilt    = joinpath(era5_dir, "SOILT",    "ERA5_soil_temps_$(date).nc")
    era5_tp       = joinpath(era5_dir, "TP",       "ERA5_total_precipitation_$(date).nc")

    function write_hourly3d(path, varname, arr)
        NCDataset(path, "c") do ds
            defDim(ds, "x", nx); defDim(ds, "y", ny); defDim(ds, "hour", hour_dim)
            defVar(ds, varname, Float64, ("x", "y", "hour"))[:] = arr
        end
    end
    function write_hourly4d(path, varname, arr)
        NCDataset(path, "c") do ds
            defDim(ds, "x", nx); defDim(ds, "y", ny)
            defDim(ds, "hour", hour_dim); defDim(ds, "soil_level", n_soil_level)
            defVar(ds, varname, Float64, ("x", "y", "hour", "soil_level"))[:] = arr
        end
    end

    write_hourly3d(era5_wind,     "wind",  wind)
    write_hourly3d(era5_temp,     "t2m",   t2m)
    write_hourly3d(era5_dewpoint, "d2m",   d2m)
    write_hourly3d(era5_press,    "sp",    sp)
    write_hourly3d(era5_strd,     "strd",  strd)
    write_hourly3d(era5_ssrd,     "ssrd",  ssrd)
    write_hourly4d(era5_soilt,    "STL1",  stl)
    write_hourly3d(era5_tp,       "tp",    tp_hourly_m)

    # ---- LAI 月度气候态 (冬 0.5 夏 4.0 正弦) ---------------------------
    lai_path = joinpath(root, "LAI", "lai_clim_01.nc")
    NCDataset(lai_path, "c") do ds
        defDim(ds, "x", nx); defDim(ds, "y", ny); defDim(ds, "month", month_dim)
        v = defVar(ds, "lai", Float64, ("x", "y", "month"))
        lai_clim = zeros(nx, ny, month_dim)
        for m in 1:month_dim
            # 1 月(谷) 0.5, 7 月(峰) 4.0
            lai_clim[:, :, m] .= 2.25 .+ 1.75 .* sin(π .* (m - 4) ./ 6)
        end
        v[:] = lai_clim
    end

    return (
        root       = root,
        static     = static_path,
        wtd        = wtd_path,
        era5_wind  = era5_wind,
        era5_temp  = era5_temp,
        era5_dewpoint = era5_dewpoint,
        era5_press = era5_press,
        era5_strd  = era5_strd,
        era5_ssrd  = era5_ssrd,
        era5_soilt = era5_soilt,
        era5_tp    = era5_tp,
        lai_clim   = lai_path,
    )
end

# ===========================================================================
# ERA5 小时强迫读取（mock stub：未来由 src/Forcings/ERA5.jl 替代）
# ===========================================================================

"""
    read_hourly_forcings(hour::Int, paths) -> NamedTuple

读取指定小时的所有 ERA5 强迫变量，返回 NamedTuple。
为保持本示例端到端可跑，**直接读 NetCDF 变量并落盘**；
未来 `src/Forcings/ERA5.jl::ERA5Forcings` 落地后可替换为本项目导出符号。
"""
function read_hourly_forcings(hour::Int, paths)
    read_hourly2d(path, varname) = begin
        NCDataset(path) do ds
            Array{Float64}(ds[varname][:, :, hour])
        end
    end

    wind  = read_hourly2d(paths.era5_wind,  "wind")
    temp  = read_hourly2d(paths.era5_temp,  "t2m")
    d2m   = read_hourly2d(paths.era5_dewpoint, "d2m")
    press = read_hourly2d(paths.era5_press, "sp")
    strd  = read_hourly2d(paths.era5_strd,  "strd")
    ssrd  = read_hourly2d(paths.era5_ssrd,  "ssrd")
    tp    = read_hourly2d(paths.era5_tp,    "tp")

    # 由露点温度近似比湿（Tetens 公式）：e = 6.112 * exp(17.67*(Td-273.15)/(Td-29.65))
    # q = 0.622 e / (p - 0.378 e)  （单位 kg/kg）
    function dewpoint_to_qair(d2m::Matrix{Float64}, press::Matrix{Float64})
        nx, ny = size(d2m)
        q = similar(d2m)
        for j in 1:ny, i in 1:nx
            Td = d2m[i, j]
            e  = 6.112 * exp(17.67 * (Td - 273.15) / (Td - 29.65)) * 100.0  # Pa
            q[i, j] = 0.622 * e / (press[i, j] - 0.378 * e)
        end
        return q
    end
    qair = dewpoint_to_qair(d2m, press)

    # 由短波+长波估算净辐射（粗略）：netrad = ssrd*(1-alb) - strd
    # 短草 albedo = 0.23
    netrad = ssrd .* (1.0 - 0.23) .- strd

    # 降水由 m/h 转 mm/h
    precip = tp .* 1000.0

    # 短波 rshort 直接取 ssrd（无晴空区分）
    rshort = ssrd

    # LAI：取当月气候态（第 1 月 → 索引 1）
    nx, ny = size(temp)
    lai = NCDataset(paths.lai_clim) do ds
        Array{Float64}(ds["lai"][:, :, 1])  # 全年 1 月份
    end

    return (; wind, temp, qair, press, netrad, rshort, precip, lai)
end

# ===========================================================================
# ERA5 路径构造（按日期）
# ===========================================================================

"""
    era5_paths_for(era5_root::String, date::String) -> NamedTuple

按「变量分目录、每日一文件」约定，为给定日期 `date`（YYYYMMDD）构造 8 个
ERA5 强迫文件路径。用于真实数据模式下的**多日区域运行**：时步循环每跨过
24 小时即推进一天，调用本函数重建当日强迫路径。

返回字段：`era5_wind / era5_temp / era5_dewpoint / era5_press / era5_strd /
era5_ssrd / era5_soilt / era5_tp`（与 `generate_mock_dataset` 保持一致）。
"""
function era5_paths_for(era5_root::String, date::String)
    return (
        era5_wind     = joinpath(era5_root, "WIND",     "ERA5_wind_speed_$(date).nc"),
        era5_temp     = joinpath(era5_root, "TEMP",     "ERA5_2m_temperature_$(date).nc"),
        era5_dewpoint = joinpath(era5_root, "DEWPOINT", "ERA5_2m_dewpoint_$(date).nc"),
        era5_press    = joinpath(era5_root, "SFCPRESS", "ERA5_surface_pressure_$(date).nc"),
        era5_strd     = joinpath(era5_root, "STRD",     "ERA5_strd_$(date).nc"),
        era5_ssrd     = joinpath(era5_root, "SSRD",     "ERA5_ssrd_$(date).nc"),
        era5_soilt    = joinpath(era5_root, "SOILT",    "ERA5_soil_temps_$(date).nc"),
        era5_tp       = joinpath(era5_root, "TP",       "ERA5_total_precipitation_$(date).nc"),
    )
end

# ===========================================================================
# 主流程
# ===========================================================================

function main(args::Vector{String})
    println("ASAP Regional 区域应用示例")
    println(repeat("=", 50))

    # 1. 解析参数
    println("\n[1/6] 解析命令行参数")
    println(repeat("-", 40))
    cfg = parse_cli_args(args)
    cfg.mock == 1 && (cfg = merge(cfg, (out = isempty(cfg.out) ?
        joinpath(tempdir(), "asap_regional_$(cfg.date).nc") : cfg.out,)))
    cfg.mock == 0 && isempty(cfg.out) && (cfg = merge(cfg, (out =
        joinpath(tempdir(), "asap_regional_$(cfg.date).nc"),)))

    @show cfg

    # 2. 准备数据
    println("\n[2/6] 准备数据")
    println(repeat("-", 40))
    if cfg.mock == 1
        println("Mock 模式：在 mktempdir() 生成合成 NetCDF 数据集")
        paths = generate_mock_dataset(cfg.date, cfg.nx, cfg.ny, cfg.nzg)
        println("  临时根目录: $(paths.root)")
        @show paths
    else
        println("真实数据模式")
        for p in (cfg.static, cfg.wtd, cfg.era5)
            isempty(p) && error("--mock 0 时必须提供 --static / --wtd / --era5")
        end
        isfile(cfg.static) || error("静态场文件不存在: $(cfg.static)")
        isfile(cfg.wtd)    || error("初始水位文件不存在: $(cfg.wtd)")
        isdir(cfg.era5)    || error("ERA5 目录不存在: $(cfg.era5)")
        paths = (; root = dirname(cfg.static),
                   static = cfg.static,
                   wtd = cfg.wtd,
                   era5_paths_for(cfg.era5, cfg.date)...,
                   lai_clim = joinpath(dirname(cfg.static), "..", "LAI", "lai_clim_01.nc"))
        for (k, v) in pairs(paths)
            (k == :root) && continue
            isfile(v) || error("ERA5 强迫文件缺失: $k → $v")
        end
    end

    # 3. 静态场初始化
    println("\n[3/6] 读取静态场并初始化土壤参数")
    println(repeat("-", 40))
    is, ie = 1, cfg.nx
    js, je = 1, cfg.ny
    nzg    = cfg.nzg
    nx, ny = ie, je

    initial = read_initial(paths.static)
    soiltxt_raw = initial.soiltxt
    topo        = initial.topo
    fdepth      = initial.fdepth
    landmask    = initial.landmask

    println("  网格 $(nx)×$(ny)，土壤层数 $nzg")
    println("  土壤类型范围: $(minimum(soiltxt_raw))..$(maximum(soiltxt_raw))")
    println("  海拔范围: $(minimum(topo))..$(maximum(topo)) m")
    println("  陆地掩码陆地格点数: $(sum(landmask))")

    # 读初始水位
    wtd_in = read_wtdnc(paths.wtd)
    println("  初始水位: $(round(mean(wtd_in), digits=3)) m")

    # 初始化土壤层深度
    z₋ₕ, dz = initializesoildepth(nzg)
    println("  土壤层深度: $(round(z₋ₕ[1], digits=3)) ~ $(round(z₋ₕ[end], digits=3)) m")

    # 计算平衡含水量（每个像元调用一次，对齐 Fortran EQSOILMOISTUREtheor）
    println("  计算平衡含水量…")
    smoieq = zeros(Float64, nzg, nx, ny)
    soiltxt = ones(Int, 2, nx, ny)        # rootdepth_main 期望 [2, nx, ny]
    for j in 1:ny, i in 1:nx
        # 陆地格点：以静态场的 soiltxt 替换（取首层）
        nsoil = max(1, min(13, Int(soiltxt_raw[i, j])))
        soiltxt[1, i, j] = nsoil
        soiltxt[2, i, j] = nsoil
        smoieq[:, i, j] .= eqsoilmoisturetheor(nzg, nsoil, z₋ₕ, dz,
                                               fdepth[i, j], wtd_in[i, j])
    end
    println("  平衡含水量: $(round(minimum(smoieq), digits=4)) ~ $(round(maximum(smoieq), digits=4))")

    # 4. 初始化状态变量（与 test_rootdepth.jl::make_rootdepth_inputs 对齐）
    println("\n[4/6] 初始化状态变量")
    println(repeat("-", 40))
    freedrain = 1
    maxinactivedays = 240
    Δt    = 3600.0
    steps = 1.0
    hveg  = fill(1.0,  nx, ny)
    veg   = fill(7.0,  nx, ny)             # 7 = 短草

    smoi         = fill(0.25, nzg, nx, ny)
    smoiwtd      = fill(0.35, nx, ny)
    wtd          = copy(wtd_in)
    waterdeficit = zeros(Float64, nx, ny)
    watext       = zeros(Float64, nzg, nx, ny)
    watextdeep   = zeros(Float64, nx, ny)
    rech         = zeros(Float64, nx, ny)
    deeprech     = zeros(Float64, nx, ny)
    et_s         = zeros(Float64, nx, ny)
    et_i         = zeros(Float64, nx, ny)
    et_c         = zeros(Float64, nx, ny)
    et_s_daily   = zeros(Float64, nx, ny)
    et_c_daily   = zeros(Float64, nx, ny)
    transptop    = zeros(Float64, nx, ny)
    intercepstore= zeros(Float64, nx, ny)
    ppacum       = zeros(Float64, nx, ny)
    pppendepth   = zeros(Float64, nx, ny)
    pppendepthold= zeros(Int8,  nx, ny)
    qlat         = zeros(Float64, nx, ny)
    qlatsum      = zeros(Float64, nx, ny)
    qsprings     = zeros(Float64, nx, ny)
    inactivedays = zeros(Int, nzg+2, nx, ny)
    floodheight  = zeros(Float64, nx, ny)
    qrf          = zeros(Float64, nx, ny)
    delsfcwat    = zeros(Float64, nx, ny)
    icefactor    = zeros(Int8,  nx, ny, nzg)
    wtdflux      = zeros(Float64, nx, ny)
    infilk       = zeros(Int8,  nx, ny)
    infilflux    = zeros(Float64, nzg, nx, ny)
    infilfluxday = zeros(Float64, nzg, nx, ny)
    infilcounter = zeros(Int16, nzg, nx, ny)
    upflux       = zeros(Float64, nzg, nx, ny)
    o18          = zeros(Float64, nzg, nx, ny)
    o18ratiopp   = zeros(Float64, nx, ny)
    tempsfc      = fill(293.0, nx, ny)
    qlato18      = zeros(Float64, nx, ny)
    transpo18    = zeros(Float64, nx, ny)

    println("  $(nx)×$(ny) 网格，$(nzg) 层土壤，状态变量已就绪")

    # 5. 时步循环
    println("\n[5/6] 时步循环 (共 $(cfg.duration) 小时)")
    println(repeat("-", 40))
    et_s_accum = zeros(Float64, nx, ny)
    et_c_accum = zeros(Float64, nx, ny)
    et_i_accum = zeros(Float64, nx, ny)
    precip_accum = zeros(Float64, nx, ny)

    try
        base_date = Date(cfg.date, dateformat"yyyymmdd")
        for hour in 1:cfg.duration
            # 多日区域运行：每跨过 24h 推进一天，hour_in_day 在 1..24 循环。
            day_offset  = div(hour - 1, 24)
            hour_in_day = mod1(hour, 24)
            cur_date    = Dates.format(base_date + Day(day_offset), dateformat"yyyymmdd")

            # 真实数据模式按滚动日期重建 ERA5 强迫路径；mock 模式复用单日合成数据。
            paths_h = (cfg.mock == 0 && day_offset > 0) ?
                merge(paths, era5_paths_for(cfg.era5, cur_date)) : paths

            # 读 ERA5 强迫（mock stub；真实数据接入完整 ERA5 模块见 README §6）
            f = read_hourly_forcings(hour_in_day, paths_h)

            # 冰冻因子：此处用土壤温度阈值 273.15 K 简单判定
            # 真实实现需读 ERA5 土壤温度（paths.era5_soilt）并按层映射
            icefactor .= 0  # 暖季无冻土

            # 累计诊断
            et_s_accum   .+= et_s
            et_c_accum   .+= et_c
            et_i_accum   .+= et_i
            precip_accum .+= f.precip

            # 主算法调用（与 test/test_rootdepth.jl 一致）
            rootdepth_main(freedrain, is, ie, js, je, nzg, z₋ₕ, dz, Δt,
                landmask, veg, hveg,
                soiltxt, f.wind, f.temp,
                f.qair, f.press, f.netrad,
                f.rshort, f.lai, f.precip,
                fill(0.0, nx, ny), smoi, smoieq,
                smoiwtd, wtd, waterdeficit,
                watext,
                watextdeep, rech,
                deeprech,
                et_s, et_i, et_c,
                intercepstore, ppacum,
                pppendepth, pppendepthold, qlat,
                qlatsum, qsprings, inactivedays,
                maxinactivedays, fdepth,
                steps, floodheight, qrf,
                delsfcwat, icefactor, wtdflux,
                et_s_daily, et_c_daily, transptop,
                infilk, infilflux, infilfluxday,
                infilcounter, hour,
                o18, o18ratiopp, tempsfc,
                qlato18, transpo18, upflux)

            @printf("  %s h=%2d  ET_s=%7.3f ET_c=%7.3f ET_i=%7.3f P=%7.3f wtd=%7.3f\n",
                cur_date, hour_in_day, mean(et_s), mean(et_c), mean(et_i),
                mean(f.precip), mean(wtd))
        end

        # 6. 输出最终状态
        println("\n[6/6] 输出最终状态")
        println(repeat("-", 40))
        println("  累计降水: $(round(sum(precip_accum), digits=3)) mm")
        println("  累计土壤蒸发: $(round(sum(et_s_accum), digits=3)) mm")
        println("  累计冠层蒸腾: $(round(sum(et_c_accum), digits=3)) mm")
        println("  累计截留蒸发: $(round(sum(et_i_accum), digits=3)) mm")
        println("  最终地下水位: $(round(mean(wtd), digits=3)) m")
        println("  最终表层含水量: $(round(mean(smoi[nzg, :, :]), digits=4))")
        println("  最终底层含水量: $(round(mean(smoi[1,   :, :]), digits=4))")
        println("  最终补给量 (rech): $(round(sum(rech), digits=3)) mm")
        println("  输出 NetCDF: $(cfg.out)")

        # 落盘最小诊断结果
        NCDataset(cfg.out, "c") do ds
            defDim(ds, "x", nx)
            defDim(ds, "y", ny)
            defVar(ds, "et_s", Float64, ("x", "y"))[:] = et_s
            defVar(ds, "et_c", Float64, ("x", "y"))[:] = et_c
            defVar(ds, "et_i", Float64, ("x", "y"))[:] = et_i
            defVar(ds, "wtd",  Float64, ("x", "y"))[:] = wtd
            defVar(ds, "rech", Float64, ("x", "y"))[:] = rech
        end
    catch e
        println("\n✗ 模型运行出现错误:")
        println("  错误类型: $(typeof(e))")
        println("  错误信息: ", sprint(showerror, e))
        if isa(e, MethodError)
            println("  提示：可能 rootdepth_main 子模块签名不匹配")
        end
        rethrow()
    end

    println("\n" * repeat("=", 50))
    println("区域应用示例运行完成")
    println(repeat("=", 50))
    if cfg.mock == 1
        println("注意：当前为 mock 模式，所用数据由 generate_mock_dataset() 在 mktempdir() 自动生成。")
        println("      真实数据准备清单见 README.md（待 P1-C 落地）。")
    end
    return nothing
end

# ===========================================================================
# 入口
# ===========================================================================

if abspath(PROGRAM_FILE) == @__FILE__
    main(ARGS)
end
