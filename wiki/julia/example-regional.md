# Regional 区域应用示例

> 源文件：`example/regional_example.jl`（653 行）
> Fortran 来源：`fortran/main.f90`（备用单步调度）+ `fortran/module_io.f90`（`READINITIAL`、`READWTDNC`）+ `fortran/module_forcings.f90`（`READFORCINGS`）+ `fortran/module_initial.f90`（`EQSOILMOISTUREtheor`）+ `fortran/module_rootdepth.f90`（`ROOTDEPTH`）
> 测试：`test/test_regional_example.jl`
> 状态：已摄取（mock + 真实两模式 + 多日滚动）

## 1. 功能概述

`example/regional_example.jl` 是 ASAP-model 的**区域端到端驱动脚本**，把模型全栈拼接成一条可执行流水线：从 NetCDF 静态场、初始水位读入开始，调用 `eqsoilmoisturetheor` 求平衡含水量，再按小时循环读 ERA5 强迫、驱动 `rootdepth_main`，最后输出累计诊断与最终状态。

脚本支持 **mock 模式**（默认，`--mock 1`，由 `generate_mock_dataset()` 在 `mktempdir()` 自动合成完整数据集，无需任何外部文件）与 **真实数据模式**（`--mock 0`，需提供 `--static` / `--wtd` / `--era5`）。多日区域运行（`--duration > 24`）时，每跨过 24 h 自动按滚动日期重建 ERA5 强迫路径，mock 模式复用单日合成数据，真实数据模式按日期切换文件。

## 2. 函数签名

```julia
# 命令行参数
parse_cli_args(args::Vector{String}) -> NamedTuple

# Mock 数据集生成
generate_mock_dataset(date::String, nx::Int, ny::Int, nzg::Int) -> NamedTuple

# 多日 ERA5 路径构造
era5_paths_for(era5_root::String, date::String) -> NamedTuple

# Mock 强迫读取（与 src/Forcings/ERA5.jl::read_hourly_forcings 互不替代）
read_mock_hourly_forcings(hour::Int, paths, month::Int=1) -> NamedTuple

# 主入口
main(args::Vector{String}) -> Nothing
```

`parse_cli_args` 关键字集合：`date` / `duration` / `nx` / `ny` / `nzg` / `mock` / `static` / `wtd` / `era5` / `out`，默认值与文档头部注释一致。

`generate_mock_dataset` 返回字段（11 个）：`root`（临时目录）、`static`（`static.nc` 路径）、`wtd`（`wtd.nc` 路径）、`era5_wind` / `era5_temp` / `era5_dewpoint` / `era5_press` / `era5_strd` / `era5_ssrd` / `era5_soilt` / `era5_tp`（8 个 ERA5 文件路径）、`lai_clim`（LAI 月气候态 NetCDF 路径）。

`read_mock_hourly_forcings` 返回 NamedTuple：`wind` / `temp` / `qair` / `press` / `netrad` / `rshort` / `precip` / `lai`（8 个字段，全部 `(nx, ny) Matrix{Float64}`）。

## 3. 算法 / 公式

### 3.1 主流程（`main`，6 步）

| 步骤 | 内容 | 对应 Fortran |
|---|---|---|
| ① 解析参数 | `parse_cli_args(ARGS)` | — |
| ② 准备数据 | mock 生成 / 真实路径校验 | `fortran/module_driver.f90::READALLINPUTS` |
| ③ 静态场初始化 | `read_initial` + `read_wtdnc` + `initializesoildepth` + `eqsoilmoisturetheor`（逐像元 2-D 循环） | `module_io.f90::READINITIAL` / `READWTDNC` + `module_initial.f90::EQSOILMOISTUREtheor` |
| ④ 状态变量初始化 | 与 `test/test_rootdepth.jl::make_rootdepth_inputs` 对齐（`smoi`、`wtd`、`et_*`、`icefactor`、`o18` 等） | `module_rootdepth.f90::ROOTDEPTH` 入口形参 |
| ⑤ 时步循环 | 每日滚动 `era5_paths_for` + `read_mock_hourly_forcings` + `rootdepth_main` | `module_rootdepth.f90::ROOTDEPTH` |
| ⑥ 输出最终状态 | 累计诊断 + 最小 NetCDF 落盘 | `module_io.f90::WRITEOUTPUT*`（未实现，写出仅示范 6 字段） |

### 3.2 mock 合成场公式

- **`soiltxt`（土壤类型 1..13）**：横向 12 阶 + 纵向 ±2 阶 + ±1 噪声；`clamp(round(xi·12 + yj·2 + rand − 0.5 + 7), 1, 13)`。
- **`topo`（地形 0..500 m）**：山脊 `250·(1 − (2xi−1)²) + 150·yj + 30·randn()`。
- **`fdepth`（根深 0.8..2.5 m）**：`2.8 − 0.15·soiltxt + 0.05·randn()`。
- **`wtd_raw`（写入正值 0.2..4.5 m）**：`max(0.2, 1.0 + 0.005·topo + 0.3·randn())`；下游 `read_wtdnc` 内部执行 `min(-raw, 0)` 翻为非正值。
- **`t2m`（2 m 气温 280..300 K）**：日正弦 `T_base + 8·sin(π·(h − 9)/12)`，`T_base = 290 + 2·randn(nx, ny)`。
- **`d2m`（露点温度）**：`t2m − ΔTd`，`ΔTd ∈ 2..6 K`。
- **`wind`（0.6..21 m/s）**：`wind_base · (1 + 0.003·topo) · (1 + 0.4·sin(2π·h/24))`，`wind_base ∈ [1, 6] m/s`；地形放大后山区可达 21 m/s。
- **`sp`（地面气压）**：压高公式 `101325·(1 − 0.0065·topo/288.15)^5.255`。
- **`ssrd`（短波 0..800 W/m²）**：仅昼间（h=6..18）`max(0, 800·sin(π·(h−6)/12))`，按云量 `0.7..1.0` 缩放。
- **`strd`（长波 ≥ 100 W/m²）**：`max(200 + 2·(t2m − 273.15), 100)`。
- **`tp`（降水 m/h）**：90 % 小时为 0；事件小时均匀分布 `1..5 mm/h → m/h`。
- **`stl1..stl4`（4 层土壤温度 K）**：分别滞后 `[0, 2, 6, 24] h`，深层第 4 层强制 285 K。
- **`lai_clim`（月气候态 0.5..4.0）**：`2.25 + 1.75·sin(π·(m−4)/6)`（1 月谷、7 月峰）。

### 3.3 派生强迫（`read_mock_hourly_forcings` 内部）

| 量 | 公式 | 单位 |
|---|---|---|
| `qair`（比湿） | Tetens：`e = 6.112·exp(17.67·(Td − 273.15)/(Td − 29.65))·100 Pa`；`q = 0.622·e/(p − 0.378·e)` | kg/kg |
| `netrad`（净辐射） | `ssrd·(1 − 0.23) − strd`（草 albedo = 0.23） | W/m² |
| `rshort`（短波） | 直接取 `ssrd`（无晴空区分） | W/m² |
| `precip`（降水） | `tp · 1000`（m/h → mm/h） | mm/h |
| `lai` | `ds["lai"][:, :, clamp(month, 1, 12)]` | m²/m² |

### 3.4 多日滚动

```julia
day_offset  = div(hour - 1, 24)        # 0, 1, 2, ...
hour_in_day = mod1(hour, 24)           # 1..24
cur_date_dt = base_date + Day(day_offset)
cur_date    = Dates.format(cur_date_dt, dateformat"yyyymmdd")
paths_h = (cfg.mock == 0 && day_offset > 0) ?
    merge(paths, era5_paths_for(cfg.era5, cur_date)) : paths
```

真实数据模式按 `era5_paths_for($era5_root, $cur_date)` 重建当日 8 个 ERA5 路径；mock 模式复用单日数据。

### 3.5 平衡含水量（每个像元一次）

```julia
smoieq[:, i, j] .= eqsoilmoisturetheor(nzg, nsoil, z₋ₕ, dz, fdepth[i, j], wtd_in[i, j])
```

调用方式保持 Fortran `EQSOILMOISTUREtheor` 的 2-D 循环语义（每个 `(i, j)` 像元重新构造闭包），与 `test/test_eqsoilmoisture.jl` 的单元测试一致。详见 [`SoilInitialization-土壤分层.md`](./SoilInitialization-土壤分层.md) §3。

### 3.6 `rootdepth_main` 调用

与 `test/test_rootdepth.jl::make_rootdepth_inputs` 完全对齐：`freedrain = 1`、`Δt = 3600.0`、`steps = 1.0`、`veg = 7.0`（短草）、`icefactor .= 0`（暖季无冻土）。所有 47 个形参按行号顺序传入；详见 [`RootDepth-主算法.md`](./RootDepth-主算法.md) §1、§2。

## 4. 关键变量与单位

| 符号 | 含义 | 单位 |
|---|---|---|
| `nx × ny` | 模型网格大小（i 方向 × j 方向） | — |
| `nzg` | 土壤层数（默认 40） | — |
| `soiltxt[2, nx, ny]` | 土壤类型（rootdepth_main 期望 `[2, nx, ny]` 形状） | —（1..13） |
| `smoieq[nzg, nx, ny]` | 平衡含水量 | m³/m³ |
| `wtd[nx, ny]` | 地下水位深度 | m（≤ 0） |
| `Δt = 3600.0` | 时间步长 | s |
| `f.wind / temp / qair / press / netrad / rshort / precip / lai` | 当小时 ERA5 强迫（8 字段） | 见 §3.3 |
| `et_s / et_i / et_c` | 土壤 / 截留 / 冠层蒸发累计 | mm |
| `icefactor[i, j, k]` | 冻结标志（本例简化置 0） | Int8（0/1） |

## 5. 与 Fortran 对应

| Fortran 子程序 / 文件 | Julia 实现 | 差异 |
|---|---|---|
| `fortran/main.f90` | `main(ARGS)` | Fortran 是单步调度器；Julia 端把同一调度序列嵌入 `main` 并包了一层 mock 数据生成。 |
| `fortran/module_io.f90::READINITIAL`（L10） | `read_initial(paths.static)` | `src/io/NetCDF.jl` 单进程版本；`fdepth < 1e-6 → 100` 夹断 + `topo < -1e5 → 0/0` 掩码已对齐；详见 [`io-NetCDF.md`](./io-NetCDF.md) §5。 |
| `fortran/module_io.f90::READWTDNC`（L309） | `read_wtdnc(paths.wtd)` | 符号约定 `min(-raw, 0)` 一致；mock 写入正值（0.2..4.5 m）由 `read_wtdnc` 翻为非正。 |
| `fortran/module_initial.f90::EQSOILMOISTUREtheor` | `eqsoilmoisturetheor(nzg, nsoil, z₋ₕ, dz, fdepth, wtd)` | 调用方式保持 Fortran 2-D 循环语义（每个像元调一次）。 |
| `fortran/module_forcings.f90::READFORCINGS`（L69-L211） | `read_mock_hourly_forcings`（mock）/ `src/Forcings/ERA5.jl::read_hourly_forcings`（生产） | mock stub 签名与生产模块**互不替代**，详见 §7 #1。 |
| `fortran/module_rootdepth.f90::ROOTDEPTH` | `rootdepth_main(...)` | 47 个形参顺序与 `test_rootdepth.jl::make_rootdepth_inputs` 对齐。 |
| `fortran/module_io.f90::WRITEOUTPUT*` | `NCDataset(cfg.out, "c") ... end` | 仅示范 6 字段（`et_s`/`et_c`/`et_i`/`wtd`/`rech`）落盘；其余写出未实现。 |

## 6. 引用

- 行号：
  - `example/regional_example.jl:L20-L41` 文档字符串（用法 / 参数 / 风格约定）
  - `example/regional_example.jl:L51-L56` `using` 列表（`ASAP` / `NCDatasets` / `Statistics` / `Dates` / `Printf` / `Random`）
  - `example/regional_example.jl:L68-L107` `parse_cli_args`
  - `example/regional_example.jl:L141-L307` `generate_mock_dataset`（含 mock 合成公式）
  - `example/regional_example.jl:L327-L374` `read_mock_hourly_forcings`
  - `example/regional_example.jl:L390-L401` `era5_paths_for`
  - `example/regional_example.jl:L407-L645` `main`（6 步主流程）
  - `example/regional_example.jl:L413-L417` mock 输出路径 fallback
  - `example/regional_example.jl:L432-L449` 真实数据模式路径校验 + `paths` 构造
  - `example/regional_example.jl:L450-L499` 静态场初始化 + 平衡含水量 2-D 循环（`eqsoilmoisturetheor` 在 L485）
  - `example/regional_example.jl:L500-L561` 状态变量初始化（47 形参顺序）
  - `example/regional_example.jl:L564-L612` 时步循环 + 滚动日期（`read_mock_hourly_forcings` 在 L562）
  - `example/regional_example.jl:L576-L599` `rootdepth_main` 调用
  - `example/regional_example.jl:L599-L607` `@printf` 日志
  - `example/regional_example.jl:L616-L636` 最终输出 NetCDF 落盘（`NCDataset(cfg.out, "c")` 在 L618）
  - `example/regional_example.jl:L651-L653` 入口 `if abspath(PROGRAM_FILE) == @__FILE__`
- 测试断言：
  - `test/test_regional_example.jl` 「mock 模式端到端」：3×3 网格 3 时步 smoke test，断言 `cfg.out` 文件存在
  - `test/test_regional_example.jl` 「静态场往返」：mock 写 soiltxt / topo / fdepth / landmask 后 `read_initial` 还原范围一致
  - `test/test_regional_example.jl` 「wtd 符号约定」：mock 写正值，`read_wtdnc` 后 ≤ 0
  - `test/test_regional_example.jl` 「generate_mock_dataset 类型断言」：`soiltxt::Matrix{Int}`、`topo::Matrix{Float64}`、`wtd_raw::Matrix{Float64}` 等
- 文档：`README.md` 「区域应用：数据准备清单」7 小节、`CLAUDE.md` §1 / §2 / §8

## 7. 已知问题与备注

1. **`read_mock_hourly_forcings` 与 `src/Forcings/ERA5.jl::read_hourly_forcings` 互不替代**（重要）：两者签名、变量名（`wind`/`strd`/`tp` vs `ws10`/`sr`/`tp & sf`）、返回字段（8 字段 NamedTuple vs `varpack(rx, ry, 24, 11) + 10 标量`）、ERA5 文件名约定（`ERA5_wind_speed_*.nc` vs `ERA5_ws10_*.nc`）均不同。脚本顶部 docstring L113-L138 列出 4 点差异；本函数**有意识保留为独立 stub**，不依赖项目 `Forcings` 模块即可端到端跑通，未来接入完整生产流水线时应切换为 `ASAP.read_hourly_forcings`。
2. **`icefactor .= 0` 简化**：真实实现需读 `paths.era5_soilt`（4 层土壤温度 K）按层映射 `icefactor[i, j, k] = T ≤ 273.15 ? Int8(1) : Int8(0)`；本示例把冻结层全部置 0，仅覆盖暖季路径。详见 [`Forcings-ERA5.md`](./Forcings-ERA5.md) §3.7。
3. **`netrad` 粗略近似**：`netrad = ssrd·(1 − 0.23) − strd` 直接以地表短波/长波代数差替代净辐射，未做晴空衰减与地表反射分离。生产实现应通过 `read_accumulated_forcings` 取得 ERA5 `rnetera` 累积量并 ÷3600。
4. **`rshort = ssrd`**：未做晴空区分，与 `read_accumulated_forcings` 中 `swdown_h = swdown / 3600`（W/m²）等价。
5. **冻结路径仅占位**：本示例未在 `rootdepth_main` 之前调 `read_soil_temps` 填充 `icefactor`；若 `--duration` 跨过秋冬，应在循环内补 ERA5 4 层 STL 读取。
6. **多日滚动依赖文件命名约定**：真实数据模式按 `$era5_root/{VARNAME}/ERA5_{VARNAME}_{YYYYMMDD}.nc` 重建路径；与 `src/Forcings/ERA5.jl::era5_hourly_path` 的命名（`ERA5_ws10_*` / `ERA5_2m_temperature_*` 等）**不一致**，二者分别对应 README §「区域应用：数据准备清单」与 README §「完整 ERA5 模块路径」两套约定。生产接入时需统一命名。
7. **输出 NetCDF 仅示范 6 字段**：`et_s` / `et_c` / `et_i` / `wtd` / `rech`；其余诊断（`smoi`、`qsrun`、`waterdeficit`、`qlat` 等）落盘未实现，按 `CLAUDE.md §6.1` 留待后续 PR。
8. **`parse_cli_args` 字符串字段直接赋值**：`String` 字段不需要 `parse(String, raw)`，已用类型分支处理（`T === String ? raw : parse(T, raw)`）。
9. **`docstring` 字符串插值转义**：原顶部 docstring 含 `$root` 触发 Julia 字符串插值，已用 `\$root` 转义（避免 `UndefVarError: root`）。
10. **依赖**：除 `ASAP` 外还直接 `using NCDatasets / Statistics / Dates / Printf / Random`；`Project.toml` test 目标已包含这 5 个包。详见 [`CLAUDE.md` §8](../../CLAUDE.md)。