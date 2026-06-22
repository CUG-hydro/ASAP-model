> 源文件：`src/Forcings/ERA5.jl`
> Fortran 来源：`fortran/module_forcings.f90`（子集：`READFORCINGS`、`READFORCINGSACC`、`READFORCINGSSNOW`、`READFORCINGSSOILT`、`READTOPOERA5`、`READLAICLIM`）
> 测试：`test/test_forcings_era5.jl`
> 状态：已摄取（P1-A 里程碑；`READO18CLIM` / 全量 `READFORCINGS` 重投影到模型网格 / NetCDF 写出均不在本期范围，见 §6.1 与 §7）

## 1. 功能概述

`module ERA5Forcings` 是 ASAP-model 的 ERA5 强迫**读取**入口，对应 Fortran `module_forcings.f90` 中 6 个最高优先级的子程序。它把全球 0.25° ERA5 / 0.1° ERA5-Land 网格上的小时尺度气象场（风、气温、气压、露点、降水、短波/净辐射、4 层土壤温度、雪融、地形、LAI 月气候态）裁剪到生产子区域 `[1067:1316, 298:587]`，并通过共享的双线性插值 helper 投影到 ASAP 的模型网格。本期**未**实现 `READO18CLIM`、把全部 ERA5 强迫重投影到模型网格的完整 `READFORCINGS` 流程、以及任何 NetCDF 写出路径。

`src/ASAP.jl:L31-L37` 通过 `include("Forcings/ERA5.jl")` 引入，并以 `using .ERA5Forcings` + 显式 `export` 的方式把 7 个对外函数（`read_hourly_forcings` / `read_accumulated_forcings` / `read_snow` / `read_snow_hour` / `read_soil_temps` / `read_topo_era5` / `read_lai_climatology` / `compute_qair`）平铺到顶层命名空间。

## 2. 函数签名

```julia
# 网格抽象（new in P1-A）
struct Grid
    name::String
    origin::Tuple{Float64,Float64}      # (xswlon°, xswlat°)
    res::Tuple{Float64,Float64}         # (dx°, dy°)
    shape::Tuple{Int,Int}               # (nx_src, ny_src)
end

const ERA5_GRID     = Grid("ERA5",      (266.5, -56.5), (0.25, 0.25), (1440, 721))
const ERA5LAND_GRID = Grid("ERA5-Land", (266.5, -56.5), (0.10, 0.10), (3600, 1801))

# 1. 24 h 小时强迫打包 → varpack(rx, ry, 24, 11)
function read_hourly_forcings(date::String, root::String;
                              lats::AbstractMatrix, lons::AbstractMatrix,
                              landmask::AbstractMatrix,
                              xrange::AbstractUnitRange = ERA5_REGION_XRANGE,
                              yrange::AbstractUnitRange = ERA5_REGION_YRANGE)

# 2. 把 varpack 中累积量双线性插值到模型网格（含 J/m² → W/m²）
function read_accumulated_forcings(hour::Int, varpack::Array{Float64,4};
                                   lats::AbstractMatrix, lons::AbstractMatrix,
                                   landmask::AbstractMatrix,
                                   grid::Grid = ERA5_GRID)

# 3. ERA5-Land 单小时雪融（smlt, m → mm, clamp ≥ 0）
function read_snow_hour(date::String, hour::Int, root::String;
                        lats::AbstractMatrix, lons::AbstractMatrix,
                        landmask::AbstractMatrix,
                        grid::Grid = ERA5LAND_GRID,
                        xrange::AbstractUnitRange = ERA5LAND_SNOW_XRANGE,
                        yrange::AbstractUnitRange = ERA5LAND_SNOW_YRANGE)

read_snow(date::String, root::String; kwargs...) :: Matrix{Float64}  # 一行 alias：read_snow_hour(date, 1, root; kwargs...)

# 4. 4 层土壤温度 → icefactor(nx, ny, 15) Int8 冻结标志
function read_soil_temps(hour::Int, varpack::Array{Float64,4};
                         nzg::Int,
                         lats::AbstractMatrix, lons::AbstractMatrix,
                         landmask::AbstractMatrix,
                         topo::AbstractMatrix, topoera::AbstractMatrix,
                         grid::Grid = ERA5_GRID)

# 5. 静态地形（z, m²/s² ÷ 9.81, north-up flip）
function read_topo_era5(path::String;
                        lats::AbstractMatrix, lons::AbstractMatrix,
                        landmask::AbstractMatrix,
                        grid::Grid = ERA5_GRID,
                        scale_factor::Real = 1.0,
                        add_offset::Real   = 0.0)

# 6. LAI 月气候态单月切片（要求源已在模型网格上）
function read_lai_climatology(path::String, month::Int;
                              target_shape::Tuple{Int,Int})

# 7. 由 d2m(K) 与 sp(Pa) 算比湿 qair (kg/kg)
compute_qair(d2m_K::AbstractArray, sp_Pa::AbstractArray) :: Array{Float64}
```

`varpack` 11 个槽位的排列顺序（与 Fortran `module_forcings.f90` L90-L210 一致）：

| 索引 | 含义 |
|---|---|
| 1 | `wind` (ws10, m/s) |
| 2 | `tempera` (t2m, K) |
| 3 | `press` (sp, Pa) |
| 4 | `rh` (相对湿度，0..1) |
| 5 | `rain` (m → mm，clamp ≥ 0) |
| 6 | `swdown` (J/m² 小时累积) |
| 7 | `rnetera` (J/m² 小时累积) |
| 8..11 | `st1..st4` (4 层土壤温度，K) |

## 3. 算法 / 公式

### 3.1 模型网格 → 源网格分数索引（`_grid_index`）

设源 SW 角为 `(xswlon, xswlat)`、分辨率为 `(dx, dy)`、源 y 维长度 `ny_src`，模型点 `(lat, lon)` 先做经度 wrap（`lon < 0 → lon + 360`），再按下式得到 1-based 分数索引：

$$
g_x = \frac{\mathrm{lon} - x_\text{swlon}}{dx} + 1,\qquad
g_y = n_{y,\text{src}} + \frac{x_\text{swlat} - \mathrm{lat}}{dy}
$$

负纬度走「北向上」约定（与 Fortran `gry = ny + (xswlat - lat)/dy` 等价）。

### 3.2 共享双线性插值（`bilinear`）

```text
rx, ry = size(src_used)         # flip_lat=true 时 src 沿 dim=2 翻转
i1 = clamp(floor(gx), 1, rx-1)  # 1-based；i2 = i1 + 1
j1 = clamp(floor(gy), 1, ry-1)  # 1-based；j2 = j1 + 1
fx = clamp(gx - i1, 0, 1)
fy = clamp(gy - j1, 0, 1)
out[i, j] = (1-fx)(1-fy) src[i1, j1] + fx(1-fy) src[i2, j1]
                + (1-fx)fy   src[i1, j2] + fx·fy    src[i2, j2]
landmask[i, j] == 0 → out[i, j] = 0   # 海洋点强制 0
```

P1-A 把原本散落在 5 个函数里的内联双线性代码折叠为这一个 helper，并通过 `flip_lat` 关键字兼容 Fortran L839-L843 的 north-up 翻转约定（仅 `read_topo_era5` 不需要再翻转，由 `flip_lat=false` 显式声明）。

### 3.3 降水拆分（`read_hourly_forcings` 内部）

$$
\mathrm{rain}_{[m]} = \max\!\big(\mathrm{tp}_{[m]} - \mathrm{sf}_{[m]},\; 0\big),\qquad
\mathrm{rain}_{[mm]} = \mathrm{rain}_{[m]} \times 1000
$$

ERA5 `tp`（total precipitation）与 `sf`（snowfall）单位均为 m；最终以 mm 进入 `varpack[:, :, :, 5]`。

### 3.4 相对湿度（Tetens 简化）

$$
e_\text{sat}(T) = 6.112 \cdot \exp\!\left(\frac{17.67\,(T - 273.15)}{T - 29.65}\right),\quad
\mathrm{rh} = \min\!\left(\frac{e_\text{sat}(T_d)}{e_\text{sat}(T)},\; 1.0\right)
$$

其中 `e_sat` 在内部以 hPa 计算；最终 `rh ∈ [0, 1]` 写入 `varpack[:, :, :, 4]`。

### 3.5 辐射通量（`read_accumulated_forcings`）

ERA5 `ssrd` / `sr` 为小时累积（J/m²）；按 `÷ 3600` 转为 W/m² 再写入 `rain_h / swdown_h / rnetera_h` 三个 `view` 后双线性插值到模型网格。

### 3.6 比湿（`compute_qair`，Tetens）

$$
e = 6.112 \cdot \exp\!\left(\frac{17.67\,(T_d - 273.15)}{T_d - 29.65}\right) \times 100\;\text{Pa},\qquad
q = \frac{0.622\,e}{\mathrm{sp} - 0.378\,e}\quad [\text{kg/kg}]
$$

输入 `d2m_K`、`sp_Pa` 同形状数组；逐元素运算，返回同形状 `Array{Float64}`。

### 3.7 土壤温度 → 冻结因子（`read_soil_temps`）

先对 4 层 `stl1..stl4` 各自双线性插值，再用 ERA5 标准温度直减率订正到模型地形：

$$
T_\text{adj} = T_\text{era} + (-0.0065)\,\big(z_\text{topo} - z_\text{topoera}\big)
$$

（`topofactor = -0.0065 * (topo - topoera)` 在 `icefactor` 内部使用；Fortran 单位 K/m。）

冻结标志填入 `icefactor[i, j, k] = T_adj ≤ 273.15 ? Int8(1) : Int8(0)`；最后一维长度为 15（对应 Fortran `icefactor(:, :, nzg-14:nzg)`）。

### 3.8 地形几何高（`read_topo_era5`）

```
zraw = ds["z"][:, :]                    # geopotential (m²/s²)
topoera[:, ny-j+1] = zraw[:, j] * scale + offset   # north-up flip
topoera /= 9.81                         # → geometric height (m)
return bilinear(topoera, ERA5_GRID, lats, lons; landmask=..., flip_lat=false)
```

`scale_factor` / `add_offset` 关键字用于与某些带 `scale_factor`/`add_offset` NetCDF 属性的静态文件对接（默认 1.0 / 0.0）。

## 4. 关键变量与单位

| 符号 | 含义 | 单位 |
|---|---|---|
| `Grid.origin` | 源网格 SW 角 `(xswlon°, xswlat°)` | °E / °N |
| `Grid.res` | 源网格分辨率 `(dx°, dy°)` | ° |
| `Grid.shape` | 源网格维度 `(nx, ny)` | - |
| `date` | yyyymmdd 字符串 | - |
| `hour` | 0..23 小时索引（1-based） | - |
| `root` | ERA5 文件根目录 | - |
| `xrange` / `yrange` | 区域裁剪范围（1-based inclusive） | - |
| `varpack[:, :, :, 1..11]` | 见 §2 表 | 见 §2 表 |
| `wind`, `tempera`, `press`, `dtemp` | 10 m 风 / 2 m 气温 / 地面气压 / 2 m 露点 | m/s / K / Pa / K |
| `tp`, `sf`, `rain` | 总降水 / 雪 / 雨（`max(tp-sf,0)*1000`） | m / m / mm |
| `swdown`, `rnetera` | 短波 / 净辐射（ERA5 小时累积） | J/m² |
| `st1..st4` | ERA5 4 层土壤温度 | K |
| `smlt` | ERA5-Land 雪融水 | m（输出 mm） |
| `icefactor` | 冻结标志 `Int8`（0 液态 / 1 冻结） | - |
| `topo`, `topoera` | 模型地形 / ERA5 地形 | m |
| `LAI` | 月气候态叶面积指数（源已在模型网格上） | - |
| `d2m_K`, `sp_Pa` | 比湿输入（露点 K、气压 Pa） | K / Pa |
| `qair` | 比湿输出 | kg/kg |

模块常量：`HOURS_PER_DAY = 24`、`N_SOIL_LEVELS = 4`、`T_FREEZE_K = 273.15`、`ERA5_REGION_XRANGE = 1067:1316`、`ERA5_REGION_YRANGE = 298:587`、`ERA5LAND_SNOW_XRANGE = 2666:3285`、`ERA5LAND_SNOW_YRANGE = 747:1466`。

## 5. 与 Fortran 对应

| Fortran 子程序 | Julia 函数 | 差异 |
|---|---|---|
| `READFORCINGS`（L69-L211） | `read_hourly_forcings` | Fortran 串行 + MPI 分发 + 同步写私有的 `varpack(iw,je,24)`；Julia 单进程 + `_crop3d` 一次性读取 11 个变量并打包。区域范围由硬编码 `varread(1067:1316, 298:587, :)` 改为 `xrange`/`yrange` kwarg，默认值保持不变。 |
| `READFORCINGSACC`（L211-L318） | `read_accumulated_forcings` | Fortran 在每个模型格点上调一次双线性；Julia 改用共享 `bilinear` helper 并把辐射 J/m² → W/m² 转换集中于此函数。 |
| `READFORCINGSSNOW`（L612-L700） | `read_snow_hour` / `read_snow` | Fortran 固定 0.1° / 2666:3285 × 747:1466 区域；Julia 通过 `grid` + `xrange`/`yrange` kwarg 参数化，新增 `read_snow` 一行 alias 兼容老调用。 |
| `READFORCINGSSOILT`（L700-L820） | `read_soil_temps` | Fortran 4 层独立循环内联双线性 + 直减率订正；Julia 折叠为 1 次共享 helper + 1 个 `topofactor` 广播。返回 `Int8` 数组形状 `(nx, ny, 15)` 与 Fortran `icefactor(:, :, nzg-14:nzg)` 对齐。 |
| `READTOPOERA5`（L820-L870） | `read_topo_era5` | Fortran `topoera(:, ydim-j+1) = zraw(:, j)` 内联；Julia 用 helper `flip_lat` 兼容，并把 `÷ 9.81` 提到循环外。新增 `scale_factor` / `add_offset` kwarg 适配属性化文件。 |
| `READLAICLIM`（L870-L920） | `read_lai_climatology` | Fortran 6 个位置形参 `(n2, n3, nw, ne, ns, nn)` 决定裁剪；Julia 简化为单 kwarg `target_shape::Tuple{Int,Int}`，源尺寸不匹配时抛 `DimensionMismatch`。 |
| `compute_qair`（无对应 Fortran） | `compute_qair` | Tetens 公式与 Fortran `compute_qair.f90`（不在 `module_forcings.f90` 中）一致；为 Julia 端独立的便捷 helper。 |
| `READO18CLIM`（P1-B） | — | **未实现**：同位素强迫留待后续 PR；参见 `wiki/julia/IsotopeTracing-同位素追踪.md`。 |
| 完整 `READFORCINGS` 重投影到模型网格 | — | **未实现**：当前只暴露 24 h 小时打包，模型网格上的时序插值由 `read_accumulated_forcings` 接管（按小时）。 |
| 全部 NetCDF 写出 | — | **未实现**：与 `io-NetCDF.md` P0 子集一致，本模块只读不写。 |

## 6. 引用

- 行号：
  - `src/Forcings/ERA5.jl:L17-L27` `Grid` struct 与常量定义
  - `src/Forcings/ERA5.jl:L30-L34` 区域范围常量（`ERA5_REGION_*`、`ERA5LAND_SNOW_*`）
  - `src/Forcings/ERA5.jl:L38-L41` `era5_hourly_path`
  - `src/Forcings/ERA5.jl:L45-L53` `_crop3d`
  - `src/Forcings/ERA5.jl:L58-L70` `_grid_index`
  - `src/Forcings/ERA5.jl:L78-L117` `bilinear`（共享双线性 helper）
  - `src/Forcings/ERA5.jl:L120-L170` `read_hourly_forcings`
  - `src/Forcings/ERA5.jl:L175-L195` `read_accumulated_forcings`
  - `src/Forcings/ERA5.jl:L200-L235` `read_snow_hour` / `read_snow`
  - `src/Forcings/ERA5.jl:L240-L290` `read_soil_temps`
  - `src/Forcings/ERA5.jl:L295-L362` `read_topo_era5`
  - `src/Forcings/ERA5.jl:L367-L385` `read_lai_climatology`
  - `src/Forcings/ERA5.jl:L388-L395` `compute_qair`
  - `src/ASAP.jl:L31-L37` `include("Forcings/ERA5.jl")` + 显式 export
- 测试断言：
  - `test/test_forcings_era5.jl` 「`compute_qair` analytic values」: 10°C/1013.25 hPa 下 qair ≈ 0.0075（公差 5e-4）
  - `test/test_forcings_era5.jl` 「`read_hourly_forcings`」: 11 个 `varpack` 槽位的形状、单位、值与 `(tp-sf)*1000` 转换
  - `test/test_forcings_era5.jl` 「`read_accumulated_forcings`」: varpack → 模型网格的形状与单位转换
  - `test/test_forcings_era5.jl` 「`read_soil_temps (ice factor)`」: 4 层冻结 / 液态切换下的 `icefactor` 取值
  - `test/test_forcings_era5.jl` 「`read_snow`」: smlt × 1000 mm 与 clamp ≥ 0
  - `test/test_forcings_era5.jl` 「`read_topo_era5 (static)`」: `z=9810 m²/s² → 1000 m` 与 north-up flip 兼容
  - `test/test_forcings_era5.jl` 「`read_lai_climatology`」: `target_shape` 一致性、`DimensionMismatch` 错误
- 文档：`docs/*.typ` 暂未涉及（强迫读取若需推导可后续补 Typst）。

## 7. 已知问题与备注

1. **`Grid` 抽象取代硬编码 SW 角（new in P1-A）**：原 Fortran 中 `xswlat = -56.5`、`xswlon = 266.5`、`dxv = dyv = 0.25` 等常量被收纳进 `ERA5_GRID` / `ERA5LAND_GRID`。所有 `read_*` 函数都通过 `grid::Grid = ERA5_GRID` 关键字接收源网格描述，调用方可传自定义 `Grid` 适配研究区域。
2. **测试用 mock 小网格**：测试套件构造 `(30, 30)` 全球网格 + 区域 `(5:14, 5:14)`（10×10），通过 `xrange`/`yrange` kwarg 传入，使每个 mock 文件 < 1 MB；生产路径不传 kwarg 即可保持原 Fortran `[1067:1316, 298:587]` 默认行为。
3. **共享双线性 helper 替换 5 处内联**：原 Fortran 中分散在 5 个子程序里的双线性代码统一为 `bilinear`，仅 `read_topo_era5` 通过 `flip_lat=false` 显式声明无需再翻转（`flip_lat=true` 兼容部分遗留静态文件的 north-up 翻转约定）。`landmask==0` 的格点输出强制为 0。
4. **`read_snow` 是 `read_snow_hour(date, 1, root; ...)` 的 alias**：单行定义，完整关键字透传；上层调用建议直接用 `read_snow_hour` 明确小时。
5. **`read_lai_climatology` 简化签名**：从 6 个位置形参 `(n2, n3, nw, ne, ns, nn)` 改为单 kwarg `target_shape::Tuple{Int,Int}`，并以 `DimensionMismatch` 显式报错（替代 Fortran 中静默越界）。
6. **`compute_qair` Tetens 公式保持不变**：输入/输出均为同形状数组（`AbstractArray`），逐元素广播；适用于 `Matrix{Float64}`、`Vector{Float64}` 以及更高维数组。
7. **未实现范围（与 §6.1 一致）**：
   - `READO18CLIM`（同位素强迫读取，留待 P1-B / `IsotopeTracing` 串联后启用）；
   - 完整的 `READFORCINGS` 把全部 ERA5 强迫重投影到模型网格（当前只暴露 `read_hourly_forcings` 24 h 小时打包 + `read_accumulated_forcings` 单小时重投影）；
   - 任何 NetCDF 写出（`WRITEOUTPUTNC*` / `WRITEHISTORYNC*`）— 与 `io-NetCDF.md` P0 子集一致。
8. **`scale_factor` / `add_offset` 仅 `read_topo_era5` 暴露**：其它 reader 假设源数据已是物理量，需要时再按同样模式扩展。
9. **MPI 行为省略**：与 `io-NetCDF.md` §7 同，所有读者均为单进程、整域读取；原 Fortran 的 MPI 分发与子区域裁剪已删除。