# ASAP 模型 Julia 版本

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://CUG-hydro.github.io/ASAP-model/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://CUG-hydro.github.io/ASAP-model/dev)
[![CI](https://github.com/CUG-hydro/ASAP-model/actions/workflows/CI.yml/badge.svg)](https://github.com/CUG-hydro/ASAP-model/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/CUG-hydro/ASAP-model/branch/main/graph/badge.svg)](https://app.codecov.io/gh/CUG-hydro/ASAP-model/tree/main)

这是 ASAP (Agricultural Systems Analysis and Prediction) 模型的 Julia 实现版本，从原始的 Fortran `module_rootdepth.f90` 翻译而来。

## 模块结构

Julia 版本将原始的单一 Fortran 模块分解为多个独立的、可测试的模块：

### 核心模块

1. **SoilParameters.jl** - 土壤参数模块
   - 定义13种土壤类型的参数
   - 计算土壤导水率
   - 初始化土壤参数

2. **Evapotranspiration.jl** - 蒸散发计算模块
   - Priestley-Taylor 方法
   - Penman-Monteith 方法
   - Shuttleworth-Wallace 双源法

3. **Interception.jl** - 截留模块
   - 植被截留计算
   - 截留蒸发计算

4. **WaterExtraction.jl** - 水分提取模块
   - 植物根系水分提取
   - 土壤水分胁迫计算

5. **SoilFluxes.jl** - 土壤水流模块
   - 土壤水分运动方程
   - 三对角矩阵求解器
   - 氧18同位素追踪

6. **WaterTableDynamics.jl** - 地下水位动态模块
   - 浅层地下水位更新
   - 侧向流影响

7. **SoilInitialization.jl** - 土壤初始化模块
   - CLM 方法土壤层设置
   - 固定层厚方法

8. **RootDepth.jl** - 主模块
   - 整合所有子模块
   - 提供主要计算接口

## 使用方法

```julia
# 获取土壤参数
soil_params = get_soil_params(5)  # 第5种土壤类型

# 计算蒸散发
pet = potevap_penman_monteith(1, 1, 298.15, 200.0, 300.0, 101325.0, 0.01, 2.0, 3.0, 5.0, 15.0)

# 计算截留
ppdrip, et_i, new_store = interception(0.01, 5.0, 3.0, 0.2, 0.5)
```
<!-- 
### 主要特性：

- **完整的物理过程**：保留了原始模型的所有物理计算

- **多种蒸散发方法**：支持三种不同的蒸散发计算方法

- **同位素追踪**：包含氧18同位素的追踪计算

- **灵活的土壤配置**：支持多种土壤层配置方案

- **数值稳定性**：改进的数值算法确保计算稳定 -->

---

## 区域应用：数据准备清单

本节说明在**真实数据模式**下运行区域端到端示例（`example/regional_example.jl`）所需的全部输入文件及其格式约定。Mock 模式（`--mock 1`）由 `generate_mock_dataset()` 在 `mktempdir()` 中自动生成等价数据集，无需任何外部文件，可用于验证安装及调试流程。

### 1. 静态场文件（`static.nc`）

由 `src/io/NetCDF.jl::read_initial` 读取，对应 Fortran `module_io.f90::READINITIAL`。

| NetCDF 变量 | 类型       | 维度   | 单位  | 含义                                            |
| ----------- | ---------- | ------ | ----- | ----------------------------------------------- |
| `STXT`      | Int16/Int32 | x, y  | —     | USDA 土壤质地类别，整数 1–13                    |
| `topo`      | Float64    | x, y   | m     | 地表高程；值 < −1×10⁵ 视为海洋（自动置 0）     |
| `F`         | Float64    | x, y   | m     | 根系深度因子（fdepth）；< 1×10⁻⁶ 自动替换为 100 |

**陆地掩膜**由读取函数内部派生：`landmask = 1`（`topo ≥ −1×10⁵`），`landmask = 0`（海洋）；无需单独提供。

> 来源：`src/io/NetCDF.jl` — `read_initial`（L63–L100）；
> `generate_mock_dataset` 中的写入示例见 `example/regional_example.jl`（L133–L141）。

### 2. 初始地下水位文件（`wtd.nc`）

由 `src/io/NetCDF.jl::read_wtdnc` 读取，对应 Fortran `module_io.f90::READWTDNC`（L309）。

| NetCDF 变量 | 类型    | 维度 | 单位 | 含义                                           |
| ----------- | ------- | ---- | ---- | ---------------------------------------------- |
| `WTD`       | Float64 | x, y | m    | 地下水位深度（**写入正值**，单位 m，地面以下） |

**符号约定（重要）**：文件中存储**正值**（深度），读取函数内部执行 `wtd = min(−raw, 0.0)`，返回值为 ≤ 0 的负数（地面以下为负，地面水位为 0）。例如，写入 `1.5` 表示地下水位在地面以下 1.5 m，读出后为 `−1.5`。此约定与 Fortran 原版一致。

> 来源：`src/io/NetCDF.jl` — `read_wtdnc`（L113–L126）；
> `example/regional_example.jl`（L144–L149）。

### 3. ERA5 小时强迫文件（真实数据模式）

ERA5 强迫数据按**变量分目录**存储，每日一文件。给定 ERA5 根目录 `$ERA5_ROOT` 和日期 `YYYYMMDD`，文件路径模板为：

```
$ERA5_ROOT/
├── WIND/       ERA5_wind_speed_YYYYMMDD.nc
├── TEMP/       ERA5_2m_temperature_YYYYMMDD.nc
├── DEWPOINT/   ERA5_2m_dewpoint_YYYYMMDD.nc
├── SFCPRESS/   ERA5_surface_pressure_YYYYMMDD.nc
├── STRD/       ERA5_strd_YYYYMMDD.nc
├── SSRD/       ERA5_ssrd_YYYYMMDD.nc
├── SOILT/      ERA5_soil_temps_YYYYMMDD.nc
└── TP/         ERA5_total_precipitation_YYYYMMDD.nc
```

各文件的变量规格如下（变量名即代码中 `ds[varname]` 读取的 NetCDF 变量名）：

| 子目录      | 文件名模板                              | NetCDF 变量 | 维度          | 单位   | 含义                          |
| ----------- | --------------------------------------- | ----------- | ------------- | ------ | ----------------------------- |
| `WIND`      | `ERA5_wind_speed_YYYYMMDD.nc`           | `wind`      | x, y, hour   | m/s    | 10 m 风速                     |
| `TEMP`      | `ERA5_2m_temperature_YYYYMMDD.nc`       | `t2m`       | x, y, hour   | K      | 2 m 气温                      |
| `DEWPOINT`  | `ERA5_2m_dewpoint_YYYYMMDD.nc`          | `d2m`       | x, y, hour   | K      | 2 m 露点温度                  |
| `SFCPRESS`  | `ERA5_surface_pressure_YYYYMMDD.nc`     | `sp`        | x, y, hour   | Pa     | 地面气压                      |
| `STRD`      | `ERA5_strd_YYYYMMDD.nc`                 | `strd`      | x, y, hour   | W/m²   | 向下长波辐射                  |
| `SSRD`      | `ERA5_ssrd_YYYYMMDD.nc`                 | `ssrd`      | x, y, hour   | W/m²   | 向下短波辐射                  |
| `SOILT`     | `ERA5_soil_temps_YYYYMMDD.nc`           | `STL1`      | x, y, hour, soil_level | K | 四层土壤温度（层 1–4） |
| `TP`        | `ERA5_total_precipitation_YYYYMMDD.nc`  | `tp`        | x, y, hour   | m/h    | 总降水（小时累计量）          |

**维度说明**：`hour` 维长度为 24（UTC 0–23 时），`soil_level` 维长度为 4（ERA5 四层：0–7 cm、7–28 cm、28–100 cm、100–289 cm）。`x`、`y` 为模型区域网格大小。

> 来源：`example/regional_example.jl` — `generate_mock_dataset`（L152–L198）及 `main` 真实数据路径块（L254–L271）；`read_hourly_forcings`（L217–L250）。

### 4. LAI 月度气候态文件

| 文件路径（相对 `static.nc` 父目录的上级）| NetCDF 变量 | 维度           | 单位   | 含义                       |
| ---------------------------------------- | ----------- | -------------- | ------ | -------------------------- |
| `../LAI/lai_clim_01.nc`                  | `lai`       | x, y, month   | m²/m²  | 月度 LAI 气候态（12 个月） |

月份索引为 1–12，代码按当月取对应切片。真实数据模式中路径约定为 `dirname(static) + "/../LAI/lai_clim_01.nc"`。

> 来源：`example/regional_example.jl`（L273、L238–L248）。

### 5. 模型内部计算的派生量

以下变量由代码从以上原始输入计算，**用户无需提供**：

| 派生量    | 计算来源                   | 公式/方法                                                          | 单位   |
| --------- | -------------------------- | ------------------------------------------------------------------ | ------ |
| `qair`    | `d2m`（K）+ `sp`（Pa）    | Tetens 公式：`e = 6.112·exp(17.67·(Td−273.15)/(Td−29.65))·100 Pa`；`q = 0.622·e/(sp−0.378·e)` | kg/kg |
| `netrad`  | `ssrd`（W/m²）+ `strd`（W/m²） | 净辐射近似：`ssrd·(1−α) − strd`，短草 `α = 0.23`              | W/m²  |
| `precip`  | `tp`（m/h）                | 单位换算：`tp × 1000`                                             | mm/h  |
| `landmask`| `topo`                     | `topo < −1×10⁵` → 0（海洋），其余 → 1（陆地）                    | —      |

> 来源：`example/regional_example.jl` — `dewpoint_to_qair`（L219–L229）、`netrad`（L232–L233）、`precip`（L236）；`src/io/NetCDF.jl`（L82–L88）。

---

### 6. 完整 ERA5 模块路径（`src/Forcings/ERA5.jl`，供未来区域运行参考）

`src/Forcings/ERA5.jl` 实现了对 Fortran `module_forcings.f90` 中 6 个核心强迫读取子程序的完整 Julia 翻译，采用**不同于示例**的文件路径约定，适用于处理原始全球 ERA5 数据的完整区域流水线。

**路径模板**：`$ROOT/{VARNAME}/ERA5_{VARNAME}_{YYYYMMDD}.nc`，目录名使用以下大写短名：

| Fortran 子程序目录名 | NetCDF 变量 | 维度                   | 单位   | 含义                             |
| -------------------- | ----------- | ---------------------- | ------ | -------------------------------- |
| `WIND`               | `ws10`      | 1440, 721, 24          | m/s    | 10 m 风速（全球 0.25° 格点）     |
| `TEMP`               | `t2m`       | 1440, 721, 24          | K      | 2 m 气温                         |
| `SFCPRESS`           | `sp`        | 1440, 721, 24          | Pa     | 地面气压                         |
| `DEWPOINT`           | `d2m`       | 1440, 721, 24          | K      | 2 m 露点温度（→ RH/qair）        |
| `PRECIP`             | `tp`, `sf`  | 1440, 721, 24          | m      | 总降水、降雪（同文件；差值取液态降水，×1000 转 mm） |
| `RADIATION`          | `ssrd`, `sr` | 1440, 721, 24         | J/m²   | 向下短波、净辐射（÷3600 转 W/m²）|
| `SOILTEMP`           | `stl1`–`stl4` | 1440, 721, 24       | K      | 四层土壤温度                     |

**全球-区域裁切**：Fortran 代码读取全球格点（1440×721），随后裁切 `varread(1067:1316, 298:587, :)` 得到 250×290 的区域子域，存入 `varpack(250, 290, 24, 11)` 数组（11 个变量通道）。Julia 端执行相同裁切（见 `src/Forcings/ERA5.jl::_read_era5_3d`）。

此外，`read_topo_era5` 读取 ERA5 静态地势高（变量 `z`，单位 m²/s²，÷9.81 转 m），用于土壤温度的高度订正（lapse rate −0.0065 K/m）；`read_lai_climatology` 读取 LAI 气候态（变量 `LAI`，维度 x, y, month=12，单位 m²/m²）。

**当前阶段范围外项目**（项目政策明确禁止，见 `CLAUDE.md §6.1`）：

- NetCDF 写出（`WRITEOUTPUTNC*` / `WRITEHISTORYNC*`）及其余读取子程序（`READLATLON` 等）
- MPI 并行化（`MPI.jl` / 域分解）
- ¹⁸O 同位素追踪串联（`lateral_isotope!` / `updatedeepwtable!`）

> 来源：`src/Forcings/ERA5.jl` — 模块文档（L1–L44）、`ERA5_GLOBAL_NX/NY`（L56–L62）、`ERA5_REGION_NX/NY`（L68–L76）、`era5_hourly_path`（L95–L109）、`_read_era5_3d`（L116–L130）、`read_hourly_forcings`（L132–L250）、`read_topo_era5`（L370–L430）、`read_lai_climatology`（L436–L480）、`compute_qair`（L488–L510）。

---

### 7. 运行命令速查

#### Mock 模式（默认，无需任何外部文件）

```bash
# 在 mktempdir() 自动生成合成数据，运行 3 小时模拟
julia --project example/regional_example.jl --mock 1 --duration 3

# 指定网格大小和土壤层数
julia --project example/regional_example.jl --mock 1 --nx 10 --ny 10 --nzg 40 --duration 24
```

#### 真实数据模式

```bash
julia --project example/regional_example.jl \
      --mock   0 \
      --date   20200101 \
      --duration 24 \
      --static /data/ASAP/static.nc \
      --wtd    /data/ASAP/wtd.nc \
      --era5   /data/ERA5 \
      --out    /tmp/asap_regional_20200101.nc
```

**命令行参数一览**：

| 参数          | 默认值       | 含义                                                  |
| ------------- | ------------ | ----------------------------------------------------- |
| `--date`      | `20200101`   | 起始日期（YYYYMMDD），决定 ERA5 文件日期               |
| `--duration`  | `3`          | 模拟小时数                                            |
| `--nx`        | `3`          | 网格 x 方向格点数（mock 模式）                        |
| `--ny`        | `3`          | 网格 y 方向格点数（mock 模式）                        |
| `--nzg`       | `40`         | 土壤层数                                              |
| `--mock`      | `1`          | `1` = 合成数据（无需外部文件），`0` = 读真实 NetCDF   |
| `--static`    | —            | 静态场 NetCDF 路径（真实模式必填）                    |
| `--wtd`       | —            | 初始水位 NetCDF 路径（真实模式必填）                  |
| `--era5`      | —            | ERA5 数据根目录（真实模式必填）                       |
| `--out`       | `/tmp/asap_regional_<date>.nc` | 输出 NetCDF 路径                    |

> 参数来源：`example/regional_example.jl` — 文件头部文档字符串（L14–L44）及 `parse_cli_args`（L53–L100）。
