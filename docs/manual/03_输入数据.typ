// docs/manual/03_输入数据.typ
// ==========================
// 章：输入数据准备
// 来源：src/io/NetCDF.jl, example/regional_example.jl,
//       src/Forcings/ERA5.jl, src/helper.jl,
//       fortran/module_io.f90, fortran/module_driver.f90
#import "@preview/modern-cug-report:0.1.3": *
#show: doc => template(doc, footer: "CUG水文气象学2026", header: "")

#import "../preamble.typ": *
#set par(leading: 1em, spacing: 1.24em, first-line-indent: 2em)

= 1 输入数据准备

本章列出运行 `example/regional_example.jl` 真实模式（`--mock 0`）所需的
全部输入。Mock 模式（`--mock 1`）由 `generate_mock_dataset()` 在
`mktempdir()` 自动合成，无需任何外部文件。

== 1.1 概览

模型启动时按 Fortran `module_driver.f90` 的驱动顺序加载下列文件。Julia
端当前已落地的子集标记为 ✓，其余为预留给未来 P1/P2 模块的输入契约。

#figure(
  caption: [输入文件总览],
  align(center)[
  #three-line-table()[
    | *文件*                  | *字段*                                            | *状态* |
    | ---                     | ---                                               | ---    |
    | `static.nc`             | STXT, topo, F（派生 landmask）                    | ✓      |
    | `wtd.nc`                | WTD（初始地下水位）                               | ✓      |
    | `veg.nc` / `hveg.nc`    | 植被类型、植被高度                                | P1     |
    | `smoieq.nc`             | 平衡含水量（现由 `eqsoilmoisturetheor` 计算）     | P1     |
    | `riverparameters.dat`   | 流向、河长、河宽、坡度、maxdepth 等（见 §静态）   | P2     |
    | `lats` / `lons` / `area`| 网格经纬度、面积、格距                            | P2     |
    | ERA5 × 8 通道           | 风、气温、湿度、气压、辐射、降水、土壤温度        | ✓      |
    | `lai_clim_01.nc`        | LAI 月度气候态                                    | ✓      |
  ]
])

== 1.2 静态场

=== 1.2.1 地形与土壤（`static.nc`）

由 `src/io/NetCDF.jl::read_initial` 读取（对应 Fortran
`module_io.f90::READINITIAL` L10）。

#figure(
  caption: [静态场变量规格（`static.nc`）],
  align(center)[
  #three-line-table()[
    | *NetCDF 变量* | *类型*        | *单位* | *含义*                                       |
    | ---           | ---           | ---    | ---                                          |
    | `STXT`        | Int16 / Int32 | —      | USDA 土壤质地类别，整数 1–13                 |
    | `topo`        | Float64       | m      | 地表高程；值 $< -1 times 10^5$ 视为海洋       |
    | `F`           | Float64       | m      | 根系深度因子（fdepth）                       |
  ]
])

派生量在 `read_initial` 内部计算：

- `landmask`：`topo < -1e5` ⇒ 0（海洋），否则 ⇒ 1（陆地）
- `topo`：海洋格点置 0.0
- `fdepth`：`F < 1e-6` 替换为 100.0 m

=== 1.2.2 初始地下水位（`wtd.nc`）

由 `src/io/NetCDF.jl::read_wtdnc` 读取。文件中 *写入正值*（深度，单位 m），
读取后翻为 $lt.eq 0$ 负值（地面以下为负）。例如写入 `1.5` 表示地下水位
在地面以下 1.5 m，读出为 `-1.5`。此约定与 Fortran 原版一致。

=== 1.2.3 流向（D8）

`fd` / `bfd` 是河流路由与洪泛漫流的核心输入，遵循 D8 八方向编码
（`src/helper.jl::flowdir` L139–L162）：

#figure(
  caption: [D8 流向编码],
  align(center)[
  #three-line-table()[
    | *编码* | *方向*   | *i 偏移* | *j 偏移* | *几何* |
    | ---    | ---      | ---      | ---      | ---    |
    | 1      | 东       | +1       |  0       | →      |
    | 2      | 东南     | +1       | -1       | ↘      |
    | 4      | 南       |  0       | -1       | ↓      |
    | 8      | 西南     | -1       | -1       | ↙      |
    | 16     | 西       | -1       |  0       | ←      |
    | 32     | 西北     | -1       | +1       | ↖      |
    | 64     | 北       |  0       | +1       | ↑      |
    | 128    | 东北     | +1       | +1       | ↗      |
  ]
])

`0` 表示汇水点（pit）或边界；其他未列编码视为无效。`bfd`（backward flow
direction）使用相同编码，指向每个格点的上游邻居。

=== 1.2.4 河流子系统参数（`riverparameters.dat`）

Fortran `module_io.f90::READFLOWDIRECTION`（L591）+ `READRIVERPARAMETERS`
（L693）从单个 big-endian 二进制直接访问文件读取多个变量，每个变量
对应一个 record（4 字节浮点，网格大小 $n_2 times n_3$）：

#figure(
  caption: [`riverparameters.dat` 记录索引],
  align(center)[
  #three-line-table()[
    | *Record* | *变量*         | *单位* | *含义*                          |
    | ---      | ---            | ---    | ---                             |
    | 1        | `fd`           | —      | 流向（D8）                   |
    | 2        | `riverlength`  | m      | 河道长度                         |
    | 3        | `riverdepth`   | m      | 河道水深（仅 restart=1）         |
    | 4        | `riverwidth`   | m      | 河道宽度                         |
    | 5        | `slope`        | —      | 河道坡度（m/m）                  |
    | 6        | `riverflow`    | m³/s   | 河道流量（仅 restart=1）         |
    | 9        | `topoflood`    | m      | 洪泛区地形                       |
    | 10       | `bfd`          | —      | 反向流向（D8）                |
    | 11       | `maxdepth`     | m      | 河道最大深度                     |
  ]
])

派生量在 Fortran 端即时计算：`riverarea = riverwidth * riverlength`、
`floodarea = max(area - riverarea, 0.)`、`riverchannel = maxdepth * riverarea`。
Julia 端 `rivers_kw_flood!` 等函数已实现（`src/modules/`），但读取该
文件的 P2 任务尚未落地。

== 1.3 时变强迫

=== 1.3.1 ERA5 小时强迫

`example/regional_example.jl::era5_paths_for` 按「变量分目录、每日一文件」
约定构造路径：

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

各文件均含 `x, y, hour` 三维（`SOILT` 额外含 `soil_level` 维，长度 4），
`hour = 0..23` 对应 UTC 0–23 时。

#figure(
  caption: [ERA5 强迫变量规格],
  align(center)[
  #three-line-table()[
    | *子目录*    | *变量*             | *单位* | *含义*                          |
    | ---         | ---                | ---    | ---                             |
    | `WIND`      | `wind`             | m/s    | 10 m 风速                       |
    | `TEMP`      | `t2m`              | K      | 2 m 气温                        |
    | `DEWPOINT`  | `d2m`              | K      | 2 m 露点温度                    |
    | `SFCPRESS`  | `sp`               | Pa     | 地面气压                        |
    | `STRD`      | `strd`             | W/m²   | 向下长波辐射                    |
    | `SSRD`      | `ssrd`             | W/m²   | 向下短波辐射                    |
    | `SOILT`     | `STL1`..`STL4`     | K      | 四层土壤温度                    |
    | `TP`        | `tp`               | m/h    | 总降水（小时累计量）            |
  ]
])

*完整 ERA5 模块*：`src/Forcings/ERA5.jl` 实现了对 Fortran
`module_forcings.f90` 6 个核心强迫读取子程序的完整翻译（全球 $1440 times 721$
裁切为 $250 times 290$ 区域子域），与示例文件的路径约定不同，适用于
处理原始全球 ERA5 数据。

=== 1.3.2 LAI 月度气候态

#figure(
  caption: [LAI 月度气候态文件],
  align(center)[
  #three-line-table()[
    | *文件*                | *变量* | *维度*       | *单位* | *含义*                  |
    | ---                   | ---    | ---          | ---    | ---                     |
    | `../LAI/lai_clim_01.nc` | `lai` | x, y, month  | m²/m²  | 月度 LAI 气候态（12 月） |
  ]
])

路径为 `dirname(static) + "/../LAI/lai_clim_01.nc"`。代码按当月切片读取。

== 1.4 派生量与缺失值

以下变量由代码从原始输入计算，*用户无需提供*：

#figure(
  caption: [派生量计算公式],
  align(center)[
  #three-line-table()[
    | *派生量* | *公式*                                                              | *单位* |
    | ---      | ---                                                                 | ---    |
    | `qair`   | Tetens 公式：$e = 6.112 dot exp(17.67 (T_d - 273.15)/(T_d - 29.65)) dot 100$ Pa；$q = 0.622 e/(p - 0.378 e)$ | kg/kg |
    | `netrad` | $"ssrd" (1 - alpha) - "strd"$，短草 $alpha = 0.23$                    | W/m²   |
    | `precip` | $"tp" times 1000$                                                   | mm/h   |
  ]
])

异常值处理：

- `topo < -1e5` ⇒ 海洋格点（`landmask = 0`，`topo = 0`）
- `F < 1e-6` ⇒ `fdepth = 100.0` m
- `fd = 0` 或未列编码 ⇒ 汇水点 / 无效，`flowdir` 返回 `(0, 0)`
- 海洋格点（`landmask = 0`）由 `rootdepth_main` 跳过

#pagebreak()

= 2 常见问题

*Q1：`ssrd` / `strd` 单位是 J/m² 累计还是 W/m² 平均？*

A：示例直接以 W/m² 平均值读取；如使用 CDS 原始下载的累计 J/m²，
需 `÷3600` 转 W/m²。`src/Forcings/ERA5.jl` 内部已完成该转换。

*Q2：`tp` 单位是 m 还是 mm？*

A：示例中 `tp` 为 *m/h*，最终由 `precip = tp * 1000` 转为 mm/h 进入模型。

*Q3：海洋格点的水位与土壤含水量如何处理？*

A：`module_driver.f90` 中 `where (landmask .eq. 0) wtd = 0.` 将其置零；
`rootdepth_main` 通过 `landmask == 0 && continue` 跳过。建议在预处理
阶段将海洋格点的 `STXT` 设为 1（砂）。

*Q4：流向文件中部分格点的 `fd = 0` 是否正常？*

A：是。`fd = 0` 表示汇水点（pit）或海洋，不参与路由传播，调用方需保证
不向 `(0, 0)` 写入。

*Q5：多日模拟时如何扩展 ERA5 数据？*

A：每跨过 24 h 自动推进一天，调用 `era5_paths_for(era5_root, cur_date)`
按新日期重建 8 个强迫文件路径。用户需为模拟时段内每一天准备对应 NetCDF。

*Q6：网格面积 `area` 从哪里来？*

A：Fortran 端由 `READLATLON` 根据经纬度计算；Julia 端 P2 阶段实现前，
可在 `regional_example.jl` 中以 `area = 1e6` m²（1 km²）作为占位。


#pagebreak()
= 3 示例

== 3.1 Mock 模式

```bash
# 默认：在 mktempdir() 自动生成合成数据，运行 3 小时模拟
julia --project example/regional_example.jl --mock 1 --duration 3

# 指定网格大小和土壤层数
julia --project example/regional_example.jl --mock 1 --nx 10 --ny 10 --nzg 40 --duration 24
```

== 3.2 真实数据模式

```bash
julia --project example/regional_example.jl \
      --mock 0 --date 20200101 --duration 24 \
      --static /data/ASAP/static.nc \
      --wtd    /data/ASAP/wtd.nc \
      --era5   /data/ERA5 \
      --out    /tmp/asap_regional_20200101.nc
```

命令行参数：`--date` 起始日期、`--duration` 模拟小时数、`--nx/--ny/--nzg`
网格与土壤层（mock 模式）、`--mock 1|0` 模式切换、`--static/--wtd/--era5`
真实数据路径、`--out` 输出 NetCDF。

= 4 参考

- 静态场读取：`src/io/NetCDF.jl::read_initial`（L63–L100）
- 初始水位读取：`src/io/NetCDF.jl::read_wtdnc`（L113–L126）
- 流向工具：`src/helper.jl::flowdir`（L139–L162）
- 河道路由：`src/modules/rivers_kw_flood!`
- ERA5 完整模块：`src/Forcings/ERA5.jl::read_hourly_forcings`（L132–L250）
- 区域示例：`example/regional_example.jl::main`
- Fortran 参考：`fortran/module_io.f90::READINITIAL`（L10）、
  `READWTDNC`（L309）、`READFLOWDIRECTION`（L591）、
  `READRIVERPARAMETERS`（L693）；`fortran/module_wtable.f90::FLOWDIR`
- 项目政策：`CLAUDE.md §6.1`（NetCDF 写出 / MPI / ¹⁸O 同位素串联暂缓）
