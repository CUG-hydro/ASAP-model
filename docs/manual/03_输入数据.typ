// docs/manual/03_输入数据.typ
// ==========================
// 章：输入数据准备
// 状态：骨架（待摄取；与 wiki/julia/Forcings-ERA5.md §7 协调）

#import "../preamble.typ": *

= 输入数据准备

== 1. 目标

为 ASAP-model 准备区域应用所需的 ERA5 / ERA5-Land 强迫数据 + NetCDF 初始场。

== 2. 前置条件

- 已注册 CDS API（https://cds.climate.copernicus.eu/）
- 已安装 Python `cdsapi`（用于下载 ERA5 数据）
- 已完成 `manual/01_安装.typ`

== 3. 操作步骤

参见 README「区域应用：数据准备清单」7 小节（按变量列表准备 7 类 ERA5 字段）。

+ `2m_temperature`（t2m）
+ `2m_dewpoint_temperature`（d2m）
+ `surface_pressure`（sp）
+ `10m_u/v_wind_component`（u10/v10）
+ `surface_solar_radiation_downwards`（ssrd）
+ `total_precipitation`（tp）
+ `snow_depth`（sde）

初始场：地形（topo）、土壤类型（soiltxt）、水位（wtd）、植被类型（veg）、植被高度（hveg）、叶面积指数（LAI climatology）、土壤温度 4 层（stl1..4）。

== 4. 示例

参见 `example/regional_example.jl::era5_paths_for` 的日期感知路径构建。

== 5. 常见问题

- *Q：ssrd 单位是 J/m² 累计还是 W/m² 平均？*
  A：CDS 原始输出为累计 J/m²；`src/Forcings/ERA5.jl` 内部除以 3600 转为 W/m²。

== 6. 参考

- 数据读取：`src/Forcings/ERA5.jl`
- 初始场读取：`src/io/NetCDF.jl::read_initial`
- 水位读取：`src/io/NetCDF.jl::read_wtdnc`
- README 数据准备清单：#link("../README.md#区域应用数据准备清单")[README §区域应用]
- Wiki 数据细节：#link("../wiki/julia/Forcings-ERA5.md")[Forcings-ERA5.md]