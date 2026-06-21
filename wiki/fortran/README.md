# ASAP 模型 Fortran 代码总览

ASAP（Advanced Soil-Atmosphere Plant）模型是一个用于全球/大尺度水文过程的 Fortran 模型，主要模拟土壤水、地下水、河流、植被蒸散及稳定同位素（O18）示踪过程。模型在 MPI 并行框架下运行，针对各大洲分别模拟。

## 模块组成

| 文件 | 行数 | 功能 |
|------|------|------|
| `module_driver.f90` | 1118 | 主程序 `program driver`，负责初始化、时步循环、I/O 调度 |
| `main.f90` | 58 | （备用入口）调用 `LATERAL` / `ROOTDEPTH` / `GW2RIVER` / `FLOODING` 等子程序 |
| `module_rootdepth.f90` | 1914 | 核心模块：垂向水通量、土壤含水量、冠层截留、根系吸水、PET 计算 |
| `module_wtable.f90` | 1417 | 地下水侧向流、河流/洪水演进 |
| `module_forcings.f90` | 1615 | 读取 ERA5 / ERA5-Land / MODIS LAI 等强迫数据并插值到模式网格 |
| `module_io.f90` | 3599 | NetCDF I/O：静态场、初始场、强迫场、模式输出与历史场 |
| `module_initial.f90` | 1059 | 初始化（`INITIALIZE` / `EQSOILMOISTUREtheor` / `ZBRENT`） |
| `module_parallel.f90` | 658 | MPI 并行域分解、halo 通信 |
| `module_nrtype.f90` | 36 | 数值类型与常量定义（`SP`/`DP`/`I4B` 等） |
| `interp_lib.f90` | 690 | RAMS 大气模式移植的水平/垂直插值库（`gdtost2` / `gdtost3` 等） |
| `soilfluxes.f90` | 609 | 1D Richards 方程求解器（`SOILFLUXES`） |

## 典型时序（每 1 小时）

1. `READFORCINGS` 读 ERA5 风、温、压、湿度、辐射、降水
2. `READFORCINGSSNOW` 读 ERA5-Land 融雪
3. `READLAICHINA` 读 LAI
4. `LATERAL` 计算地下水侧向流
5. `GW2RIVER` 地下水补给河道
6. `ROOTDEPTH` 计算蒸散、土壤水、截留、补给
7. `FLOODING` / `RIVERS_KW_FLOOD` 河流与洪水演算
8. 累计日 / 月输出，写 NetCDF

## 强迫与参数数据

- 强迫：ERA5（小时）、ERA5-Land（融雪）、LAI（8 天）
- 静态：FAO 土壤质地、MERIT DEM、MODIS 植被与高度
- 13 种 USDA 土壤类型参数定义于 `module_rootdepth.f90`
