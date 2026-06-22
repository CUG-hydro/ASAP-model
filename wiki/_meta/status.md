# 摄取状态表

> 截至 2026-06-22 全部已建 wiki 页面的状态登记。
> 状态枚举：**已摄取** / **已摄取（含悬空 export 备注）** / **待复核** / **未建（待补）**。

## 1. 顶层结构（阶段 A 新建）

| 路径 | 行数 | 状态 | 备注 |
|---|---|---|---|
| `wiki/README.md` | 67 | 已摄取 | Wiki 总览、目录树、导航 |
| `wiki/index.md` | 80 | 已摄取 | 按子系统类别索引 |
| `wiki/log.md` | 28 | 已摄取 | 时间线（仅 2026-06-21 ingest 一条） |
| `wiki/conventions.md` | 110 | 已摄取 | 命名、单位、文档模板 |
| `wiki/julia/ASAP-主入口.md` | 145 | 已摄取 | 含悬空 export `updatewtd_qlat` 备注 |
| `wiki/julia/RootDepth-主算法.md` | 175 | 已摄取 | 含同位素追踪代码被注释的 §5 详述 |
| `wiki/_meta/status.md` | — | 已摄取 | 本文件 |
| `wiki/_meta/cross-refs.md` | — | 已摄取 | ≥ 10 条跨页面引用 |

## 2. Julia 模块页面

| 路径 | 对应源文件 | 状态 | 已知问题 |
|---|---|---|---|
| `julia/SoilParameters-土壤参数.md` | `src/SoilParameters.jl` | 已摄取 | `θ_wilt` 写死；`init_soil_param.fieldcp` 已实现 Fortran 公式（见 §5 #4） |
| `julia/SoilInitialization-土壤分层.md` | `src/SoilInitialization.jl` | 已摄取 | `slz` 长度 `nzg+1` 而 `dz` 仅 `nzg`；`DataFrames` 悬空 import 已清理（见 §5 #9） |
| `julia/SoilFluxes-土壤水运动.md` | `src/SoilFluxes.jl` | 已摄取 | 氧 18 注释段（按 §6.1 暂缓）；`transpdeep`/`qlat`/`qrf`/`flood` 形参未消费；`f_ice` 裸引用 |
| `julia/extraction-根系吸水.md` | `src/extraction.jl` | 已摄取（含悬空 export 备注） | `kroot = k - 1` 已澄清语义（见 §5 #5）；`dθ_deep` 始终 0 |
| `julia/Evapotranspiration-蒸散发.md` | `src/Evapotranspiration.jl` | 已摄取 | `potevap_shuttleworth_wallace` 别名已新增（见 §5 #8）；未 export 任何符号；`BIOPARMS` 索引 18-21/3-6 重复 |
| `julia/Interception-截留.md` | `src/Interception.jl` | 已摄取 | `minpprate` 已澄清为模块常量（见 §5 #7）；`interception` 已 export（L1） |
| `julia/updatewtd_shallow-浅层水位.md` | `src/updatewtd_shallow.jl` | 已摄取 | `freedrain` 形参未消费；`flag` 已显式初始化 + 注释加固（见 §5 #6） |
| `julia/lateral_flow-侧向地下水流.md` | `src/modules/lateral_flow.jl` | 已摄取 | 陆面掩膜 halo 未保护；`Δt` 未参与计算 |
| `julia/gw2river-地下水河流交换.md` | `src/modules/gw2river.jl` | 已摄取 | `apply_specific_diversions` 硬编码网格坐标（Taquari 4498/4535）；`length` 形参命名冲突 ⚠️（见 §5 #10） |
| `julia/rivers_kw_flood-运动波路由.md` | `src/modules/rivers_kw_flood.jl` | 已摄取 | `length` 已重命名为 `river_length`（见 §5 #10）；曼宁 `n=0.03` 硬编码 |
| `julia/rivers_dw_flood-扩散波路由.md` | `src/modules/rivers_dw_flood.jl` | 已摄取 | `width==0` 兜底 `√area` 经验式；与 KW 存在功能重叠；`length` 已重命名为 `river_length`（见 §5 #10） |
| `julia/flooding-洪泛漫流.md` | `src/modules/flooding.jl` | 已摄取 | `ntsplit=1` 固定；`delsfcwat<0` 仅源端约束；`bfd` 悬空参数 |
| `julia/IsotopeTracing-同位素追踪.md` | `src/modules/Tracing/IsotopeTracing.jl` | 已摄取 | 4 方向不含对角线；`*sum` 数组持续累加；`RootDepth.jl` 调用已禁用（见 §6.1） |
| `julia/modules-水位模块聚合.md` | `src/modules/Modules.jl` | 已摄取 | `wtable!`/`updatewtd!` export 已删除（见 §5 #2）；`helper.jl` 注释掉 |
| `julia/io-NetCDF.md` | `src/io/NetCDF.jl` | 已摄取（仅 P0 子集） | 仅实现 `read_initial` / `read_wtdnc`；P1/P2 子程序（`READLATLON`/`READVEG`/`READHISTORYNC`/全部写出）未实现 |
| `julia/Forcings-ERA5.md` | `src/Forcings/ERA5.jl` | 已摄取（含 Grid 重构备注） | 5 份插值收敛到 `bilinear`；`read_snow` 是 alias；6 形参 LAI 收为 1 Tuple |

## 3. Fortran 原版页面

| 路径 | 对应源文件 | 状态 | 已知问题 |
|---|---|---|---|
| `fortran/README.md` | `fortran/` 目录总览 | 已摄取 | — |
| `fortran/module_driver.md` | `fortran/module_driver.f90` (43 KB) | 已摄取 | — |
| `fortran/module_forcings.md` | `fortran/module_forcings.f90` (53.7 KB) | 已摄取 | — |
| `fortran/module_initial.md` | `fortran/module_initial.f90` (34.5 KB) | 已摄取 | — |
| `fortran/module_io.md` | `fortran/module_io.f90` (111.4 KB) | 已摄取 | — |
| `fortran/module_parallel.md` | `fortran/module_parallel.f90` (23.5 KB) | 已摄取 | — |
| `fortran/module_rootdepth.md` | `fortran/module_rootdepth.f90` (65.5 KB) | 已摄取 | — |
| `fortran/module_wtable.md` | `fortran/module_wtable.f90` (51.8 KB) | 已摄取 | — |
| `fortran/main.f90.md` | `fortran/main.f90` (2.5 KB, 58 行) | 已摄取 | 备用单步调度脚本：LATERAL→GW2RIVER→ROOTDEPTH→FLOODING→RIVERS_KW_FLOOD；无 module 包裹；与 `module_driver.f90` 主循环等价 |
| `fortran/soilfluxes.f90.md` | `fortran/soilfluxes.f90` (21 KB, 609 行) | 已摄取 | 1D Richards 求解器 `SOILFLUXES`（Crank-Nicolson 三对角 + ¹⁸O 同位素段）；对应 `src/SoilFluxes.jl::soilfluxes`；Julia 端同位素段已注释 |
| `fortran/interp_lib.f90.md` | `fortran/interp_lib.f90` (20 KB, 690 行) | 已摄取 | RAMS v4.3.0.2 插值库：13 个子程序（TRNCL1/2、INTRP、INTRRAP、BINOM、GDTOST/2/3、WEIGHTS、HTINT/2/HTINTCP、AWTCMP）；Julia 端无对应 |
| `fortran/module_nrtype.f90.md` | `fortran/module_nrtype.f90` (1.9 KB, 36 行) | 已摄取 | 数值类型与常量（I4B/SP/DP/EULER 等）+ 稀疏矩阵派生类型；Julia 端无对应 |

## 4. 跨语言映射（阶段 B）

| 路径 | 状态 | 规模 |
|---|---|---|
| `mapping/julia-fortran-对照.md` | 已摄取 | 17 Julia ↔ 11 Fortran 行级映射 |
| `mapping/algorithm-索引.md` | 已摄取 | 18 类物理过程反向索引 |

## 5. 全局问题汇总

> **状态枚举**：✅ 已修复 / ⏸️ 按 §6.1 暂缓 / 🟡 部分修复 / ⚠️ 待后续 PR 解决。
> 与 `CLAUDE.md §7` 一一对应；2026-06-22 lint 后更新。

| 优先级 | 问题 | 涉及页面 | 状态 | 修复位置 / 备注 |
|---|---|---|---|---|
| 高 | `RootDepth.jl` 同位素追踪注释掉 | `RootDepth-主算法.md` §5、`SoilFluxes-土壤水运动.md` §7 | ⏸️ 按 §6.1 暂缓 | `src/SoilFluxes.jl:347-L351`；恢复需用户明确重新授权 |
| 高 | `Modules.jl` 悬空 export `wtable!`/`updatewtd!` | `modules-水位模块聚合.md` §7 | ✅ 已修复 | `src/modules/Modules.jl:10-L16` 仅注释保留解释，`export` 已删 |
| 高 | `ASAP.jl` 悬空 export `updatewtd_qlat` | `ASAP-主入口.md` §7 | ✅ 已修复 | `src/ASAP.jl:17` `# include(...)` + export 行整体移除 |
| 中 | `SoilParameters.init_soil_param.fieldcp` 占位 0 | `SoilParameters-土壤参数.md` §7 | ✅ 已修复 | `src/SoilParameters.jl:73-L78` 实现 Fortran 公式 `θ_fc = θ_sat·(ψ_sat/-3.366)^(1/b)` |
| 中 | `extraction.kroot = k - 1` 待与 Fortran 复核 | `extraction-根系吸水.md` §7 | ✅ 已修复 | `src/extraction.jl:65-L70` 注释澄清「最后活跃层索引」语义 |
| 中 | `updatewtd_shallow.flag` 未初始化 | `updatewtd_shallow-浅层水位.md` §7 | ✅ 已修复 | `src/updatewtd_shallow.jl:35,39-L41` 显式 `flag = 0` + 注释 |
| 中 | `Interception.minpprate` 函数体未使用 | `Interception-截留.md` §7 | ✅ 已修复 | `src/Interception.jl:3` `const minpprate = 0.01`；L33/L44 实际消费 |
| 低 | `Evapotranspiration` 拼写错误 `shutteworth` | `Evapotranspiration-蒸散发.md` §7 | ✅ 已修复 | `src/Evapotranspiration.jl:7-L8` 双 t 旧名 + 拼写正确别名 `potevap_shuttleworth_wallace` 共存 |
| 低 | `SoilInitialization.DataFrames` 悬空 import | `SoilInitialization-土壤分层.md` §7 | ✅ 已修复 | `src/SoilInitialization.jl` 与 `Project.toml` 均已移除 `DataFrames` |
| 低 | `rivers_*` 中 `length` 与 `Base.length` 同名 | `rivers_kw_flood`/`rivers_dw_flood`/`gw2river` | 🟡 部分修复 | `rivers_kw_flood.jl:18`/`rivers_dw_flood.jl:13` 已重命名为 `river_length`；`gw2river.jl:23` 形参 `length::M` 仍存在，⚠️ 待后续 PR |

## 6. 子系统状态（examples）

| 路径 | 状态 | 备注 |
|---|---|---|
| `example/regional_example.jl` | 已摄取（mock+真实两模式+多日） | 区域端到端驱动；mock 模式自给自足，真实模式按滚动日期重建 ERA5 路径 |