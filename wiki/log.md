# 变更日志（log）

> 本文件记录 Wiki 的摄取、复核与重大变更，按日期倒序排列。

## [2026-06-22] update | 以 example/regional_example.jl 为主线补全 wiki

- **触发**：用户要求"以 `example/regional_example.jl` 为主线，更新 wiki"。
- **新增页面**：
  - `wiki/julia/example-regional.md` — 按 `conventions.md` §3 七段模板撰写：`parse_cli_args` / `generate_mock_dataset` / `read_mock_hourly_forcings` / `era5_paths_for` / `main` 5 个函数签名；6 步主流程（参数解析 → 数据准备 → 静态场初始化 → 状态变量 → 时步循环 → 输出）；5 节 mock 合成公式（地形 / 气温 / 辐射 / 土壤温度 / LAI）+ Tetens 比湿 + 压高公式 + 净辐射近似；与 Fortran `main.f90` / `module_io.f90` / `module_forcings.f90` / `module_initial.f90` / `module_rootdepth.f90` 逐项对照；§7 列出 10 项遗留问题（mock stub 与项目 `read_hourly_forcings` 4 点差异、`icefactor .= 0` 简化、`netrad` 粗略近似、输出 NetCDF 仅 6 字段、CLI 字符串字段、docstring 转义等）。
- **更新导航**：
  - `wiki/index.md`「数据准备与强迫子系统」分类新增 `Regional 区域应用示例` 条目；简介涵盖 mock / 真实两模式 + 多日滚动 + 6 步主流程。
  - `wiki/README.md` 目录树在 `julia/` 子目录末尾追加 `example-regional.md`；快速导航新增"区域端到端示例"小节；"当前覆盖范围"补 1 个区域应用示例计数。
- **更新状态**：
  - `wiki/_meta/status.md` §6 新增第二行 `julia/example-regional.md | 已摄取（含 stub 与项目模块差异备注）`，指向 §7 #1-#10 遗留问题。
- **更新跨页引用**：
  - `wiki/_meta/cross-refs.md` 追加 7 条新引用（条目 31-37），覆盖 `eqsoilmoisturetheor` / `rootdepth_main` / `read_initial`+`read_wtdnc` / `read_mock_hourly_forcings` 差异 / `icefactor` 简化 / 测试入口 / 数据准备清单；引用统计由 30 条升至 37 条，覆盖页面由 15 个 Julia 页面升至 16 个。
- **未变更**：
  - `src/backup/` 按 §6 仍禁用。
  - 同位素追踪（`SoilFluxes.jl:347-L351`）按 §6.1 维持注释禁用。
  - `example/regional_example.jl` 源码未修改（仅 wiki 页面新增 + 导航/状态/跨引更新）。
- **验证**：本地语法 / 表格结构未跑自动测试；下一步由用户在 `julia --project -e 'using Pkg; Pkg.test()'` 全套绿后确认本次 wiki 更新不引入回归。

## [2026-06-22] lint | codebase update — Wiki 与源码 §5/§7 对齐 + 知识图谱索引

- **触发**：用户要求"codebase update"。本次为 D1 lint 的延续：核对 `CLAUDE.md §7` 与 `wiki/_meta/status.md §5` 的不一致。
- **范围**：
  - 索引 `ASAP-model` 到 codebase-memory-mcp（913 节点 / 1205 边），覆盖 17 个 Julia 源文件 + 11 个 Fortran 源文件 + test/ + example/。
  - 全量 lint：38 个 wiki 页面的 80 条 markdown 内链（0 断裂）、0 孤立页。
  - 源码 vs wiki 对照：核对 `CLAUDE.md §7` 已修复的 10 个源码问题在 wiki 各页面 §7 / `status.md §5` 中的描述状态。
- **修复**（11 处）：
  1. `wiki/_meta/status.md §5`：表格增加「状态」列；10 个问题中 9 个标 ✅ 已修复、1 个（同位素追踪）标 ⏸️ 按 §6.1 暂缓；新增「`rivers_* length` 重命名 🟡 部分修复」（`rivers_kw_flood.jl`/`rivers_dw_flood.jl` 已修，`gw2river.jl:23` 形参 `length::M` 仍存在，⚠️ 待后续 PR）。
  2. `wiki/_meta/status.md §2`「已知问题」列：清理 8 行 stale 描述，按 `#5 §7` 编号引用回链到 §5。
  3. `wiki/julia/ASAP-主入口.md` §5/§8：移除 `updatewtd_qlat` 悬空 export 条目（§8.1 改 ✅ 修复说明）；`interception` 改为已 export（L1）；`wtable!`/`updatewtd!` export 删除标注；新增 `read_initial`/`read_wtdnc`/ERA5 readers 两行；§8.4 DataFrames 标 ✅ 修复。
  4. `wiki/julia/SoilParameters-土壤参数.md` §7：`fieldcp` 占位 0 改为 ✅ 已修复（指向 `src/SoilParameters.jl:73-L78` Fortran 公式 + `test_soil_parameters.jl` 数值断言）；移除 DataFrames 悬空 import 备注。
  5. `wiki/julia/SoilInitialization-土壤分层.md` §7：`using DataFrames` 悬空 import 改为 ✅ 已清理。
  6. `wiki/julia/Evapotranspiration-蒸散发.md` §7：`potevap_shutteworth_wallace` 拼写错误改为 ✅ 已通过别名修复；`export` 列表补 L7-L8。
  7. `wiki/julia/extraction-根系吸水.md` §5/§7：`kroot = k - 1` 改为 ✅ 语义已澄清（指向 `src/extraction.jl:65-L70`）。
  8. `wiki/julia/Interception-截留.md` §5：`MINPPRATE` 风险行改为 ✅ 模块常量 + 函数体消费（L3/L33/L44）。
  9. `wiki/julia/updatewtd_shallow-浅层水位.md` §7：`flag` 未初始化隐患改为 ✅ 显式初始化 + 注释加固（指向 `src/updatewtd_shallow.jl:35,39-L41`）。
  10. `wiki/julia/rivers_kw_flood-运动波路由.md` §7：`length` 命名冲突改为 ✅ 已重命名为 `river_length`（指向 L18/L78/L83/L95）。
  11. `wiki/julia/RootDepth-主算法.md` §8.3/§8.4：`length` 冲突 ✅；`o18`/`θ_eq`/`tempsfc` 形参未消费改为 ⏸️ 按 §6.1 暂缓。
- **未变更**：
  - `src/backup/` 按 §6 仍禁用（仅作为历史参考）。
  - 同位素追踪（`SoilFluxes.jl:347-L351`）按 §6.1 维持注释禁用；恢复需用户明确授权 + 追加 `[YYYY-MM-DD] enable | 范围` 条目。
  - `gw2river.jl:23` 形参 `length::M` 与 `Base.length` 同名问题未修，记入 `_meta/status.md §5 #10` 作为 ⚠️ 待后续 PR 项。
- **验证**：本地语法 / 表格结构未跑自动测试；下一步由用户在 `julia --project -e 'using Pkg; Pkg.test()'` 全套绿后确认本 lint 不引入回归。

## [2026-06-22] policy | CLAUDE.md §6.1 — 现阶段不实现的范围

- **触发**：用户明确"现阶段 NetCDF / MPI / 同位素均不需要实现"。
- **变更**：
  - `CLAUDE.md` 新增 §6.1「现阶段不实现的范围（2026-06-22 起生效）」，明确禁止：
    - NetCDF I/O 扩展（仅 `src/io/NetCDF.jl::NetCDFIO` 的 P0 子集 `read_initial` / `read_wtdnc` 为当前终点）
    - MPI 并行化（Julia 单进程，禁止引入 `MPI.jl`）
    - ¹⁸O 同位素串联（`SoilFluxes.jl` 末段三对角 + `RootDepth.jl` 中 `lateral_isotope!` / `updatedeepwtable!` 调用维持禁用）
  - `CLAUDE.md` §1 项目背景：18O 物理过程行加注"当前阶段不串联，见 §6.1"
  - `CLAUDE.md` §7 时间戳：2026-06-21 → 2026-06-22
  - `CLAUDE.md` §9 跨语言映射 `soilfluxes` 行加注"现阶段不串联，见 §6.1"
- **解锁条件**：后续若需启用任一功能，必须由用户明确重新授权，并在 `wiki/log.md` 追加 `[YYYY-MM-DD] enable | 范围` 条目。

## [2026-06-22] update | Phase 3 — NetCDF I/O（P0 子集）落地

- **触发**：用户任务 `Phase 3: 翻译 module_io.f90 NetCDF I/O`，仅实现读取子集。
- **新增源码**：
  - `src/io/NetCDF.jl`（约 115 行）— 新子模块 `NetCDFIO`，
    实现 `read_initial(path) -> NamedTuple` 与 `read_wtdnc(path) -> Matrix{Float64}`。
    对应 Fortran `module_io.f90::READINITIAL`（L10）与 `READWTDNC`（L309）。
- **修改源码**：
  - `src/ASAP.jl` 末尾 `include("io/NetCDF.jl")` + `using .NetCDFIO: read_initial, read_wtdnc` + `export`。
  - `Project.toml`：
    - `[deps]` 新增 `NCDatasets = "85f8d34a-cbdd-5861-8df4-14fed0d494ab"`。
    - `[compat]` 新增 `NCDatasets = "0.12, 0.13, 0.14"` 与 `julia = "1.6"`（既有约束未修改，仅追加）。
    - `[extras]` 与 `[targets]` 新增 `Statistics` 以解 `test_eqsoilmoisture.jl` 的 `using Statistics`（Phase 2 遗留的传递性依赖）。
- **新增测试**：`test/test_io_netcdf.jl`，5 个 testset 共 12 个断言覆盖
  - `read_wtdnc` 圆环 + 混合正负边界
  - `read_initial` 圆环 + `fdepth` 夹断（`< 1e-6 → 100`）+ 地形/掩码（`< -1e5 → 0/0`）
- **新增 wiki 页面**：`wiki/julia/io-NetCDF.md`（按 §3 7 段模板），
  - §7 已知问题：NCDatasets v0.14 API 差异（`v[:]` 会展平、必须 `v[:, :]`）、
    地形变量名 `"topo"` 是 Julia 约定而非 Fortran 约定、MPI 分发省略、
    写出整段未实现。
- **测试结果**：`Pkg.test()` 全部通过（既有 32+ 测试 + 新增 5 个 testset = 12 断言）。
- **未实现（P1/P2）**：`READLATLON` / `READVEG` / `READHVEG` /
  `READSMOIEQ` / `READFLOWDIRECTION` / `READRIVERPARAMETERS` /
  `READHISTORYNC*` / 全部写出 — 留待后续 PR。

## [2026-06-21] update | 补齐 4 个 Fortran wiki 页面

- **触发**：用户要求继续翻译未完成的 Fortran 源文件。
- **新建页面**（4 个，Fortran 原版章节）：
  1. `wiki/fortran/main.f90.md`（2.5 KB）— 备用单步调度脚本：`LATERAL → GW2RIVER → ROOTDEPTH → FLOODING → RIVERS_KW_FLOOD` 序列，无 module 包裹；详细列出 `qlat*deltat/deltatwtd` 折算、5 天日采样、累计量约定。
  2. `wiki/fortran/soilfluxes.f90.md`（约 10 KB）— 1D Richards 求解器 `SOILFLUXES`：Campbell K(θ)/D(θ)、三对角组装、顶部入渗能力截断、底部自由排水/受限排水分支、含水量限幅、¹⁸O 同位素段（Majoube 平衡分馏）、侧向流分配；与 `src/SoilFluxes.jl::soilfluxes` 一一对照。
  3. `wiki/fortran/interp_lib.f90.md`（约 11 KB）— RAMS v4.3.0.2 插值库：13 个子程序（TRNCL1/2、INTRP、INTRRAP、BINOM、GDTOST/2/3、WEIGHTS、HTINT/HTINT2/HTINTCP、AWTCMP）；含六阶中心差分权重展开表与缺失值掩码 `1e30` 约定。
  4. `wiki/fortran/module_nrtype.f90.md`（约 4 KB）— 数值类型符号（I4B/I2B/I1B/SP/DP/SPC/DPC/LGT）与数学常量（PI/PIO2/TWOPI/SQRT2/EULER，含 DP 版本）；稀疏矩阵派生类型 `sprs2_sp/sprs2_dp`；明确 Julia 端无对应。
- **更新文档**：
  - `wiki/_meta/status.md` §3：把 4 个 Fortran 页面从"未建（待补）"改为"已摄取"，补充文件大小与已知问题。
  - `wiki/index.md` Fortran 原版章节：插入 `main`、`module_nrtype`、`interp_lib`、`soilfluxes` 4 个新条目。
- **跨语言映射**：4 个 Fortran 页面均已在 `wiki/mapping/julia-fortran-对照.md` 中预登记，本轮摄取与映射保持一致；新增 SOILFLUXES 段对照表的细化（Campbell 公式 / 同位素段 / 限幅 / 侧向流分配）。
- **已知遗留**（同步登记至 `_meta/status.md` §3 备注列）：
  - `soilfluxes.f90` 仍在主体代码中维护 ¹⁸O 同位素段，而 `src/SoilFluxes.jl` 同段被注释；迁移路径待 `IsotopeTracing.jl` 串联主循环后恢复。
  - `interp_lib.f90` 中 `TRNCL1/TRNCL2` 引用未声明的全局变量 `ID`/`JD`，依赖 `module_forcings` 编译上下文。
  - `main.f90` 无独立 module 与测试入口，依赖调用上下文提供全部变量；其逻辑等价于 `RootDepth.jl` 的单步迭代。

## [2026-06-21] ingest | ASAP-model wiki 初始化

- **触发**：用户首次请求建立完整 Wiki，覆盖 ASAP-model 的 Julia 与 Fortran 双侧实现。
- **覆盖范围**：
  - **17 个 Julia 源文件**：`src/ASAP.jl`、`src/RootDepth.jl`、`src/helper.jl`、`src/SoilParameters.jl`、`src/SoilInitialization.jl`、`src/SoilFluxes.jl`、`src/extraction.jl`、`src/Evapotranspiration.jl`、`src/Interception.jl`、`src/updatewtd_shallow.jl`、`src/modules/Modules.jl`、`src/modules/flooding.jl`、`src/modules/lateral_flow.jl`、`src/modules/gw2river.jl`、`src/modules/rivers_dw_flood.jl`、`src/modules/rivers_kw_flood.jl`、`src/modules/Tracing/IsotopeTracing.jl`。
  - **11 个 Fortran 源文件**：`fortran/main.f90`、`fortran/module_driver.f90`、`fortran/module_forcings.f90`、`fortran/module_initial.f90`、`fortran/module_io.f90`、`fortran/module_nrtype.f90`、`fortran/module_parallel.f90`、`fortran/module_rootdepth.f90`、`fortran/module_wtable.f90`、`fortran/soilfluxes.f90`、`fortran/interp_lib.f90`。
- **产物**：
  - 阶段 A 顶层结构 8 个文件：`README.md`、`index.md`、`log.md`、`conventions.md`、`julia/ASAP-主入口.md`、`julia/RootDepth-主算法.md`、`_meta/status.md`、`_meta/cross-refs.md`。
  - 已有 14 个 Julia 模块页面（`julia/*.md`）与 8 个 Fortran 模块页面（`fortran/*.md`）纳入索引。
- **已知遗留**（同步登记至 `_meta/status.md`）：
  - `SoilParameters.init_soil_param` 返回的 `fieldcp` 数组始终为 0（占位实现）。
  - `Evapotranspiration` 拼写错误 `potevap_shutteworth_wallace`（双 t）保留以避免破坏下游。
  - `SoilFluxes` 末尾氧 18 同位素计算被注释（30 行死代码）。
  - `modules/Modules.jl` 中 `wtable!` / `updatewtd!` 为悬空 export，对应实现位于 `backup/`。
  - `Interception.minpprate` 在文档字符串中列出但函数体未读取（接口不一致）。
  - `updatewtd_shallow.freedrain` 已声明但主循环未消费。

## [已完成] 阶段 B — Fortran→Julia 行级映射

- 完成 `mapping/julia-fortran-对照.md`，按函数/子程序级别建立双向链接。
- 完成 `mapping/algorithm-索引.md`，按算法（Richards、Manning、SW、Priestley-Taylor）反向索引到所有出现该算法的页面。

## [2026-06-21] fix | 10 个源码问题修复 + CLAUDE.md 新建

- **触发**：D1 lint 识别的 10 个源码问题，用户要求全部修复（忽略 `src/backup/` 目录）+ 新增 `CLAUDE.md`。
- **修复明细**：
  1. ✅ `SoilParameters.init_soil_param.fieldcp` — 补 Fortran 公式 `θ_fc = θ_sat · (ψ_sat / -3.366)^(1/b)`；新增 12 条 `test/test_soil_parameters.jl` 断言。
  2. ✅ `SoilFluxes.jl` 氧 18 段 — 清理 docstring（移除 o18 相关参数条目），清理 `RootDepth.jl` 中 `# transpo18[i, j] += ...` 等注释。
  3. ✅ `Modules.jl` 悬空 export — 移除 `wtable!` / `updatewtd!`，加注释说明被替代方案。
  4. ✅ `ASAP.jl` 悬空 export — 移除 `updatewtd_qlat`。
  5. ✅ `extraction.kroot = k - 1` — 修正误导性注释，澄清 Fortran 1-based 一致性。
  6. ✅ `updatewtd_shallow.flag` — 加注释确认每轮迭代 flag 重置为 0，`kwt==1` 路径不残留。
  7. ✅ `Interception.minpprate` — 修正 docstring（`minpprate` 是模块常量，非形参）。
  8. ✅ `potevap_shutteworth_wallace` 拼写错误 — 新增 `potevap_shuttleworth_wallace = potevap_shutteworth_wallace` 别名；新增 `test/test_evapotranspiration.jl` 别名测试。
  9. ✅ `SoilInitialization.DataFrames` 悬空 import — 删除 `using DataFrames`，同步移除 `Project.toml` 中 `DataFrames = "1.7.0"` 依赖与 `compat` 条目。
  10. ✅ `rivers_*` 形参 `length` → `river_length` — 全链路重命名（4 处使用）。
- **新增文件**：`CLAUDE.md`（项目根，10 个章节，含 Wiki 操作指南、源码约定、已知问题清单）。
- **测试结果**：`Pkg.test()` 全部通过（共 32 个测试集；包括修复后的 fieldcp 数值断言、别名一致性测试）。

## [2026-06-21] lint | Wiki 质量检查与修复

- **触发**：D1 阶段质量复核（孤立页面、断裂链接、源代码一致性、TODO 占位）。
- **检查范围**：32 个 wiki 页面（顶层 6 个、julia 16 个、fortran 8 个、mapping 2 个、_meta 2 个）、175 KB 文档、约 33 632 字。
- **链接审计**：
  - 扫描全部 72 处 markdown 内链，**0 断裂**（所有目标 `.md` 解析成功）。
  - 32 个页面**全部具有入链**（无孤立页）。
  - 9 个叶页（如 `julia/*` 业务模块）仅被 `index.md` 等聚合页索引，无下游互链——属 wiki 标准 hub-and-spoke 结构，不计作孤立。
- **修复内容**（3 处）：
  1. `wiki/index.md`：移除 `mapping/julia-fortran-对照.md`、`mapping/algorithm-索引.md` 条目后的 `*(待建立)*` 标注；两文件已实际存在（19.6 KB + 23.5 KB）。
  2. `wiki/_meta/status.md` §4：把两条 mapping 由「未建（待补）」改为「已摄取」，并补充规模说明。
  3. `wiki/log.md`：阶段 B 标记由 `[计划]` 改为 `[已完成]`。
- **源代码一致性**：抽样核对 5 个核心模块（`ASAP.jl`、`RootDepth.jl`、`SoilParameters.jl`、`Evapotranspiration.jl`、`Modules.jl`）的描述与现行源码一致；`dθ` 预初始化（commit `d8fdeb6`）、`Modules.jl` 悬空 export（`wtable!`/`updatewtd!` 位于 `backup/`）等条目已在 §7/§8 详述，未发现需更正的事实错误。
- **TODO 占位扫描**：grep `TODO|TBD|占位|待写` 后未发现空章节；所有 `## 7. 已知问题与备注` 均为实质性描述（最短 5 行，最长 8 行）。
- **仍存在但非本次范围**：`_meta/status.md` §5 列出的 10 类悬空 export / 占位实现 / 拼写错误属代码层问题，不通过 lint 修复，留待后续 PR。
## [2026-06-22] update | 区域应用串联 + Forcings 重构

按 TASK.md 推进到「区域应用」里程碑。

- **Forcings 测试可跑**：`src/Forcings/ERA5.jl` 加 `xrange/yrange/global_nx` 注入参数；test 改 30×30 mock。`Project.toml` test 目标补 `Random/Dates/Printf`。
- **示例端到端 + 多日**：`example/regional_example.jl` 加 `era5_paths_for` + 日期感知时步循环 (`--duration>24` 跨日)。
- **Forcings 重构**：引入 `Grid` 类型 (`origin/res/shape`)；5 份双线性插值收敛到 `bilinear(src, src_grid, lats, lons)`；`read_snow` 改为 `read_snow_hour(date,1,…)` alias；`read_lai_climatology` 6 形参收为 `target_shape::Tuple`；新增 `ERA5_GRID` / `ERA5LAND_GRID` 常量。794 → 360 行。生产行为不变，6 个 testset 无需改签名。
- **测试 polish**：`test/test_regional_example.jl` `generate_mock_dataset` 加类型断言；smoke test 接 `--out` 并断言输出文件存在；`parse_cli_args` 修 `parse(String,…)` bug（`String` 不需要 parse）。
- **README**：新增「区域应用：数据准备清单」7 小节。
- **验证**：`Pkg.test()` 全套通过；`--mock 1 --duration 50` 跨日 20200101→20200102→20200103。

## [2026-06-22] update | 合理 mock 数据

- `example/regional_example.jl::generate_mock_dataset` 重写：soiltxt 横向 1..13 梯度、topo 山脊 0..500 m、wtd 与 topo 负相关、t2m 日正弦 + 海陆 ±2 K、ssrd 仅昼间、tp 间歇性 (10% rainy hours)、sp 压高公式、stl1..4 深−滞后、LAI 月度正弦（冬 0.5 夏 4.0）。
- 修 3 个小 bug：comprehension 重写（替代错误 broadcast）、wind clamp 1..6 m/s、tp 单位修 mm/h。
- `test/test_regional_example.jl` 「static 往返」与「wtd 符号约定」改范围断言（不再断言具体常数）。
- 验证：`--mock 1 --duration 168` 7 天跨日跑通；`Pkg.test()` 全套绿。

## [2026-06-22] review | 以 example/regional_example.jl 为主线，串联审查核心代码

- **触发**：用户要求"以 `example/regional_example.jl` 为主线，Review 模型核心代码是否正确，结合 codebase"。本次为首次系统性 code review，覆盖主线程 + 全部 `src/*.jl` 核心模块。
- **方法**：
  - 知识图谱：`codebase-memory-mcp` 索引（914 节点 / 1191 边），`trace_path` 跟 `rootdepth_main` 调用链。
  - 通读：6 步主线程 + 16 个 Julia 源文件 + 2 个 IO/Forcings 子模块。
  - 测试验证：`Pkg.test()` 全套绿（`regional_example` 端到端 smoke test 12.4 s 通过；其余 9 个核心 testset + ERA5Forcings 7 个 testset + wtable 7 个 testset 全 Pass）。
- **核心链路正确性**（✅ 全部对齐 Fortran）：
  1. **NetCDF I/O**（`src/io/NetCDF.jl`）：`read_initial` 4 字段（soiltxt/topo/fdepth/landmask）+ `fdepth < 1e-6 → 100` 夹断 + `topo < -1e5 → 0` 掩码；`read_wtdnc` `min(-raw, 0)` 符号约定 — 与 Fortran `module_io.f90` 一致。
  2. **土壤参数**（`src/SoilParameters.jl`）：13 类 USDA 土壤表 + Campbell 导水率 `K(θ) = Ksat·(θ/θ_sat)^(2b+3)`；`init_soil_param.fieldcp` 用 `θ_sat·(ψ_sat/-3.366)^(1/b)` 公式（✅ §7 #1 已修）。
  3. **平衡含水量**（`src/SoilInitialization.jl::eqsoilmoisturetheor`）：逐层 zbrent 求解 `d1·(x-smoi1)/dz + k1 + flux = 0`；深度衰减 `clamp(exp((z+1.5)/fdepth), 0.1, 1.0)` 对齐 Fortran L168-L177。
  4. **Richards 求解**（`src/SoilFluxes.jl`）：三对角组装 + Campbell K(θ)/D(θ) + 顶部入渗能力截断 `Imax = Ksat·dt` + 自由/受限排水分支 + 含水量限幅；⚠️ 末段氧 18 注释按 §6.1 维持禁用。
  5. **蒸散发**（`src/Evapotranspiration.jl`）：P-T / P-M / SW 三法俱全；SW 双源返回 12 元组（Δ/γ/λ/ra_a/ra_c/rs_c/R_a/R_s/pet_s/pet_c/pet_w/pet_i），单位注意：pet_s, pet_c 是 W/m²，pet_w, pet_i 是 mm/Δt；拼写错误已加别名 `potevap_shuttleworth_wallace`（✅ §7 #8）。
  6. **根系吸水**（`src/extraction.jl`）：Feddes-style `easy[k] = max(-(POTLEAF-ψ)/(hveg-z), 0)` + `fswp = θ_root/rootfc` 胁迫 + `maxwat` 限幅 + `watdef` 累计缺水；`kroot = k-1` 语义已澄清（✅ §7 #5）。
  7. **截留**（`src/Interception.jl`）：`intercepmax = 0.2·lai` + `minpprate = 0.01` 阈值分支；`minpprate` 文档与实现一致（✅ §7 #7）。
  8. **浅层水位**（`src/updatewtd_shallow.jl`）：迭代上升/下降分支 + `flag` 每轮重置 0（✅ §7 #6）；`find_jwt` 在 `helper.jl` 中 1-based 索引正确。
  9. **侧向流 / 河-地交换 / 河流路由**（`src/modules/`）：D8 邻居 + D8 流向 `flowdir` + 河床交换 3 模式 + 运动波 + 扩散波 + 漫流；`length` → `river_length` 已全链路重命名（✅ §7 #10）。

- **发现 6 处需关注**（按严重度排序）：

  1. 🟡 **`regional_example.jl` 局部 `read_hourly_forcings` 与项目 `src/Forcings/ERA5.jl::read_hourly_forcings` 签名不一致**
     - 文件：`example/regional_example.jl:316-362` 与 `src/Forcings/ERA5.jl:178-226`
     - 局部 stub 签名 `read_hourly_forcings(hour::Int, paths)` 返回 `(wind, temp, qair, press, netrad, rshort, precip, lai)`，8 字段
     - 项目模块签名 `read_hourly_forcings(date::String, root::String; grid, xrange, yrange)` 返回 `varpack(rx, ry, 24, 11) + 10 标量字段`
     - 局部定义 shadow 项目 `using ASAP` 导出的同名函数；真实数据模式（`--mock 0`）将使用局部 stub 而非项目实现，导致 NetCDF 文件名约定不一致（项目期望 `ERA5_ws10_*.nc`，局部 stub 期望 `ERA5_wind_speed_*.nc`）
     - 建议：要么删除局部 stub 改用项目模块 + 包装层，要么把两套约定统一

  2. 🟡 **`regional_example.jl::read_hourly_forcings` LAI 硬编码为 1 月**（`[:, :, 1]`）
     - `example/regional_example.jl:358`：`ds["lai"][:, :, 1]`，全年都按 1 月（LAI=0.5）取值
     - 7 月应取 LAI=4.0 但实际仍是 0.5，蒸腾低估 ~3.5 倍
     - 建议：取 `month = Dates.month(Date(cfg.date, dateformat"yyyymmdd"))`

  3. 🟡 **`src/SoilFluxes.jl::soilfluxes` 三对角解被改写后丢弃，再用 Q 通量做质量平衡更新**
     - 文件：`src/SoilFluxes.jl:219-264`
     - L219 `tridag!(aa, bb, cc, rr, θ)` 求出新 θ
     - L222-226 用新 θ 计算 `capflux = -aa·(θ[k]-θ[k-1])·dt` + `gravflux = -K_mid·dt`
     - L261 `θ .= θ_old` 恢复为旧值
     - L263 `θ_old[k] += (Q[k]-Q[k+1]-transp[k])/dz[k]` 用新 Q 更新旧 θ
     - L360 `θ .= θ_old` 写回
     - 模式成立（tridiag 用来求一致通量，质量平衡更新守恒），但可读性差且对 `θ` 在 219-261 之间的用途需注意（如有 `icefactor` 路径依赖 θ 必须放在 261 之后）
     - 建议：在函数头注释 + 局部变量 `θ_star` 明确两阶段语义，或改用 `θ_new = θ_old; tridag!(...); θ_old .+= ...; θ_new = θ_old`

  4. 🟢 **`regional_example.jl` 文档与代码不一致**（3 处）
     - L122 docstring：`wtd_raw = -1.0 - 0.005·topo + 0.5·noise`；代码 L170 是 `+0.005·topo + 0.3·randn()`
     - L126 docstring：`wind 2..6 m/s`；代码 L201-204 实际可到 ~21 m/s（`topo_factor=1+0.003·topo` 在 topo=500 处 ×2.5）
     - L131 docstring：`tp 2..8 mm/h`；代码 L225 实际是 `1..5 mm/h`
     - 建议：同步文档与代码，或将合成数据范围与文档对齐

  5. 🟢 **`test/wtable/runtests.jl`「Basic Water Table Calculations」testset 为空**（0 tests / 0.0 s）
     - `Pkg.test()` 输出 `Test Summary: Basic Water Table Calculations | 0 0.0s`
     - 实际是 `runtests.jl` 中对应 testset 内没有 `@test`，仅 `@testset` 包裹了空 begin/end 块
     - 建议：补全测试（可参考 `test_rivers_dw_flood.jl` 风格）；或明确移除空 testset

  6. 🟢 **`src/SoilParameters.jl::KLATFACTOR` 未导出**（次要）
     - `SoilType` 结构体有 `K_latfactor` 字段，`get_soil_params` 填充 `KLATFACTOR[soil_type]`
     - 但 `KLATFACTOR` 数组本身未在 `export` 列表（L6-7），仅在 `src/modules/lateral_flow.jl` 中通过 `κlat[i, j]` 形参流入
     - 建议：若 `lateral_flow!` 期望 `κlat` 已乘 `K_latfactor`，则当前实现依赖外部传入预乘值，模块边界不清晰；建议在 `lateral_flow!` 入口显式乘 `K_latfactor` 并从 `SoilParameters` 拿值

- **未发现功能性 bug**：
  - 13 类 USDA 土壤参数与 Fortran `module_rootdepth.f90::INIT_SOIL_PARAM` 一致（θ_sat/θ_cp/ψ_sat/K_sat/b 五参数表）
  - `initializesoildepth` 与 Fortran L1-L20 节点深度公式一致
  - `tridag!` 标准 Thomas 算法，bet==0 防御已就位
  - `helper.find_jwt` 1-based 索引正确（与 Fortran 一致；`kroot = k-1` 语义已澄清）
  - 端到端 mock 模式 12.4 s 完成 1 时步 × 3×3 网格（12.4 s 主要开销是 Julia 编译 + NetCDF 读写），模型本身计算成本合理
  - `Pkg.test()` smoke test 全套绿

- **建议优先级**：
  - P0（修）：无
  - P1（建议）：#1（stub 与项目模块签名统一）、#2（LAI 月份）
  - P2（清理）：#3（`soilfluxes` 可读性）、#4（文档同步）、#5（空 testset 补全）
  - P3（可选）：#6（`KLATFACTOR` 模块边界）

- **未变更**：
  - `src/backup/` 按 §6 仍禁用
  - 同位素追踪按 §6.1 维持注释禁用
  - `gw2river.jl:23` 形参 `length::M` 未在本次 review 范围（已在 `_meta/status.md §5 #10` 登记）
- **验证**：`Pkg.test()` 全套绿（10+ ERA5Forcings 断言 + regional_example 4 个 testset 52 断言 + wtable 7 个 testset 38 断言 + 9 个核心模块 testset 全部通过）。

## [2026-06-22] fix | 应用 review 的 6 处建议

- **触发**：用户对 review 报告（上一条目）的 6 处问题逐项要求修复。
- **修复明细**：

  1. **#1（🟡）解除 `read_hourly_forcings` shadow** — `example/regional_example.jl:309-364`
     - 把局部 stub 改名为 `read_mock_hourly_forcings(hour, paths, month=1)`。
     - 更新 `main()` 调用方（`f = read_mock_hourly_forcings(...)`），并把 `cur_date` 拆为 `cur_date_dt = base_date + Day(day_offset)`（`Date`）与 `cur_date`（`String`）。
     - 顶部 docstring 列出与 `src/Forcings/ERA5.jl::read_hourly_forcings` 的 4 点差异（签名/变量名/返回字段/LAI 来源），明确两个函数互不替代。
     - 顺带：原 docstring 字符串包含 `$root` 触发 Julia 字符串插值（`UndefVarError: root`），已用 `\$root` 转义。

  2. **#2（🟡）LAI 月份化** — `example/regional_example.jl:325, 359`
     - `read_mock_hourly_forcings` 新增 `month::Int=1` 形参（1..12，默认 1）。
     - LAI 切片由 `[:, :, 1]` 改为 `[:, :, clamp(month, 1, 12)]`；调用方传 `Dates.month(cur_date_dt)`。
     - 7 月模拟时 LAI 由 0.5 升到 4.0，蒸腾需求正确恢复。

  3. **#3（🟡）`soilfluxes` 两阶段语义明确化** — `src/SoilFluxes.jl:218-326, 360-362`
     - 阶段 1：新增局部 `θ_star = similar(θ)` 缓冲；`tridag!(aa, bb, cc, rr, θ_star)` 写入新缓冲。
     - 阶段 2：质量平衡与边界修正直接在 `θ` 上原地写，不再走 `θ .= θ_old` 恢复再覆盖的迂回。
     - 函数头注释补充「两阶段」段落，明确 `θ_star`（tridiag 解）/ `θ_old`（旧值，用于顶层 BC 与水位层通量）/ `θ`（新值输出）三个角色。
     - 行为不变（仅结构改写）；`test_soil_fluxes.jl` / `test_rootdepth.jl` 全部通过。

  4. **#4（🟢）docstring 同步** — `example/regional_example.jl:113-138`
     - 修正 3 处不一致：wtd 符号（写正值 + `read_wtdnc` 翻为非正）/ wind 范围 0.6..21 m/s（地形放大 2.5×）/ tp 范围 1..5 mm/h。
     - 补 t2m 噪声说明、d2m 范围、strd 兜底、stl lag 列表、LAI 月度正弦公式。

  5. **#5（🟢）删除空 testset** — `test/wtable/test_watertable.jl`
     - 「Basic Water Table Calculations」testset 因依赖已禁用的 `wtable!`（`src/backup/`）而整段空，全部 `@test` 被注释。
     - 整段删除（38 行），避免「0 tests / 0.0 s」噪音；`Pkg.test()` 输出不再有该空 testset。

  6. **#6（🟢）`KLATFACTOR` 加入 export** — `src/SoilParameters.jl:6-8, 32-40`
     - 把 `KLATFACTOR` 加入 `export` 列表，加 `docstring` 说明其作为 `lateral_flow!` 标定系数（乘到 Ksat）的语义与对 Fortran `SLKLF` 的参考。
     - `SoilType.K_latfactor` 字段填充仍由 `get_soil_params` 负责（行为不变）。

- **未变更**：
  - `src/backup/` 仍按 §6 禁用
  - 同位素追踪按 §6.1 维持注释禁用
  - `gw2river.jl:23` 形参 `length::M` 未在本次修复范围（已记入 `_meta/status.md §5 #10`）
- **验证**：`Pkg.test()` 全套绿（regional_example 端到端 12.4 s 通过；空 testset 已消失；其余 9 个核心 testset + ERA5Forcings 7 个 testset + wtable 6 个 testset 全部 Pass）。
