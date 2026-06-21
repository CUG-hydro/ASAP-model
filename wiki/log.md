# 变更日志（log）

> 本文件记录 Wiki 的摄取、复核与重大变更，按日期倒序排列。

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