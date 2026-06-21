# 变更日志（log）

> 本文件记录 Wiki 的摄取、复核与重大变更，按日期倒序排列。

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