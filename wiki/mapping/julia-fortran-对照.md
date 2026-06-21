# Julia ↔ Fortran 双向对照表

> 用途：将 ASAP-model 的 **17 个 Julia 源文件** 与 **11 个 Fortran 模块/主程序** 中的同名/同义过程做行级映射；差异、缺失项、悬空 export 全部显式标注。
> 状态：基于 `wiki/julia/*.md` 与 `wiki/fortran/*.md` 既有摄取页，并直接核对了 `src/`、`fortran/` 当前源码行号。
> 路径约定：所有路径均为绝对路径；行号以 `wc -l` 与源码实测为准。

## 0. 模块概览

| 维度 | Julia（17 文件 / 3 文件夹） | Fortran（11 模块） |
|---|---|---|
| 主入口 | `src/ASAP.jl`（`module ASAP`，L5-L27） | `fortran/module_driver.f90:1 program driver` |
| 物理过程实现 | `src/*.jl`（12 文件）+ `src/modules/*.jl`（6 文件，含 `Tracing/`） | `fortran/module_rootdepth.f90`、`module_wtable.f90`、`soilfluxes.f90`、`module_initial.f90` |
| 模块聚合 | `src/modules/Modules.jl:1-L14` | `module_rootdepth`、`module_wtable`、`module_initial`（独立 module） |
| 工具 / 插值 | `src/helper.jl`（`find_jwt`、`flowdir`） | `fortran/interp_lib.f90`（`TRNCL1/2`、`INTRP`、`GDTOST2/3`、`HTINTCP` 等 11 个子程序） |
| 并行 / I/O / 强迫 | 无直接对应（单进程设计） | `module_parallel.f90`、`module_io.f90`、`module_forcings.f90` |
| 数据类型 | `module_nrtype.f90`（无对应） | 仅 Fortran |

## 1. 物理过程双向对照

### 1.1 冠层界面（截留 / PET / 根系吸水）

| 物理过程 | Julia 函数 | Julia 位置 | Fortran 子程序 | Fortran 位置 | 差异说明 |
|---|---|---|---|---|---|
| Priestley-Taylor PET | `potevap_priestly_taylor` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/Evapotranspiration.jl:70-L101` | `POTEVAP_Priestly_Taylor` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/fortran/module_rootdepth.f90:546` | Julia 用局部 `CP_LOCAL` 覆盖模块常量；Fortran 取 `(i,j)` 网格索引 |
| Penman-Monteith PET | `potevap_penman_monteith` | `src/Evapotranspiration.jl:103-L166` | `POTEVAP_Penman_Monteith` | `fortran/module_rootdepth.f90:569` | 参数顺序不同：Fortran 先传湿度，Julia 用比湿+压强计算 `e` |
| Shuttleworth-Wallace PET | `potevap_shutteworth_wallace` | `src/Evapotranspiration.jl:168-L288` | `POTEVAP_Shutteworth_Wallace` | `fortran/module_rootdepth.f90:653` | Julia 返回 12 元组（含 `Δ, γ, λ, ra_a, ra_c, rs_c, R_a, R_s, pet_s, pet_c, pet_w, pet_i`）；Fortran 仅返回 4 个分量 |
| 生物物理参数表 | `const BIOPARMS` / `RL` / `ZDIS` | `src/Evapotranspiration.jl:15-L43` | 模块级 `RL`/`ZDIS` 数组 | `fortran/module_rootdepth.f90`（隐式，含 USGS 30 类） | Julia 31 类（USDA+IGBP 重复），Fortran 30 类 |
| 冠层截留 | `interception` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/Interception.jl:18-L54` | `INTERCEPTION` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/fortran/module_rootdepth.f90:277` | `minpprate` 在 Julia 文档列出但函数体未读取，存在**接口不一致**风险；Fortran 全局参数判定阈值 |
| 根系吸水与 PET | `extraction` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/extraction.jl:36-L187` | `EXTRACTION` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/fortran/module_rootdepth.f90:307` | Julia 引入 `fdepth` 深度因子；新增 `maxeasy×0.001` 阈值剪枝（L99）与非活跃层抑制；`kroot` 返回 `k-1` 而非 `k`（待复核） |
| 主算法时步入口 | `rootdepth_main` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/RootDepth.jl:90-L278` | `ROOTDEPTH` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/fortran/module_rootdepth.f90:32` | Julia 内部串联：`PET → INTERCEPTION → EXTRACTION → SOILFLUXES → UPDATESHALLOWWTD/UPDATEWTDQLAT`；Fortran 在 `ROOTDEPTH` 一次性调度所有子程序 |

### 1.2 土壤水运动（Richards / 三对角 / 边界）

| 物理过程 | Julia 函数 | Julia 位置 | Fortran 子程序 | Fortran 位置 | 差异说明 |
|---|---|---|---|---|---|
| 1D Richards 求解 | `soilfluxes` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/SoilFluxes.jl:35-L376` | `SOILFLUXES` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/fortran/soilfluxes.f90:2` / `fortran/module_rootdepth.f90:994` | Julia 用 kwargs（`freedrain::Bool`）替代 Fortran 整型标志；返回元组替代 COMMON 块；删除了氧18同位素过程（末 30 行已注释） |
| Thomas 三对角算法 | `tridag!` | `src/SoilFluxes.jl:377-L401` | `tridag` | `fortran/module_rootdepth.f90:1868` | Julia 原地修改右端向量 `r` 与解 `u`；省略 Fortran 的 `bet` 工作数组（用 `gam` 替代） |
| Campbell 导水率 `K(θ)` | `cal_K` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/SoilParameters.jl:94-L96` | `khyd` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/fortran/module_rootdepth.f90:1906` | 公式一致；Julia 显式 `Float64` 类型；Fortran 通过 `(smoi, nsoil)` 二参 |
| 含水量限幅 | 内联于 `soilfluxes` | `src/SoilFluxes.jl:300-L328` | `theta_limits`（已内联） | `fortran/soilfluxes.f90` 中段 | Julia 在主循环末尾统一校正超饱和/过干；Fortran 同名段在 `SOILFLUXES` 中段 |
| 入渗能力 | `Imax = Ksat * dt` 内联 | `src/SoilFluxes.jl:26-L42` | `infil_capacity`（已内联） | `fortran/soilfluxes.f90` | 不再单独子程序化；逻辑相同 |
| Green-Ampt 改进 | `g = 16/(π+4)` 系数 | `docs/下渗_黄土高原.typ` 对应 | 无独立子程序 | `fortran/soilfluxes.f90` | Fortran 当前**未实现**椭圆面积修正；Julia 在 `SoilFluxes.jl` 内联实现 |

### 1.3 土壤参数化与初始化

| 物理过程 | Julia 函数 / 常量 | Julia 位置 | Fortran 子程序 / 数组 | Fortran 位置 | 差异说明 |
|---|---|---|---|---|---|
| 13 类土壤参数表 | `const θSAT/SOILCP/SLBS/KSAT/ΨSAT/KLATFACTOR` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/SoilParameters.jl:25-L31` | 模块级 `θSAT/SOILCP/SLBS/KSAT/ΨSAT/KLATFACTOR` | `fortran/module_rootdepth.f90`（共享 USE 段） | Julia 使用 Unicode 标识符；Fortran 数组索引从 1，值一致 |
| `SoilType` 数据结构 | `struct SoilType` | `src/SoilParameters.jl:13-L21` | 无对应（Fortran 直接用裸数组） | — | Julia 提供类型化包装 |
| 土壤参数查询 | `get_soil_params` | `src/SoilParameters.jl:43-L57` | 直接数组访问 `slmsts(nsoil)` 等 | `fortran/module_rootdepth.f90` | Julia 返回 `SoilType`；Fortran 直接按数组下标取 |
| 凋萎点计算 | `init_soil_param` | `src/SoilParameters.jl:68-L79` | `init_soil_param`（私有 FUNCTION） | `fortran/module_rootdepth.f90:1891` | 公式相同；Julia 返回 `(fieldcp, θ_wilt)` 元组；当前 `fieldcp` 始终为 0 |
| CLM 指数分层 | `initializesoildepth_clm` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/SoilInitialization.jl:26-L52` | `INITIALIZESOILDEPTHCLM` | `fortran/module_rootdepth.f90:937` | 公式一致；Julia 显式生成 `dz` 与 `slz`；存在 `slz` 多 1 项的 bug |
| 固定厚度分层 | `initializesoildepth` | `src/SoilInitialization.jl:65-L80` | `INITIALIZESOILDEPTH` | `fortran/module_rootdepth.f90:973` | Julia 显式校验 `nzg ≤ length(DZ2)` |
| 平衡含水量初始化 | **无对应** | — | `EQSOILMOISTUREtheor` / `EQSOILMOISTUREiter` / `SOILFLUXES_EQSMOI` | `fortran/module_initial.f90:141 / 321 / 444` | Fortran 使用 Brent 求根（`zbrent`、`func1/2/3`）逐层求解，Julia 未移植 |
| Brent 求根 | **无对应** | — | `zbrent` / `func` / `func1/2/3` | `fortran/module_initial.f90:638 / 621 / 734+ / 744+` | 仅 Fortran；Julia 用简单插值替代 |
| `wtd → jwt` 查找 | `find_jwt` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/helper.jl:4-L29` | 内联查找（在 `SOILFLUXES` 等多处） | `fortran/soilfluxes.f90` | Julia 抽成工具函数；Fortran 多处复制 |

### 1.4 地下水位（浅层 / 深层 / 侧向）

| 物理过程 | Julia 函数 | Julia 位置 | Fortran 子程序 | Fortran 位置 | 差异说明 |
|---|---|---|---|---|---|
| 浅层水位更新 | `updatewtd_shallow` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/updatewtd_shallow.jl:26-L102` | `UPDATESHALLOWWTD` | `fortran/module_rootdepth.f90:1603` | 命名采用小写下划线；内联 `find_jwt`；输出改为多返回值 |
| 深度因子工具 | `cal_factor` | `src/updatewtd_shallow.jl:2-L25` | 内联 `slmsts(k) = ...` | `fortran/module_rootdepth.f90` | Julia 显式导出工具；Fortran 内联 |
| 深层水位更新 | `updatedeepwtable!` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/modules/Tracing/IsotopeTracing.jl:197-L279` | `UPDATEDEEPWTABLE` | `fortran/module_wtable.f90:185` | Julia 引入 `bottomflux` 显式回写 |
| 简化版深层更新 | `updatewtd_simple` | `src/modules/Tracing/IsotopeTracing.jl:280-L337` | **无对应** | — | 新增简化版本，原 Fortran 未单列子程序 |
| 8 方向侧向地下水流（D8） | `lateral_flow!` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/modules/lateral_flow.jl:8-L58` | `LATERAL` / `LATERALFLOW` / `LATERALFLOW4` | `fortran/module_wtable.f90:138 / 269 / 344` | Julia 8 邻居 + 角度因子；Fortran 旧版 `LATERAL` 5 邻居，`LATERALFLOW4` 4 邻居加 O18 |
| 8 方向侧向 + O18 | `lateral_isotope!` | `src/modules/Tracing/IsotopeTracing.jl:20-L72` | `LATERALFLOW4` | `fortran/module_wtable.f90:344` | 4 方向 + 纬度 cos 修正，O18 同步 |
| 侧向 + 深层 + O18 复合 | `lateralflow_with_isotope!` | `src/modules/Tracing/IsotopeTracing.jl:73-L196` | 内联于 `LATERALFLOW4` 流程 | `fortran/module_wtable.f90:344-L467` | Julia 将复合流程封装为单函数；Fortran 由 `WTABLE` 调度 |
| D8 流向编码 | `flowdir` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/helper.jl:32-L55` | `FLOWDIR` | `fortran/module_wtable.f90:1384` | 工具函数；逻辑一致 |
| 地下水补给更新 | **无对应**（仅在 `updatewtd_qlat` 历史实现中） | `src/backup/updatewtd_qlat.jl` | `UPDATEWTDQLAT` | `fortran/module_rootdepth.f90:1725` | Julia 当前实现迁入 `backup/`，活跃代码未引用 |
| 主 `wtable!` 调度 | **无对应**（悬空 export） | `src/modules/Modules.jl:11` | `WTABLE` | `fortran/module_wtable.f90:13` | Julia 中 `wtable!` / `updatewtd!` 已迁移到 `backup/`，`Modules.jl` 仍 export，会触发 `UndefVarError` |
| 深层水位补给 | `UPDATEDEEPWTABLE` 简化版 | — | `UPDATEWTD`（以 `totwater = qlat + qspring + rech` 驱动） | `fortran/module_wtable.f90:468` | Fortran 同时维护 `UPDATEDEEPWTABLE` 与 `UPDATEWTD` 两个版本 |

### 1.5 河流-地下水交换

| 物理过程 | Julia 函数 | Julia 位置 | Fortran 子程序 | Fortran 位置 | 差异说明 |
|---|---|---|---|---|---|
| 地下水→河流通量 | `gw2river!` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/modules/gw2river.jl:19-L72` | `GW2RIVER` | `fortran/module_wtable.f90:744` | 同名移植，类型参数化；保留 `K_rb` 指数衰减与 50 mm/day 上限（Miguez-Macho 2007） |
| 河流 `qrf` 下游重分配 | `moveqrf!` | `src/modules/gw2river.jl:73-L97` | `MOVEQRF` | `fortran/module_wtable.f90:1325` | 直接移植小河流 → 下游重分配逻辑 |
| 硬编码人工分流 | `apply_specific_diversions` | `src/modules/gw2river.jl:98-L125` | 内联硬编码分流（Taquari 河等） | `fortran/module_wtable.f90`（RIVERS_KW_FLOOD 内联） | 抽成独立函数，按 `(i,j)` 索引加减下游格点流量 |
| 运动波河流路由 | `rivers_kw_flood!` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/modules/rivers_kw_flood.jl:11-L101` | `RIVERS_KW_FLOOD` | `fortran/module_wtable.f90:800` | 接口完全对齐；将 `deltat`/`dtlr` 拆为 `Δt`/`δt` 参数语义化 |
| 扩散波河流路由 | `rivers_dw_flood!` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/modules/rivers_dw_flood.jl:8-L100` | `RIVERS_DW_FLOOD` | `fortran/module_wtable.f90:996` | 接口完全对齐；`g0` 提升为模块级 `const` |
| 洪泛区漫流 | `flooding!` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/modules/flooding.jl:6-L70` | `FLOODING` | `fortran/module_wtable.f90:1226` | 主体移植；`ntsplit` 显式赋 1（Fortran 中可外部传入）；保留 D8 邻格遍历顺序 |
| 河床连通性 | `apply_specific_diversions` + `moveqrf!` 联合 | `src/modules/gw2river.jl` | 通过 `fd`/`bfd` 在 `GW2RIVER`、`RIVERS_KW_FLOOD` 内隐式判断 | `fortran/module_wtable.f90:744 / 800` | Julia 抽成显式函数调用；Fortran 在子程序内分支处理 |

### 1.6 同位素（¹⁸O）追踪

| 物理过程 | Julia 函数 | Julia 位置 | Fortran 子程序 | Fortran 位置 | 差异说明 |
|---|---|---|---|---|---|
| 侧向流同位素 | `lateral_isotope!` | `src/modules/Tracing/IsotopeTracing.jl:20-L72` | `LATERALFLOW4`（参数 `o18wtd`） | `fortran/module_wtable.f90:344` | 4 方向而非 8 方向；累积数组 `*sum` 每次调用末尾乘以 `1e3`（m→mm） |
| 土壤层 ¹⁸O 对流 | 内联于 `SoilFluxes.jl`（已注释） | `src/SoilFluxes.jl` 末 30 行（注释） | `SOILFLUXES`（含 `transpo18`/`o18`） | `fortran/module_rootdepth.f90:994` | Julia 当前**未实现**；仅保留注释骨架 |
| 同位素平衡分馏 | **无对应** | — | 内联 `alpha = 1/exp(1137/T²-0.4156/T-0.0020667)`（Majoube 公式） | `fortran/module_rootdepth.f90` | Fortran 用 Majoube 经验公式 |
| O18 降水气候态读取 | **无对应** | — | `READO18CLIM` | `fortran/module_forcings.f90:1174` | 仅 Fortran |

### 1.7 强迫 / I/O / 并行（Fortran 独占）

| 物理过程 | Julia 函数 | Julia 位置 | Fortran 子程序 | Fortran 位置 | 差异说明 |
|---|---|---|---|---|---|
| 时步主循环 | `rootdepth_main` | `src/RootDepth.jl:90` | `program driver` | `fortran/module_driver.f90:1` | Julia 单进程；Fortran MPI 并行 |
| 小时 ERA5 强迫 | **无对应** | — | `READFORCINGS` | `fortran/module_forcings.f90:14` | Fortran 专用 |
| 累积辐射/降水强迫 | **无对应** | — | `READFORCINGSACC` | `fortran/module_forcings.f90:337` | Fortran 专用 |
| 净辐射分通道 | **无对应** | — | `READFORCINGSRNET` | `fortran/module_forcings.f90:404` | Fortran 专用 |
| 融雪强迫 | **无对应** | — | `READFORCINGSSNOW` | `fortran/module_forcings.f90:557` | Fortran 专用 |
| 土壤温度 + 冰冻因子 | **无对应** | — | `READFORCINGSSOILT` | `fortran/module_forcings.f90:692` | Fortran 专用；Julia 中 `icefactor` 由调用方外部传入 |
| ERA5 高程 | **无对应** | — | `READTOPOERA5` | `fortran/module_forcings.f90:792` | Fortran 专用 |
| LAI 读取 | **无对应** | — | `READLAI` / `READLAICLIM` / `READLAICHINA` | `fortran/module_forcings.f90:913 / 1028 / 1100` | Fortran 专用 |
| NetCDF 静态场读取 | **无对应** | — | `READINITIAL` / `READWTDNC` / `READVEG` / `READFLOWDIRECTION` / `READRIVERPARAMETERS` 等 12 个 | `fortran/module_io.f90:10 / 309 / 378 / 591 / 693` | Fortran 专用 |
| Restart / History 读取 | **无对应** | — | `READHISTORYNC` / `READHISTORYVARNC` / `READHISTORYVAR3DNC` / `READHISTORYNCblock` | `fortran/module_io.f90:823 / 1050 / 1173 / 1245` | Fortran 专用 |
| NetCDF 输出 | **无对应** | — | `WRITEOUTPUTNC_par` / `WRITEHISTORYNC_par` / `WRITEOUTPUTNC_INFIL_par` / `WRITEOUTPUTNC_DAILY_par` | `fortran/module_io.f90:2115 / 2415 / 3311 / 3095` | Fortran 专用 |
| MPI 域分解 | **无对应** | — | `INITIALIZEDOMAIN` / `DIVIDEDOMAIN` | `fortran/module_parallel.f90:36 / 252` | Fortran 专用 |
| Halo 通信 | **无对应** | — | `SENDBORDERS4` / `SENDBORDERS4blocking` / `SENDBORDERSFLOOD4` | `fortran/module_parallel.f90:558 / 522 / 601` | Fortran 专用 |
| 自定义 MPI 数据类型 | **无对应** | — | `domblock` / `arraysection` / `updownhalo` / `leftrighthalo` 等 11 个 | `fortran/module_parallel.f90`（CONTAINS 前） | Fortran 专用 |
| 网格插值（双线性 / 三次） | **无对应** | — | `gdtost2` / `gdtost3` / `TRNCL1/2` / `INTRP` / `HTINT2` 等 11 个 | `fortran/interp_lib.f90` | Fortran 专用 |
| 饱和水汽混合比（Tetens） | **无对应** | — | `rslf` | `fortran/module_forcings.f90`（私有） | Julia 用浮点公式等价内联 |

### 1.8 模块聚合

| 物理过程 | Julia 文件 | Fortran 模块 | 差异 |
|---|---|---|---|
| 水位模块统一入口 | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/modules/Modules.jl:1-L14` | 无独立 module；功能由 `module_wtable` / `module_rootdepth` 承担 | Julia 集中 `include` + `export`；Fortran 直接 `use` 各 module |
| 11 个 `!` 函数 export | `lateral_flow!`、`lateral_isotope!`、`updatedeepwtable!`、`wtable!`（悬空）、`updatewtd!`（悬空）、`gw2river!`、`rivers_kw_flood!`、`rivers_dw_flood!`、`flooding!`、`moveqrf!`、`apply_specific_diversions` | 无对应概念（Fortran 子程序均 PUBLIC） | Julia 通过 `!` 后缀表示就地修改；Fortran 子程序均默认 PUBLIC |
| 类型稳定性注释 | 无对应 | 无对应 | Julia 类型参数化 `M<:Matrix{T}` |

## 2. 未移植 / 未实现清单（Julia 缺失）

以下 Fortran 功能在当前 Julia 代码库中**没有对应实现**或仅在 `backup/` 中：

| Fortran 来源 | 功能描述 | Julia 状态 |
|---|---|---|
| `module_parallel.f90` 全部 | MPI 域分解、halo 通信、自定义数据类型 | **无对应**（Julia 单进程） |
| `module_io.f90` 全部 | NetCDF 输入/输出、history、restart | **无对应** |
| `module_forcings.f90` 全部 | ERA5 / LAI / O18 / 融雪 / 土壤温度读取 | **无对应** |
| `interp_lib.f90` 全部 | 网格插值工具 | **无对应** |
| `module_initial.f90:141-444` | `EQSOILMOISTUREtheor` / `EQSOILMOISTUREiter` / `SOILFLUXES_EQSMOI` 平衡含水量初始化 | **无对应**（Julia 跳过此步骤） |
| `module_initial.f90:621-826` | Brent 求根 `zbrent` / `func1/2/3` | **无对应** |
| `module_rootdepth.f90:1891 init_soil_param`（内部） | Fortran 内部初始化 | Julia 公开同名 `init_soil_param` 但公式可能不同 |
| `module_rootdepth.f90` 内嵌 `alpha` 平衡分馏 | Majoube 公式 | **无对应**（Julia 土壤层 ¹⁸O 已注释） |
| `soilfluxes.f90` Green-Ampt 椭圆修正 | 王文焰 2003 系数 16/(π+4) | **无对应**（仅在 `docs/下渗_黄土高原.typ` 中作为方案列出） |
| `module_wtable.f90:13 WTABLE` | 侧向流 + 深层 + 河流全流程调度 | Julia 中 `wtable!` 悬空 export（迁移至 `backup/`） |

## 3. 悬空 / 风险项

| 项 | 位置 | 风险 |
|---|---|---|
| `minpprate` 未读 | `src/Interception.jl:5`（签名/函数体） | 调用方若未提供会触发 `UndefVarError` |
| `wtable!` / `updatewtd!` 悬空 export | `src/modules/Modules.jl:11` | 调用会触发 `UndefVarError` |
| `updatewtd_qlat.jl` | `src/backup/` | 活跃代码未引用；活跃版本仅 `updatewtd_shallow` |
| `helper.jl` 在 `Modules.jl` 中已注释 | `src/modules/Modules.jl:1` | `find_jwt`、`flowdir` 通过 `ASAP.jl:9` 顶层 include 仍可用 |
| `slz` 多 1 项 bug | `src/SoilInitialization.jl:23-L24, L48` | `slz` 长度 `nzg+1`，`dz` 仅 `nzg` 项 |
| `fieldcp` 始终为 0 | `src/SoilParameters.jl:71` | 占位实现 |
| `kroot = k - 1` 待复核 | `src/extraction.jl:67-L72` 注释 | 可能跳过根区最浅层 |
| `icefactor` 全局注入 | `src/SoilFluxes.jl:52`（裸引用 `f_ice`） | 来自父模块作用域，建议统一为 `icefactor[k]` |
| `DataFrames` 悬空 import | `src/SoilInitialization.jl:1` | 当前实现未使用 |

## 4. 引用

- Julia 源文件根目录：`/mnt/z/GitHub/jl-pkgs/ASAP-model/src/`
  - `ASAP.jl`、`Evapotranspiration.jl`、`extraction.jl`、`helper.jl`、`Interception.jl`、`RootDepth.jl`、`SoilFluxes.jl`、`SoilInitialization.jl`、`SoilParameters.jl`、`updatewtd_shallow.jl`
  - `modules/Modules.jl`、`modules/flooding.jl`、`modules/gw2river.jl`、`modules/lateral_flow.jl`、`modules/rivers_dw_flood.jl`、`modules/rivers_kw_flood.jl`、`modules/Tracing/IsotopeTracing.jl`
- Fortran 源文件根目录：`/mnt/z/GitHub/jl-pkgs/ASAP-model/fortran/`
  - `module_driver.f90`、`module_forcings.f90`、`module_initial.f90`、`module_io.f90`、`module_nrtype.f90`、`module_parallel.f90`、`module_rootdepth.f90`、`module_wtable.f90`、`soilfluxes.f90`、`interp_lib.f90`、`main.f90`
- 既有 Julia wiki 页面（14 篇）：`/mnt/z/GitHub/jl-pkgs/ASAP-model/wiki/julia/`
- 既有 Fortran wiki 页面（7 篇）：`/mnt/z/GitHub/jl-pkgs/ASAP-model/wiki/fortran/`