# Wiki 索引

> 按子系统类别组织的页面清单，每行一条：标题 + 简介 + 链接。

## 概览

- [README](./README.md) — Wiki 总览、项目背景、目录树、导航入口
- [conventions](./conventions.md) — 命名约定、单位约定、文档页面模板
- [index](./index.md) — 本页（按类别组织的完整索引）
- [log](./log.md) — 摄取与变更时间线
- [ASAP-主入口](./julia/ASAP-主入口.md) — `module ASAP` 的 include 顺序、`@reexport` 与 export 列表
- [RootDepth-主算法](./julia/RootDepth-主算法.md) — `rootdepth_main` 泛型签名、主循环 6 步流程、同位素追踪遗留问题

## 土壤子系统

- [SoilParameters-土壤参数](./julia/SoilParameters-土壤参数.md) — 13 种土壤类型参数表、`SoilType` struct、Campbell 导水率 `cal_K`、凋萎点 `init_soil_param`
- [SoilInitialization-土壤分层](./julia/SoilInitialization-土壤分层.md) — CLM 指数分层 `initializesoildepth_clm` 与固定层厚 `initializesoildepth`
- [SoilFluxes-土壤水运动](./julia/SoilFluxes-土壤水运动.md) — 一维 Richards 方程 Crank-Nicolson 离散、Thomas 三对角求解、自由/受限排水边界
- [extraction-根系吸水](./julia/extraction-根系吸水.md) — Penman-Monteith 蒸腾、根系活性权重、各层水分提取 `dθ` 与水分亏缺 `watdef`

## 蒸散发子系统

- [Evapotranspiration-蒸散发](./julia/Evapotranspiration-蒸散发.md) — Priestley-Taylor、Penman-MonteITH、Shuttleworth-Wallace 双源（唯一显式 `module` 包裹）
- [Interception-截留](./julia/Interception-截留.md) — 冠层截留容量 `0.2·LAI`、穿透降水与截留蒸发分支判定

## 地下水位子系统

- [updatewtd_shallow-浅层水位](./julia/updatewtd_shallow-浅层水位.md) — 浅层 WTD 更新、`cal_factor` 深度衰减、`find_jwt` 水位层定位
- [lateral_flow-侧向地下水流](./julia/lateral_flow-侧向地下水流.md) — D8 八邻居达西通量、有效导水率按水位埋深分支
- [IsotopeTracing-同位素追踪](./julia/IsotopeTracing-同位素追踪.md) — ¹⁸O 侧向流传输、深层水位更新、4 方向含纬度修正

## 河流子系统

- [gw2river-地下水河流交换](./julia/gw2river-地下水河流交换.md) — `gw2river!` 双向通量、`moveqrf!` 下游重分配、人工分流
- [rivers_kw_flood-运动波路由](./julia/rivers_kw_flood-运动波路由.md) — Manning 公式河道-洪泛演算、非洪水态回退
- [rivers_dw_flood-扩散波路由](./julia/rivers_dw_flood-扩散波路由.md) — 扩散波隐式更新、frozen-coefficient 线性化、CFL 钳制
- [flooding-洪泛漫流](./julia/flooding-洪泛漫流.md) — 8 邻域最低高程漫流、主河道扣除、地表水守恒

## 模块聚合

- [modules-水位模块聚合](./julia/modules-水位模块聚合.md) — `src/modules/Modules.jl` 的 include 序列与 export 表（含悬空 export 注释）

## Fortran 原版

- [README](./fortran/README.md) — Fortran 源目录总览
- [module_driver](./fortran/module_driver.md) — 主驱动模块
- [module_forcings](./fortran/module_forcings.md) — 强迫数据处理
- [module_initial](./fortran/module_initial.md) — 初始化与土壤参数
- [module_io](./fortran/module_io.md) — 输入/输出
- [module_parallel](./fortran/module_parallel.md) — 并行化
- [module_rootdepth](./fortran/module_rootdepth.md) — 根系深度主算法
- [module_wtable](./fortran/module_wtable.md) — 地下水位与河流耦合

## 跨语言映射

- [julia-fortran-对照](./mapping/julia-fortran-对照.md) — 17 个 Julia 文件 ↔ 11 个 Fortran 文件 的逐函数映射
- [algorithm-索引](./mapping/algorithm-索引.md) — 按算法（Richards、Manning、Shuttleworth-Wallace 等）反向索引

## 元信息

- [_meta/status](./_meta/status.md) — 全部页面摄取状态表（已摄取 / 待复核）
- [_meta/cross-refs](./_meta/cross-refs.md) — 跨页面符号引用登记表（至少 10 条）