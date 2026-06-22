# ASAP-model Wiki

> 项目根路径：`/mnt/z/GitHub/jl-pkgs/ASAP-model`
> Wiki 版本：v1.0（2026-06-21 初始化）
> 维护者：Dongdong Kong

## 项目背景

ASAP（Agricultural Systems Analysis and Prediction）模型最初以 Fortran `module_rootdepth.f90` 实现，是一套面向农业生态系统的陆面水文模型，覆盖土壤水运动、根系吸水、蒸散发、地下水位动态与河流-地下水耦合等关键过程。本仓库是其 **Julia 1.x 翻译版本**，将原始单一 Fortran 模块拆分为多个独立、可测试、可复用的子模块（参数、初始化、PET、截留、Richards 求解、侧向流、洪水路由等），并采用泛型类型参数实现 Float64 / Float32 双精度支持。

本 Wiki 是模型的**完整技术档案**：每一份源文件都对应一张知识卡片，详述其函数签名、关键公式、跨语言映射、已知缺陷与引用行号。读者可按子系统（土壤/蒸散发/地下水位/河流）或按 Fortran → Julia 映射两种路径阅读。

## Wiki 目录树

```
wiki/
├── README.md                          # 本文件
├── index.md                           # 全部页面索引
├── log.md                             # 摄取/变更时间线
├── conventions.md                     # 命名、单位、文档模板
├── _meta/
│   ├── status.md                      # 摄取状态表
│   └── cross-refs.md                  # 跨页面引用登记
├── julia/
│   ├── ASAP-主入口.md                  # 主模块 include/export
│   ├── RootDepth-主算法.md             # rootdepth_main 泛型主循环
│   ├── SoilParameters-土壤参数.md
│   ├── SoilInitialization-土壤分层.md
│   ├── SoilFluxes-土壤水运动.md
│   ├── extraction-根系吸水.md
│   ├── Evapotranspiration-蒸散发.md
│   ├── Interception-截留.md
│   ├── updatewtd_shallow-浅层水位.md
│   ├── lateral_flow-侧向地下水流.md
│   ├── gw2river-地下水河流交换.md
│   ├── rivers_kw_flood-运动波路由.md
│   ├── rivers_dw_flood-扩散波路由.md
│   ├── flooding-洪泛漫流.md
│   ├── IsotopeTracing-同位素追踪.md
│   ├── io-NetCDF.md
│   ├── Forcings-ERA5.md
│   ├── modules-水位模块聚合.md
│   └── example-regional.md             # example/regional_example.jl 端到端脚本
├── fortran/                           # Fortran 原版 8 个模块
└── mapping/                           # 跨语言对照表
```

## 快速导航

- **新读者**：先读 [`index.md`](./index.md) 与 [`conventions.md`](./conventions.md)，再按需深入子系统页面。
- **核心算法入口**：[`julia/RootDepth-主算法.md`](./julia/RootDepth-主算法.md) 描述 `rootdepth_main` 泛型主循环。
- **区域端到端示例**：[`julia/example-regional.md`](./julia/example-regional.md) 描述 `example/regional_example.jl` 的 6 步主流程、mock 数据集合成、ERA5 多日滚动与最终输出。
- **移植对照**：[`fortran/`](./fortran/) 子目录对应 8 个原始 Fortran 模块；`mapping/` 提供 Fortran→Julia 行级映射。
- **缺陷追踪**：[`_meta/status.md`](./_meta/status.md) 列出每个页面的"已摄取 / 待复核"状态及悬空 import、悬空 export 等遗留问题。
- **跨页引用**：[`_meta/cross-refs.md`](./_meta/cross-refs.md) 记录至少 10 条符号级别的跨页面引用，便于追溯公式与变量在多个模块间的传播。

## 当前覆盖范围

- 17 个 Julia 源文件（`src/*.jl` + `src/modules/*.jl` + `src/modules/Tracing/*.jl`）
- 1 个区域应用示例（`example/regional_example.jl`，含 mock / 真实两模式 + 多日滚动）
- 11 个 Fortran 源文件（`fortran/*.f90`）
- 阶段 A 顶层结构 8 个文件（README / index / log / conventions / ASAP-主入口 / RootDepth-主算法 / status / cross-refs）

更多历史变更见 [`log.md`](./log.md)。