# Modules.jl — 水位与同位素模块聚合

> 源文件：src/modules/Modules.jl:L1-L16
> 测试：test/wtable/test_watertable.jl, test/wtable/test_isotope_trouting.jl
> 状态：已摄取（2026-06-21 更新：移除 2 个悬空 export）

## 1. 功能概述

`Modules.jl` 是水位与河流-地下水耦合相关子模块的统一入口与导出表。它按依赖顺序 include 浅层/侧向/河-含水层等文件，并把对外公开的 `!` 修饰符函数集中导出，便于在 `ASAP` 主模块中直接 `using .Modules`。

## 2. 函数签名

本文件不定义函数，仅做模块聚合与 `export`。具体函数签名见：
- `lateral_flow!` — `src/modules/lateral_flow.jl`
- `lateral_isotope!` / `updatedeepwtable!` — `src/modules/Tracing/IsotopeTracing.jl`
- `rivers_kw_flood!` / `rivers_dw_flood!` / `flooding!` — 河流-地下水相互作用

## 3. 算法 / 公式

无独立算法。本文件只负责 `include` 与 `export`，不实现数值方法。

## 4. 关键变量与单位

| 项 | 含义 | 备注 |
|---|---|---|
| `Modules.jl` 包含文件 | `flooding.jl` / `lateral_flow.jl` / `gw2river.jl` / `rivers_dw_flood.jl` / `rivers_kw_flood.jl` / `Tracing/IsotopeTracing.jl` | `helper.jl` 已被注释掉 |
| `export` 列表 | 6 个 `!` 函数 | 见 L8-L11 |

## 5. 与 Fortran 对应

| Fortran 模块 | Julia 模块聚合 | 差异 |
|---|---|---|
| `WaterTableModule`（旧） | `Modules.jl` | 旧版函数 `wtable!` / `updatewtd!` 已迁移到 `backup/`，不再 export |
| `LateralFlowModule` | `lateral_flow.jl` | 直接对应 |
| `IsotopeModule` | `Tracing/IsotopeTracing.jl` | 新增子目录结构 |

## 6. 引用

- 行号：src/modules/Modules.jl:L1-L16
- `include` 序列：L2-L9
- `export` 列表：L8-L11
- 注释掉的 include：L1 `# include("helper.jl")`
- 测试入口：test/wtable/test_watertable.jl、test/wtable/test_isotope_tracing.jl

## 7. 已知问题与备注

**2026-06-21 修复**：
- `wtable!` / `updatewtd!` 已从 export 列表移除（2026-06-21）。
  - 原因：实现位于 `src/backup/`（已废弃），功能已被 `updatewtd_shallow` + `lateral_flow!` + `gw2river!` 组合替代。
- `gw2river!` / `moveqrf!` 由 `gw2river.jl` 自身导出，不在本模块 export 列表中。

其他备注：
- 河流-地下水组 export 已分别建独立 wiki 页（`gw2river-地下水河流交换.md` 等）。
- `helper.jl` 被注释，确认无外部引用后方可删除。
- `Tracing/` 子目录的 include 使用相对路径 `Tracing/IsotopeTracing.jl`，要求该子模块文件存在且自洽。
