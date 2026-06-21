# Modules.jl — 水位与同位素模块聚合

> 源文件：src/modules/Modules.jl:L1-L20
> 测试：test/wtable/test_watertable.jl, test/wtable/test_isotope_tracing.jl
> 状态：已摄取（含悬空 export 备注）

## 1. 功能概述

`Modules.jl` 是水位与河流-地下水耦合相关子模块的统一入口与导出表。它按依赖顺序 include 浅层/侧向/河-含水层等文件，并把对外公开的 `!` 修饰符函数集中导出，便于在 `ASAP` 主模块中直接 `using .Modules`。

## 2. 函数签名

本文件不定义函数，仅做模块聚合与 `export`。具体函数签名见：
- `lateral_flow!` — `src/modules/lateral_flow.jl`
- `lateral_isotope!` / `updatedeepwtable!` — `src/modules/Tracing/IsotopeTracing.jl`
- `gw2river!` — `src/modules/gw2river.jl`
- `rivers_kw_flood!` / `rivers_dw_flood!` / `flooding!` / `moveqrf!` — 河流-地下水相互作用
- `wtable!` / `updatewtd!` — 悬空 export（见下）

## 3. 算法 / 公式

无独立算法。本文件只负责 `include` 与 `export`，不实现数值方法。

## 4. 关键变量与单位

| 项 | 含义 | 备注 |
|---|---|---|
| `Modules.jl` 包含文件 | `flooding.jl` / `lateral_flow.jl` / `gw2river.jl` / `rivers_dw_flood.jl` / `rivers_kw_flood.jl` / `Tracing/IsotopeTracing.jl` | `helper.jl` 已被注释掉 |
| `export` 列表 | 11 个 `!` 函数 | 见 L8-L11 |

## 5. 与 Fortran 对应

| Fortran 模块 | Julia 模块聚合 | 差异 |
|---|---|---|
| `WaterTableModule`（旧） | `Modules.jl` | 旧版函数 `wtable!` / `updatewtd!` 已迁移到 `backup/`，当前未 include |
| `LateralFlowModule` | `lateral_flow.jl` | 直接对应 |
| `IsotopeModule` | `Tracing/IsotopeTracing.jl` | 新增子目录结构 |

## 6. 引用

- 行号：src/modules/Modules.jl:L1-L20
- `include` 序列：L2-L9
- `export` 列表：L8-L11
- 注释掉的 include：L1 `# include("helper.jl")`
- 测试入口：test/wtable/test_watertable.jl、test/wtable/test_isotope_tracing.jl

## 7. 已知问题与备注

**悬空 export**（重要）：
- `wtable!` 与 `updatewtd!` 在 L8 被 `export`，但当前 `Modules.jl` 中**没有对应的 include**，其实现位于 `backup/` 目录。若调用方依赖这两个名字会触发 `UndefVarError`。
- 修复建议：要么将 `backup/` 中的同名文件 include 进来（注意接口兼容），要么从 `export` 列表中移除这两个名字直至正式复活。

其他备注：
- 河流-地下水组 4 个 export 未在本 wiki 单独建页，需后续补页。
- `helper.jl` 被注释，确认无外部引用后方可删除。
- `Tracing/` 子目录的 include 使用相对路径 `Tracing/IsotopeTracing.jl`，要求该子模块文件存在且自洽。
