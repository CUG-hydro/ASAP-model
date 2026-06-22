# ASAP — 主入口模块

> 源文件：`src/ASAP.jl`（29 行）
> 依赖：[`Project.toml`](../../Project.toml) — `Reexport = "189a3867-3050-52da-a836-e630ba90ab69"`
> 状态：已摄取

## 1. 功能概述

`module ASAP` 是 ASAP-model Julia 版本的顶层包入口。它的职责有三：

1. **按依赖顺序 include** 所有 `src/*.jl` 子模块；
2. 把 `Evapotranspiration` 子模块 `@reexport` 给上层用户；
3. **统一 export** 主算法函数 `rootdepth_main` 与浅层水位函数 `updatewtd_shallow`、`updatewtd_qlat`。

外部调用方执行 `using ASAP` 后即可直接访问 `potevap_shutteworth_wallace`、`rootdepth_main`、`updatewtd_shallow` 等符号。

## 2. 完整源码（含行号）

```julia
"""
ASAP模型主模块
整合所有子模块，提供完整的水文模型功能
"""
module ASAP

using Reexport

include("helper.jl")
include("SoilParameters.jl")
include("Evapotranspiration.jl")
include("Interception.jl")
include("extraction.jl")
include("SoilInitialization.jl")
include("SoilFluxes.jl")
include("updatewtd_shallow.jl")
# include("updatewtd_qlat.jl")
include("RootDepth.jl")

# 导入新的水位模块
include("modules/Modules.jl")

@reexport using .Evapotranspiration


export updatewtd_shallow, updatewtd_qlat
export rootdepth_main

end # module ASAP
```

## 3. include 顺序与依赖关系

`include` 顺序即源码编译顺序，必须满足"被 include 文件用到的类型/函数在此之前已出现"。ASAP.jl 的顺序严格遵循如下依赖链：

| 顺序 | 文件 | 行号 | 依赖的上一级 |
|---|---|---|---|
| 1 | `helper.jl` | L9 | 无（导出 `find_jwt`、`flowdir`） |
| 2 | `SoilParameters.jl` | L10 | 无 |
| 3 | `Evapotranspiration.jl` | L11 | `SoilParameters`（USGS BIOPARMS 表） |
| 4 | `Interception.jl` | L12 | 无（仅基础算术） |
| 5 | `extraction.jl` | L13 | `SoilParameters`（POTWILT、θ_sat、ψsat）、`Evapotranspiration`（petfactor_s/c、λ、Δ、γ） |
| 6 | `SoilInitialization.jl` | L14 | 无 |
| 7 | `SoilFluxes.jl` | L15 | `SoilParameters`（cal_K、θ_sat、ψsat、b） |
| 8 | `updatewtd_shallow.jl` | L16 | `helper.jl`（`find_jwt`） |
| 9 | `updatewtd_qlat.jl` | L17 | **被注释掉**（备份在 `src/backup/updatewtd_qlat.jl`） |
| 10 | `RootDepth.jl` | L18 | 以上所有（含 `modules/Modules.jl`） |
| 11 | `modules/Modules.jl` | L21 | `src/modules/*.jl` 聚合 |

> **注意**：`modules/Modules.jl` 的 include 在 `RootDepth.jl` 之后，但 `RootDepth.jl` 内部并不调用水位聚合模块（仅使用 `helper.jl` 中已 include 的 `find_jwt`）。水位聚合模块为下游耦合调用方预留，独立于 `rootdepth_main` 主循环。

## 4. @reexport 与 export 列表

### 4.1 `@reexport using .Evapotranspiration`

`Evapotranspiration.jl` 是 `src/` 中**唯一显式以 `module Evapotranspiration ... end` 包裹**的文件（参见 [`Evapotranspiration-蒸散发.md`](./Evapotranspiration-蒸散发.md) 第 1 节与第 7 节）。该模块本身**没有任何 `export` 语句**——其内部 `potevap_priestly_taylor`、`potevap_penman_monteith`、`potevap_shutteworth_wallace` 都未导出。

`@reexport using .Evapotranspiration` 的作用是：让 ASAP 用户 `using ASAP` 时**等价于自动 `using Evapotranspiration`**，并把 `.Evapotranspiration` 模块的内部符号直接引入 ASAP 命名空间，从而无需前缀即可调用三个 PET 函数。`@reexport` 来自 `Reexport` 包（已在 `Project.toml` 声明）。

### 4.2 export 列表

```julia
export updatewtd_shallow, updatewtd_qlat   # L26
export rootdepth_main                       # L27
```

| 符号 | 来源 | 备注 |
|---|---|---|
| `updatewtd_shallow` | `src/updatewtd_shallow.jl` | 浅层 WTD 更新（已被 `RootDepth.jl` 调用） |
| `rootdepth_main` | `src/RootDepth.jl` | 主算法入口（详见 [`RootDepth-主算法.md`](./RootDepth-主算法.md)） |
| `read_initial` / `read_wtdnc` | `src/io/NetCDF.jl` | P0 子集（仅 NetCDF 读取，详见 §6.1） |
| `read_hourly_forcings` / `read_accumulated_forcings` / `read_snow` / `read_snow_hour` / `read_soil_temps` | `src/Forcings/ERA5.jl` | ERA5 强制场读取（详见 [`Forcings-ERA5.md`](./Forcings-ERA5.md)） |

## 5. 未显式 export 但可访问的符号

由于 `module ASAP` 内部按 `include` 顺序把每个 `src/*.jl` 的全部顶层定义都装入自己的命名空间，下列符号**不需要 export 也可被 `ASAP.xxx` 访问**：

- `find_jwt`、`flowdir`（`helper.jl` L1 export）
- `SoilType`、`get_soil_params`、`init_soil_param`、`cal_K`（`SoilParameters.jl` export）
- `initializesoildepth_clm`、`initializesoildepth`（`SoilInitialization.jl`）
- `extraction`（`extraction.jl`）
- `interception`（`Interception.jl` L1 export）
- `cal_factor`（`updatewtd_shallow.jl`，未 export，依赖模块隐式可见）
- `lateral_flow!`、`rivers_kw_flood!`、`rivers_dw_flood!`、`flooding!`、`lateral_isotope!`、`updatedeepwtable!`（`modules/Modules.jl` L14-L16 export，已 include）

## 6. 与 Fortran 对照

| Fortran | Julia | 差异 |
|---|---|---|
| `module_rootdepth.f90` 单一巨型模块 | `module ASAP` 聚合 14+ 个 include | Julia 拆分更细，Fortran 内全部子程序共享全局 COMMON 块 |
| `PROGRAM ASAP7` 主程序 | 顶层模块（非可执行文件） | Julia 通过 `example/complete_example.jl` 演示调用流程 |
| `use IsoCType` 等 Fortran 90 模块 | `using Reexport` | Julia 仅引入 1 个元包，依赖收敛 |
| `module_nrtype.f90` 数值类型 | 浮点参数化 `where T<:Real` | Julia 通过泛型统一 `Float32`/`Float64` |

## 7. 引用

- 完整源码：`src/ASAP.jl` L1-L29
- 依赖声明：`Project.toml` L8（Reexport = "189a3867-3050-52da-a836-e630ba90ab69"）
- 子模块 include 行：L9-L21
- `@reexport` 行：L23
- `export` 行：L26-L27

## 8. 已知问题与备注

1. **`updatewtd_qlat` 已从主入口移除**（✅ 2026-06-22 修复）：`src/ASAP.jl` L17 的 `# include("updatewtd_qlat.jl")` 注释保留以作历史说明；对应 export 已整体删除。备份实现仍保留在 `src/backup/updatewtd_qlat.jl`，按 §6 规定不再启用。
2. **`helper.jl` 与 `modules/Modules.jl` 重复定义**：`helper.jl` 在 `ASAP.jl` L9 被 include；`Modules.jl` L1 又把 `# include("helper.jl")` 注释掉，避免重复 include 报错，但 `find_jwt`、`flowdir` 仍由主模块路径提供。
3. **`Evapotranspiration.jl` 内部无 export**：三个 PET 函数仅通过 `@reexport` 暴露，对外部 `include("src/Evapotranspiration.jl")` 的脚本（非包用户）必须显式 `import .Evapotranspiration` 才能访问。
4. **`DataFrames` 悬空 import 已清理**（✅ 2026-06-22 修复）：`SoilInitialization.jl` L1 的 `using DataFrames` 与 `Project.toml` 中的 `DataFrames` 依赖均已移除，详见 `_meta/status.md §5 #9`。
5. **`wtable!` / `updatewtd!` 悬空 export 已清理**（✅ 2026-06-22 修复）：`src/modules/Modules.jl` L14-L16 仅导出当前活跃实现（`lateral_flow!` / `rivers_kw_flood!` / `rivers_dw_flood!` / `flooding!` / `lateral_isotope!` / `updatedeepwtable!`）；L10-L13 注释保留对 `backup/` 旧实现的替代说明。