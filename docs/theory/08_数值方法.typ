// docs/theory/08_数值方法.typ
// ==========================
// 章：数值方法（综述）
// 状态：骨架（待摄取）

#import "../preamble.typ": *

= 数值方法

== 1. 概述

本章汇总 ASAP-model 各物理过程的离散格式与求解算法，侧重稳定性、精度与守恒性分析。

== 2. 时间离散

- 主循环：`src/RootDepth.jl:L131-L274`（6 步顺序执行）
- 子步：`steps = 1.0`（默认），`Δt/steps` 控制 CFL

== 3. 空间离散

- 1D 土壤水：中心差分（详见 `theory/02_Richards方程.typ`）
- 2D 侧向流：D8 流向分配（详见 `theory/05_侧向流.typ`）
- 1D 河流：显式时间步 + 河道-洪泛区分摊（详见 `theory/06_河流路由.typ`）

== 4. 求解器

- Thomas 三对角：`src/SoilFluxes.jl`（Crank-Nicolson 解 Richards）
- `zbrent`：标量根求解（`src/ASAP.jl` 已 export）

== 5. 类型参数化

`src/RootDepth.jl:L17-L42` 的泛型签名支持 `Float64`/`Float32` 双精度运行：

```julia
function rootdepth_main(...) where {
  T<:Real, V<:Vector{T}, M<:Matrix{T}, A3<:Array{T,3}
}
```

== 6. 局限与已知问题

- 类型参数化对 `Int8`/`Int16`/`Int` 离散数组保持具体类型（`icefactor`、`infilk`），
  调用方需保证类型匹配。
- 数值精度切换测试在 `test/test_rootdepth.jl` 断言 `Float32`/`Float64` 一致至 `1e-4`。
- 数值守恒测试在 `test/test_extraction.jl:L118` 验证 `sum(dθ) ≈ pet_c*1e-3`。