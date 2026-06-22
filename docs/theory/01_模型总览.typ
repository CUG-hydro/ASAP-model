// docs/theory/01_模型总览.typ
// ==========================
// 章：模型总览与主循环
// 状态：骨架（待摄取）
// 来源：src/RootDepth.jl::rootdepth_main

#import "../preamble.typ": *

= 模型总览与主循环

== 1. 物理过程

ASAP-model 是一个一维 + 网格耦合的全球陆面水文模型，每个网格格点按 6 步流程
顺序执行完整的水量平衡计算。

== 2. 控制方程

待补：水量平衡守恒方程 + 6 步主循环公式描述。

$ "ΔS" = P - "ET" - R - Q_"gw" $ <eq:water_balance>

== 3. 数值离散与求解

待补：时间步离散（`Δt = 3600 s`）、子步 `steps` 控制、CFL 约束。

== 4. 参数与单位

待补：状态变量清单与单位。

== 5. 代码实现

- 主入口：`src/RootDepth.jl:L131-L274`（`rootdepth_main` 6 步主循环）
- 类型参数化：`src/RootDepth.jl:L17-L42`（`T<:Real, V<:Vector{T}, M<:Matrix{T}, A3<:Array{T,3}`）
- 测试断言：`test/test_rootdepth.jl`

== 6. 局限与已知问题

待补。