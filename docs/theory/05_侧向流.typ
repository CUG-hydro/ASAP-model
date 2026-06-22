// docs/theory/05_侧向流.typ
// ========================
// 章：侧向地下水流（D8 算法）
// 状态：骨架（待摄取）
// 来源：src/modules/lateral_flow.jl
// Fortran 参考：fortran/module_wtable.f90::LATERALFLOW4

#import "../preamble.typ": *

= 侧向地下水流

== 1. 物理过程

ASAP 使用 8 邻居（D8）算法将每层地下侧向流分配到相邻 8 个网格格点，
驱动量为水头差，落入对应高程更低邻居。

== 2. 控制方程

待补：D8 流向分配 + 水头差驱动的通量公式。

== 3. 数值离散与求解

待补：D8 索引约定、`flowdir` 工具函数、halo 边界处理。

== 4. 参数与单位

待补。

== 5. 代码实现

- 主函数：`src/modules/lateral_flow.jl`（`lateral_flow!`）
- 流向工具：`src/helper.jl`（`find_jwt`、`flowdir`）
- 公式对照：`fortran/module_wtable.f90:LATERALFLOW4`
- 测试断言：`test/wtable/test-river.jl`

== 6. 局限与已知问题

- 陆面掩膜 halo 未保护，需在调用方确保 `is+1:ie-1 / js+1:je-1` 内层循环。
- `Δt` 当前未参与计算（流量单位等价于瞬时），与 `rootdepth_main` 的
  `Δt/steps` 切片需在调用层对齐。