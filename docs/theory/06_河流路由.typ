// docs/theory/06_河流路由.typ
// =========================
// 章：河流路由（运动波 KW / 扩散波 DW）
// 状态：骨架（待摄取；现有 docs/汇流_扩散波.typ 讲稿保持不动）
// 来源：src/modules/rivers_kw_flood.jl、src/modules/rivers_dw_flood.jl

#import "../preamble.typ": *

= 河流路由

== 1. 物理过程

ASAP-model 提供两种河流路由方案：

- 运动波（KW）：忽略惯性项，适用于陡峭山区
- 扩散波（DW）：保留压强梯度项，适用于河漫滩

== 2. 控制方程

待补：从 `docs/汇流_扩散波.typ` 现有讲稿迁入（保留原文推导）。

=== 2.1 圣维南方程组

== 3. 数值离散与求解

待补：显式时间步 + Manning 公式 + 冻结系数线性化。

== 4. 参数与单位

待补：河道-洪泛区重分配规则（`riverarea·maxdepth` 阈值）。

== 5. 代码实现

- KW：`src/modules/rivers_kw_flood.jl`（`rivers_kw_flood!`）
- DW：`src/modules/rivers_dw_flood.jl`（`rivers_dw_flood!`）
- 公式对照：`fortran/module_wtable.f90:RIVERS_KW_FLOOD` / `RIVERS_DW_FLOOD`
- 测试断言：`test/wtable/test_rivers_dw_flood.jl` / `test_river_routing.jl`

== 6. 局限与已知问题

- `n=0.03` 曼宁糙率硬编码，未在参数列表暴露。
- KW 与 DW 存在功能重叠，应在缓坡态由 DW 接管。
- `rivers_kw_flood.jl` 的 `length` 形参已重命名为 `river_length`（2026-06-22 修复）。