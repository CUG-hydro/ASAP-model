// docs/theory/07_同位素.typ
// ========================
// 章：¹⁸O 同位素追踪
// 状态：骨架（按 CLAUDE.md §6.1 暂缓在主循环启用；本章仅描述理论）
// 来源：src/modules/Tracing/IsotopeTracing.jl、src/SoilFluxes.jl
// Fortran 参考：fortran/module_rootdepth.f90::Majoube 平衡分馏 α 公式

#import "../preamble.typ": *

= ¹⁸O 同位素追踪

== 1. 物理过程

待补：¹⁸O 作为水文示踪剂的基本原理、Majoube 平衡分馏 @majoube1971oxygen。

== 2. 控制方程

待补：

- α 平衡分馏公式
- δ¹⁸O 的 Craig-Gordon 模型
- 蒸发富集方程

== 3. 数值离散与求解

待补：三对角组装 ¹⁸O（已被注释，按 §6.1 暂缓）。

== 4. 参数与单位

待补：δ¹⁸O 单位（‰ VSMOW）。

== 5. 代码实现

- 子模块：`src/modules/Tracing/IsotopeTracing.jl`
  - `lateral_isotope!`
  - `lateralflow_with_isotope!`
  - `updatedeepwtable!`
  - `updatewtd_simple`
- 注释段：`src/SoilFluxes.jl:L347-L351`（30 行死代码，三对角组装）
- 单元测试保留（独立验证）

== 6. 局限与已知问题

- 按 `CLAUDE.md §6.1` 规定，*当前阶段不实现 ¹⁸O 串联*：
  - `src/SoilFluxes.jl` 末段被注释的 ¹⁸O 三对角组装维持禁用；
  - `src/RootDepth.jl` 禁止调用 `lateral_isotope!` / `updatedeepwtable!`；
  - `fortran/module_rootdepth.f90` 的 Majoube α 公式禁止实现。
- 解锁条件：用户明确重新授权 + 在 `wiki/log.md` 追加
  `#raw[[YYYY-MM-DD] enable | 范围]` 条目。
- `IsotopeTracing.jl` 4 方向不含对角线，仅保留单测。