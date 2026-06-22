// docs/theory/03_根系吸水.typ
// ==========================
// 章：根系吸水与冠层蒸腾
// 状态：骨架（待摄取）
// 来源：src/extraction.jl::extraction
// Fortran 参考：fortran/extraction.f90::extraction

#import "../preamble.typ": *

= 根系吸水与冠层蒸腾

== 1. 物理过程

`extraction` 是 ASAP 模型的冠层-土壤耦合核心：根据 Penman-Monteith 方程
计算潜在冠层蒸腾（`pet_c`）与土壤蒸发（`pet_s`），并按 Feddes 水分胁迫函数
从各土壤层提取水分。

== 2. 控制方程

待补：

- Feddes 水分胁迫函数 $f_"swp"(theta)$ @feddes1978simulation
- Penman-Monteith 蒸腾公式
- 根系活性权重 $w_k(z)$

== 3. 数值离散与求解

待补：层间权重、阈值剪枝、非活跃日自适应抑制。

== 4. 参数与单位

待补：完整参数表（见 wiki/julia/extraction-根系吸水.md §4 已成型）。

== 5. 代码实现

- 主函数：`src/extraction.jl:L45-L186`（`extraction` 全函数签名）
- 公式对照：`fortran/extraction.f90:extraction`
- 测试断言：`test/test_extraction.jl`

== 6. 局限与已知问题

- `kroot = k - 1` 语义已在 `src/extraction.jl:L65-L70` 注释澄清。
- `dθ_deep` 始终返回 0.0，地下水直接耦合尚未启用。
- `POTLEAF == POTWILT = -153.0` 未暴露为参数。