// docs/theory/04_蒸散发.typ
// =========================
// 章：潜在蒸散发（三种 PET 方法）
// 状态：骨架（待摄取）
// 来源：src/Evapotranspiration.jl

#import "../preamble.typ": *

= 潜在蒸散发

== 1. 物理过程

ASAP-model 提供三种潜在蒸散发（PET）方法：

- Priestley-Taylor（P-T）：简化的能量平衡法 @priestley1972assessment
- Penman-Monteith（P-M）：经典组合公式
- Shuttleworth-Wallace（S-W）：冠层-土壤双源 @shuttleworth1985evaporation

== 2. 控制方程

待补：三种方法的控制方程与适用场景。

=== 2.1 Priestley-Taylor

$ "PET"_"PT" = alpha frac(Delta (R_n - G), Delta + gamma) $ <eq:pet_pt>

=== 2.2 Penman-Monteith

$ "PET"_"PM" = frac(
    Delta (R_n - G) + rho_a c_p (e_s - e_a) / r_a,
    Delta + gamma (1 + r_s / r_a)
  ) $ <eq:pet_pm>

=== 2.3 Shuttleworth-Wallace

$ "ET" = C_c "PM"_c + C_s "PM"_s $ <eq:pet_sw>

返回 12 元组覆盖 SW 四分量（水面/土壤/冠层/截留）。

== 3. 数值离散与求解

待补：能量平衡项处理、空气动力阻抗与表面阻抗的差分。

== 4. 参数与单位

待补：完整 BIOPARMS 参数表与辐射项单位换算。

== 5. 代码实现

- 主入口：`src/Evapotranspiration.jl`
- 三函数：
  - `potevap_priestly_taylor`
  - `potevap_penman_monteith`
  - `potevap_shutteworth_wallace`（拼写错误，`src/Evapotranspiration.jl:L7`）
  - `potevap_shuttleworth_wallace`（拼写修正的别名，L8）
- 测试断言：`test/test_evapotranspiration.jl`
- 调用入口：`src/RootDepth.jl:L131-L142`（步骤 ①）

== 6. 局限与已知问题

- 拼写错误 `shutteworth` 保留以避免破坏下游调用，新代码推荐使用别名。
- Priestley-Taylor 默认步长假设为日，逐小时调用需注意 `rad*24*3600*1e-6` 转换。
- `BIOPARMS` 索引 18-21 与 3-6 重复（USDA/IGBP 双分类并存），调用前需 `clamp(veg_int, 1, 31)`。