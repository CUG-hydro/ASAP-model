// docs/manual/05_调参指南.typ
// ==========================
// 章：调参指南
// 状态：骨架（待摄取）

#import "../preamble.typ": *

= 调参指南

== 1. 目标

为不同气候带与下垫面条件提供 ASAP-model 参数调整建议。

== 2. 前置条件

已熟悉 `theory/01_模型总览.typ` 主循环与 `manual/04_输出诊断.typ` 输出字段。

== 3. 操作步骤

+ 土壤参数（`src/SoilParameters.jl`）：
  - `θ_wilt`（凋萎点）默认 0；可通过 `init_soil_param` 回填。
  - Campbell `b` 指数：USDA 13 类土壤查表（`src/SoilParameters.jl:L25-L31`）。
+ 植被参数（BIOPARMS）：
  - USDA（索引 1-12）/ IGBP（13-31）双分类，注意 `theory/04_蒸散发.typ §6` 的重复索引问题。
  - 调参前先 `clamp(veg_int, 1, 31)`。
+ 时间步参数：
  - `Δt = 3600 s`：默认小时步
  - `steps`：子步倍率（默认 1.0），CFL 受限时可调小

== 4. 示例

待补：3 个典型场景（湿润区 / 干旱区 / 冻土区）的参数表。

== 5. 常见问题

- *Q：调参后 ET 偏大/偏小？*
  A：先排查 LAI 与 BIOPARMS 索引匹配；然后检查辐射单位换算。

== 6. 参考

- `src/SoilParameters.jl §3 Campbell 公式`
- `src/Evapotranspiration.jl BIOPARMS 表`
- `wiki/julia/SoilParameters-土壤参数.md`
