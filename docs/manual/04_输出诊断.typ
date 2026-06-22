// docs/manual/04_输出诊断.typ
// ==========================
// 章：输出诊断
// 状态：骨架（待摄取；按 CLAUDE.md §6.1 阶段不实现 NetCDF 写出）

#import "../preamble.typ": *

= 输出诊断

== 1. 目标

读懂 ASAP-model 主循环产生的状态变量与诊断字段。

== 2. 前置条件

已完成 `manual/03_输入数据.typ` 的准备并跑通至少一个时间步。

== 3. 操作步骤

+ 关键输出字段清单（详见 `theory/01_模型总览.typ §4`）：
  - `et_s` / `et_c` / `et_i`：土壤/冠层/截留蒸发累计
  - `qsrun`：地表径流累计
  - `rech` / `deeprech`：浅层/深层补给
  - `wtd`：地下水位深度
  - `waterdeficit`：水分亏缺累计
  - `infilk`：入渗深度层索引
  - `icefactor`：冰冻抑制因子
+ 单位约定：所有通量 `m³/m³ → mm` 需 ×1e3；`m³/s → mm/h` 需 ×3.6e6。

== 4. 示例

```julia
# 读出最后一帧 wtd 与累计 et_s
@show wtd[end, end]
@show sum(et_s) * 1e3  # mm
```

== 5. 常见问题

- *Q：如何判断数值守恒？*
  A：`sum(θ_new) ≈ sum(θ_old) + (precip - et - runoff) * Δt`，详见 `test/test_extraction.jl:L118`。

== 6. 参考

- `src/RootDepth.jl §4 关键变量与单位`
- 按 §6.1，*NetCDF 写出（`WRITEOUTPUTNC*` / `WRITEHISTORYNC*`）当前阶段不实现*。