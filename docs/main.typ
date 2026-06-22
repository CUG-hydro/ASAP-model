// docs/main.typ — ASAP-model 文档编排入口
// ============================================
// 编译命令：
//   typst compile docs/main.typ
//   typst watch   docs/main.typ    # 开发期实时预览
//
// 输出：ASAP-model 原理与用户手册.pdf

#import "preamble.typ": *

#show: doc => template(doc, footer: "ASAP-model 原理与用户手册", header: "")

// ===================== 封面 =====================
#titlepage(
  title: "ASAP-model 原理与用户手册",
  subtitle: [
    Agricultural Systems Analysis and Prediction \
    一个全球陆面水文模型的 Julia 实现
  ],
  author: "ASAP-model 开发组",
  date: datetime(year: 2026, month: 6, day: 22),
)

// ===================== 目录 =====================
#outline(depth: 2, indent: auto)

// ===================== 摘要 =====================
= 摘要

ASAP-model 是一个全球陆面水文模型，原始实现基于 Fortran，本项目用 Julia 重写，
将单一 Fortran 模块拆分为 17 个可独立测试的源文件。

本文档分两部分：

- *原理篇*（`theory/`）：覆盖模型每个物理过程的数学推导、数值方法、参数说明，
  与源代码 `src/*.jl` 行级映射。
- *手册篇*（`manual/`）：覆盖安装、快速开始、输入数据准备、输出诊断、调参指南，
  与 `example/*.jl` 实例对应。

// ===================== 原理篇 =====================
#counter(heading).update(0)
#part("原理篇")
#include "theory/01_模型总览.typ"
#include "theory/02_土壤水运动.typ"
#include "theory/03_根系吸水.typ"
#include "theory/04_蒸散发.typ"
#include "theory/05_侧向流.typ"
#include "theory/06_河流路由.typ"
#include "theory/07_同位素.typ"
#include "theory/08_数值方法.typ"

// ===================== 手册篇 =====================
#counter(heading).update(0)
#part("手册篇")
#include "manual/01_安装.typ"
#include "manual/02_快速开始.typ"
#include "manual/03_输入数据.typ"
#include "manual/04_输出诊断.typ"
#include "manual/05_调参指南.typ"

// ===================== 附录 =====================
#counter(heading).update(0)
#part("附录")

== 源码 ↔ 文档对照表

参见 `wiki/mapping/julia-fortran-对照.md`（Fortran ↔ Julia 行级映射）与
`wiki/mapping/algorithm-索引.md`（18 类物理过程反向索引）。

== 参考文献

#bibliography("references.bib", title: none)
