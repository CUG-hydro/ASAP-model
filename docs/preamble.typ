// docs/preamble.typ — ASAP-model 文档全局设置
// =====================================================
// 此文件由 docs/main.typ include；现有讲稿式 .typ 文件（土壤水运动_ASAP.typ /
// 汇流_扩散波.typ / 下渗_黄土高原.typ）通过以下命令单独编译时也需引入：
//
//   typst compile --input preamble=docs/preamble.typ <file>.typ
//
// 或在该文件首行加：#include "preamble.typ"

// ---------- 字体（CJK 优先于西文） ----------
// 使用 Typst 内置字体 + 系统 CJK 字体回退；缺哪个跳过哪个（不阻断编译）。
#set document(
  title: "ASAP-model 原理与用户手册",
  author: "ASAP-model 开发组",
)

#set text(
  // font: (
  //   "New Computer Modern",
  //   "Noto Serif CJK SC",
  //   "Source Han Serif SC",
  //   "Noto Sans CJK SC",
  //   "AR PL UKai CN",
  // ),
  size: 11pt,
  lang: "zh",
  region: "CN",
)

// ---------- 数学 ----------
// Typst 原生数学模式；math 宏来自 physica 包（提供 pdv / tensor / bra-ket 等）。
// https://typst.app/universe/package/physica/
#import "@preview/physica:0.9.8": *

// Markdown 风格表格：Typst 原生不渲染 `| a | b |`，借助 tablem 包处理。
// 用法：#tablem[\n| 表头 | … |\n| --- | --- |\n| 单元 | … |\n]
// 该 import 经 `#import "preamble.typ": *` 被各 theory/manual 文件继承。
#import "@preview/tablem:0.3.0": tablem, three-line-table

// 现有讲稿式文件（土壤水运动_ASAP.typ / 汇流_扩散波.typ）使用 mitex 兼容 LaTeX：
//   #import "@preview/mitex:0.2.4": *
//   #mitex(`$$ ... $$`)
// 两套风格可在不同文件中并存，本文档（theory/ + manual/）统一用 physica。

// ---------- heading 样式 ----------
#set heading(numbering: "1.1")
#show heading.where(level: 1): set text(size: 16pt, weight: "bold")
#show heading.where(level: 2): set text(size: 13pt, weight: "bold")
#show heading.where(level: 3): set text(size: 11.5pt, weight: "bold")

// ---------- 行内代码 / 代码块 ----------
#show raw.where(block: false): box.with(
  fill: luma(240),
  inset: (x: 3pt, y: 0pt),
  outset: (y: 2pt),
  radius: 2pt,
)
#show raw.where(block: true): block.with(
  fill: luma(248),
  inset: 8pt,
  radius: 3pt,
  width: 100%,
)

// ---------- 引用颜色 ----------
// 已有讲稿式文件使用 `markhl(color:)` 高亮变量；以下为极简实现。
// 完整版可改用 `#import "@preview/annotate:0.1.0": *`。
#let markhl(color: yellow, body) = {
  box(
    fill: color.lighten(60%),
    outset: (y: 2pt),
    radius: 1pt,
    body,
  )
}
#let markhl-red = markhl.with(color: red)
#let markhl-yellow = markhl.with(color: yellow)
#let markhl-blue = markhl.with(color: blue)

// ---------- 警示框 ----------
#let box-red(body) = block(
  fill: red.lighten(85%),
  stroke: 1pt + red,
  inset: 8pt,
  radius: 3pt,
  body,
)
#let box-note(body) = block(
  fill: blue.lighten(90%),
  stroke: 1pt + blue,
  inset: 8pt,
  radius: 3pt,
  body,
)

// ---------- template：讲稿页眉页脚 ----------
// 兼容现有 .typ 文件的 `#show: doc => template(doc, footer: "...", header: "...")` 约定。
#let template(doc, footer: "", header: "") = {
  set page(
    paper: "a4",
    margin: (x: 2.2cm, y: 2.5cm),
    header: [
      #set text(size: 9pt, fill: gray)
      #grid(
        columns: (1fr, auto),
        align: (left, right),
        [#header], [#context counter(page).display()],
      )
      #line(length: 100%, stroke: 0.4pt + gray)
    ],
    footer: [
      #set text(size: 9pt, fill: gray)
      #line(length: 100%, stroke: 0.4pt + gray)
      #grid(
        columns: (1fr, auto, 1fr),
        align: (left, center, right),
        [#footer], [], [#context [#counter(page).display() / #counter(page).final().first()]],
      )
    ],
  )
  doc
}

// ---------- titlepage：主文档封面 ----------
// 由 main.typ 的 `#titlepage(...)` 调用。
#let titlepage(title: "", subtitle: none, author: "", date: none) = {
  set page(
    margin: (x: 3cm, y: 4cm),
    header: none,
    footer: none,
  )
  v(4cm)
  align(center)[
    #set text(size: 28pt, weight: "bold")
    #title
    #v(0.5em)
    #set text(size: 14pt, weight: "regular", fill: gray)
    #if subtitle != none { subtitle }
    #v(2em)
    #set text(size: 12pt)
    #author
    #v(1em)
    #if date != none [
      #set text(size: 11pt, fill: gray)
      #date.display("[year] 年 [month] 月 [day] 日")
    ]
  ]
  pagebreak()
}

// ---------- part：篇章分隔页（强制换页 + 居中标题） ----------
#let part(title) = {
  pagebreak()
  v(1fr)
  align(center)[
    #set text(size: 24pt, weight: "bold")
    #title
  ]
  v(1fr)
  pagebreak()
}

// ---------- 数学宏补充（physica 未覆盖的） ----------
// 数值离散符号（Typst 原生已有 Delta，以下补小写）
#let dt = math.upright("dt")
// 平衡态标记
#let eq = math.op("eq")

// ---------- 引用源码约定 ----------
// docs/ 中引用源码使用以下格式（一致于 wiki）：
//   src/SoilFluxes.jl:L42-L98
//   src/RootDepth.jl:L131-L142（步骤 ① SW PET）
