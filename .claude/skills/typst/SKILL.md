---
name: typst
description: >
  Typst 0.15.0 学术文档编写规则与避坑指南。当用户编写或调试 .typ 文件、
  遇到编译报错（unknown variable、cannot reference equation、path would escape）、
  或需要使用以下包时触发：physica:0.9.8（偏导数/梯度/散度）、
  modern-cug-report:0.1.3（学术报告模板）、tablem:0.3.0（Markdown 风格表格）。
  关键词：#import、#show、#set、math equation、tablem、physica、pdv、
  underbrace、@label、context、#include、--root 编译参数。
---

**Typst 编写技能指南**

版本：Typst 0.15.0。适用于使用 `physica:0.9.8` + `modern-cug-report:0.1.3` + `tablem:0.3.0` 的学术文档。

---

## 1 编译工作流

```bash
# 单文件（若 #import 跨越了目录，必须加 --root）
typst c --root . docs/manual/03_输入数据.typ /tmp/out.pdf

# 主文档（推荐：所有子文件 #include 进主文档）
typst c --root . docs/main.typ /tmp/out.pdf 2>&1 | grep -E '^error' -A6

# 检查 warning（字体类可忽略，语法类必修）
typst c --root . docs/main.typ /tmp/out.pdf 2>&1 | grep 'warning:' | sort | uniq -c
```

`--root .` 原因：`#import "../preamble.typ": *` 会尝试访问项目根目录之外的路径，不加 `--root` 时报 `path would escape the project root`。

---

## 2 数学模式：多字母标识符必须加引号

**这是最常见的 error，不报 warning，改了才知道。**

Typst 数学模式中，连续字母串被解析为"变量名"。`ssrd`、`tp`、`qair` 等非单字母符号在 math 模式下触发 `unknown variable: ssrd`。

```typ
// ❌ 报错：unknown variable: ssrd / tp / netrad
$ssrd (1 - alpha) - strd$
$tp times 1000$

// ✅ 变量加引号 → 渲染为直立文本
$"ssrd" (1 - alpha) - "strd"$
$"tp" times 1000$

// ✅ 或用 upright()（physica 风格）
$upright("ssrd") (1 - alpha) - upright("strd")$
```

自查命令（找出 math 模式内两字母以上的未加引号标识符）：

```bash
rg -n '\$[^$]*[a-z]{3,}[^$"]*\$' docs/
```

---

## 3 数学下标/上标：用 `(…)` 不是 LaTeX 的 `{…}`

花括号 `{}` 在 Typst math 里渲染为**字面量**（不是分组），不报错但结果错。

```typ
// ❌ LaTeX 写法，Typst 把 {…} 当字面量显示
theta_{k-1}^{n+1}
Q_{k+1/2}

// ✅ Typst 写法
theta_(k-1)^(n+1)
Q_(k+1/2)
```

另一陷阱：下标无括号时只吃紧邻的一个 token。`Q_k+1/2` 只把 `k` 当下标，
渲染成 `Q_k + 1/2` → 必须写 `Q_(k+1/2)`。

自查：`rg -n '[{}]' docs/ --include='*.typ'`（排除 `#if`、代码块中的花括号）。

---

## 4 粗体/斜体：与 Markdown 相反

| 效果     | Markdown | Typst |
| -------- | -------- | ----- |
| **粗体** | `**x**`  | `*x*` |
| _斜体_   | `*x*`    | `_x_` |

`**x**` 在 Typst 里被解析为"空强调 + 文本"，warning `no text within stars`，且不加粗。

> 缺字体只是 warning，可忽略：`#set text(font: (...))` 是回退列表，缺哪个跳过哪个，不阻断编译。

---

## 5 physica:0.9.8 常用命令

通过 `#import "@preview/physica:0.9.8": *` 引入（preamble.typ 已全局导入）。

| 命令           | 输出    | 用途     |
| -------------- | ------- | -------- |
| `pdv(f, x)`    | ∂f/∂x   | 偏导数   |
| `pdv(f, x, 2)` | ∂²f/∂x² | 二阶偏导 |
| `dv(f, x)`     | df/dx   | 全导数   |
| `div A`        | ∇·A     | 散度     |
| `grad f`       | ∇f      | 梯度     |
| `laplacian f`  | ∇²f     | 拉普拉斯 |
| `dd(x)`        | dx      | 微分符号 |
| `abs(x)`       | \|x\|   | 绝对值   |
| `norm(x)`      | ‖x‖     | 范数     |

```typ
// Richards 方程示例
$ pdv(theta, t) = pdv(, z)[D(theta) pdv(theta, z) - K(theta)] - S $

// 通量
$ Q_(k+1/2) = -K_(k+1/2) - D_(k+1/2) (theta_(k+1)^(n+1) - theta_k^(n+1)) / (Delta z_(k+1/2)) $
```

---

## 6 modern-cug-report:0.1.3 模板使用

用于独立编译的学术章节（讲稿 / 单章报告）。若文件被 `#include` 进 main.typ，则 **不要** 加 `#show` rule（否则每个 include 块都会独立套一层 header/footer/页码）。

```typ
// 独立编译的文件（如 theory/*.typ 单独输出 PDF）
#import "@preview/modern-cug-report:0.1.3": *
#show: doc => template(
  doc,
  size: 12pt,          // 正文字号
  footer: "CUG水文气象学2026",   // 页脚文字
  header: "",          // 页眉文字（空字符串 = 无页眉）
)
#import "../preamble.typ": *   // 必须在 #show 之后导入，否则 tablem 等宏不可用
```

**顺序关键**：`#show: doc => template(...)` 必须在 `#import "../preamble.typ": *` 之前，否则 preamble 中的 `#set` 规则优先级可能被覆盖。

---

## 7 tablem:0.3.0 Markdown 表格

Typst 原生不渲染 `| a | b |` 风格表格，须包进 `tablem`（preamble.typ 已全局导入）。

```typ
// 普通有框表格
#figure(
  caption: [],
  align(center)[
  #tablem()[
    | *列头 A* | *列头 B* | *单位* |
    | -------- | -------- | ------ |
    | `STXT`   | 土壤质地 | —      |
    | `topo`   | 高程     | m      |
  ]
])

// 学术三线表（无竖线，顶线/底线/中线）
#figure(
  caption: [],
  align(center)[
  #three-line-table()[
    | 符号       | 含义       | 单位  |
    | ---------- | ---------- | ----- |
    | $theta$    | 体积含水量 | m³/m³ |
    | $K(theta)$ | 导水率     | m/s   |
  ]
])
```

列对齐：`| --- |` 左对齐（默认），`| :---: |` 居中，`| ---: |` 右对齐。

> tinymist IDE 对经 preamble 再导出的包符号会误报 `unknown variable`；以 `typst c` CLI 结果为准。

---

## 8 公式编号与引用

```typ
// 文件顶部声明（必须在引用任何 @label 之前）
#set math.equation(numbering: "(1)", supplement: none)
// supplement: none → 引用显示「式 (3)」而非「式 Equation 3」

// 给公式打标签
$ Q_(k+1/2) = ... $ <eq_flux>

// 引用
将 @eq_flux 代入连续方程 ...
```

若未设 `numbering` 就用 `@label`，报 `cannot reference equation without numbering`。

---

## 9 `context` 只覆盖紧跟其后的表达式

`counter(page).final()`、`measure()` 等"间接值"只能在 `context` 块内求值。

```typ
// ❌ context 只罩右半，左边的 display() 在 context 外报错
[#counter(page).display() / #context counter(page).final().first()]

// ✅ 整段包进一个 context
[#context [#counter(page).display() / #counter(page).final().first()]]
```

---

## 10 静默误渲染（不报错，须主动 grep）

| 误写（Markdown） | Typst 正确写法       | 自查                        |
| ---------------- | -------------------- | --------------------------- |
| `## 标题`        | `== 标题`            | `rg -n '^\s*#{2,}\s' docs/` |
| `[文字](url)`    | `#link("url")[文字]` | `rg -n '\]\(' docs/`        |
| `> 引用块`       | 自定义 block         | `rg -n '^\s*> ' docs/`      |
| `theta_{k}`      | `theta_(k)`          | `rg -n '[{}]' docs/`        |
| `**粗体**`       | `*粗体*`             | `rg -n '\*\*' docs/`        |

一键扫描：

```bash
rg -n '^\s*#{2,}\s|\]\(|^\s*> |\*\*|[{}]' docs/ --include='*.typ'
```

---

## 11 #include 中的 #show 规则作用域

Typst 中 `#show` 规则的作用域是**当前 content block**。在 `#include "chapter.typ"` 中定义的 show rule 只影响该文件内部，不泄漏到主文档。

因此，standalone 编译用的 `#show: doc => template(...)` 放在被 include 的文件里也是安全的——但仍建议只在真正独立编译的文件中加，include 场景统一由 main.typ 控制模板。

---

## 12 分数与常见数学写法速查

| LaTeX                | Typst                             |
| -------------------- | --------------------------------- |
| `\frac{a}{b}`        | `a/b`                             |
| `\partial`           | `diff`                            |
| `\Delta`             | `Delta`                           |
| `\nabla`             | `nabla`                           |
| `\cdot`              | `dot`                             |
| `\times`             | `times`                           |
| `\leq`               | `lt.eq`                           |
| `\geq`               | `gt.eq`                           |
| `\infty`             | `infinity`                        |
| `\text{abc}`         | `"abc"`                           |
| `\mathrm{d}`         | `dd(x)`（physica）或 `upright(d)` |
| `\underbrace{x}_{y}` | `underbrace(x, y)`                |
| `\left( \right)`     | 自动大括号，无需写                |
| `\begin{cases}`      | `cases(...)`                      |

```typ
// 联立方程
$
cases(
  a + b = 1,
  c - d = 2,
)
$
```
