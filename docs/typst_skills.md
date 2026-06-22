# Typst 排错经验（ASAP-model docs）

Typst 0.15.0。编译：`typst compile docs/main.typ`。
theory/、manual/ 各文件靠 `#include` 进 main.typ，依赖其全局 `#set` 与 preamble 宏，**不单独编译**。

## 1. `context` 只覆盖紧跟其后的表达式

`counter(page).display()`、`.final()`、`measure()`、`locate()` 等"间接值"只能在 `context` 内求值。
陷阱：`#context` 不向左覆盖同行前面的内容。

```typ
// ❌ context 只罩住右半
[#counter(page).display() / #context counter(page).final().first()]
// ✅ 整段包进一个 context
[#context [#counter(page).display() / #counter(page).final().first()]]
```

## 2. 粗体是单星 `*…*`，不是 Markdown 的 `**…**`

`**x**` 被解析成"空强调+文本" → warning `no text within stars`，且不加粗。

| | Markdown | Typst |
|---|---|---|
| 粗体 | `**b**` | `*b*` |
| 斜体 | `*i*` | `_i_` |

## 3. 缺字体只是 warning，可忽略

`#set text(font: (...))` 是回退列表，缺哪个跳过哪个，不阻断编译。

## 4. 静默误渲染的 Markdown 残留（不报错也不报警，须主动 grep）

| Markdown | Typst 正确写法 |
|---|---|
| `## 标题` | `== 标题` |
| `[文字](url)` | `#link("url")[文字]` |
| `> 引用` | 自定义 block |

```bash
rg -n '^\s*#{2,}\s|\]\(|^\s*> |\*\*' docs/
```

## 5. 数学：下标分组用 `(…)`，不是 LaTeX 的 `{…}`

`theta_{k-1}^{n+1}` 在 Typst 里把花括号**当字面量渲染**（显示成 `θ^{n+1}_{k-1}`），不报错。

```typ
// ❌ theta_{k-1}^{n+1}   Q_{k+1/2}   s^{-1}
// ✅ theta_(k-1)^(n+1)   Q_(k+1/2)   s^(-1)
```

另：`Q_k+1/2` 只把 `k` 当下标 → 写 `Q_(k+1/2)`。自查 `rg -n '[{}]' file.typ`（排除 `#if {}` 等代码块）。

## 6. 公式引用：先开编号，再去 supplement

`@label` 引用公式前必须 `#set math.equation(numbering: "(1)")`，否则报
`cannot reference equation without numbering`。中文文档加 `supplement: none`，
引用才显示「式 (3)」而非「式 Equation 3」。

```typ
#set math.equation(numbering: "(1)", supplement: none)
$ ... $ <eq_2>
... 将式 @eq_2 带入 ...
```

## 7. Markdown 表格：Typst 原生不渲染，须用 tablem 包

`| a | b |` 会被当普通文本（一行带竖线），须包进 tablem：

```typ
#import "@preview/tablem:0.3.0": tablem, three-line-table
#three-line-table[          // 学术三线表；普通边框用 #tablem[…]
| 符号 | 含义 | 单位 |
| --- | --- | --- |
| $theta$ | 体积含水量 | m³/m³ |
]
```

> tinymist(IDE) 对 tablem 等**经 preamble 再导出**的包符号会误报
> `unknown variable` / `unclosed delimiter`；以 `typst compile` CLI 结果为准。

## 8. 工作流

```bash
typst compile docs/main.typ /tmp/o.pdf 2>&1 | grep -E '^error' -A6   # 先修 error，逐个重编
typst compile docs/main.typ /tmp/o.pdf 2>&1 | grep warning: | sort | uniq -c   # 字体类可忽略，语法类必修
```
再补查第 4 节静默项。
