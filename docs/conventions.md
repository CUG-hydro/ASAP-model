# docs/conventions.md — Typst 写作规范

> 本文件定义 ASAP-model 文档的写作约定，供新增章节、迁移讲稿式 .typ 时参考。
> 适用范围：`docs/theory/*.typ`、`docs/manual/*.typ`；现有讲稿式 .typ
> （`土壤水运动_ASAP.typ` 等）保持不动，仅在新文档中遵循本规范。

---

## 1. 命名与目录

| 层级       | 命名规则                         | 示例                                  |
| ---------- | -------------------------------- | ------------------------------------- |
| 章文件     | `NN_中文标题.typ`                | `02_Richards方程.typ`                 |
| 手册章文件 | `NN_中文标题.typ`                | `03_输入数据.typ`                     |
| 图文件     | `Figure_<模块>.typ`（cetz）      | `cetz/Figure_Soil_ASAP.typ`           |
| 引用键     | `fig:xxx`、`eq:xxx`、`tbl:xxx`   | `fig:soil_column`、`eq:richards_1d`   |
| 源文件引用 | `src/<file>.jl:L<start>-L<end>`  | `src/SoilFluxes.jl:L42-L98`           |

**禁止**：
- ❌ 在 `docs/` 根新建散落的 .typ 文件（统一归入 `theory/` 或 `manual/`）。
- ❌ 修改 `cetz/` 子目录外的现有 .typ 文件（讲稿风格，单独编译）。

## 2. 数学

使用 Typst **原生数学模式**：

```typst
$
  pdv(theta, t) = pdv(, z)[D(theta) pdv(theta, z) - K(theta)] - S
$
```

**避免**：
- ❌ LaTeX 行内命令（`\frac{}{}`），除非要兼容 LaTeX 输出。
- ✅ 优先用 `#import "@preview/mitex:0.2.4": *` + `#mitex("$$...$$")` 复用 LaTeX 笔记。

常用算子（`preamble.typ` 已定义）：

| 宏       | 含义     | 渲染 |
| -------- | -------- | ---- |
| `pdv`    | 偏导算子 | ∂    |
| `Delta`  | 增量     | Δ    |
| `dt`     | 时间步长 | dt   |

## 3. 代码块

```typst
```julia
# 注释说明
function rootdepth_main(...) :: Nothing
  ...
end
```
```

约定：
- Julia 代码块默认 `lang: "julia"`。
- 行号引用跟在代码块下方，不内嵌。
- 公式与代码的对应关系用 `#figure()` 包裹居中。

## 4. 公式引用

```typst
$ Q = K + D frac(pdv(theta), pdv(z)) $ <eq:richards_flux>

据式 #[@eq:richards_flux] 可得三对角系数：
```

- 公式标签：`#[@eq:key]`
- 图标签：`#[@fig:key]`
- 表标签：`#[@tbl:key]`
- 文献标签：`#[@bibkey]`

## 5. 章节模板（理论篇）

每章遵循以下 6 段结构（与 wiki 一致但视角不同）：

```typst
= 章标题

== 1. 物理过程
（一段话讲清楚本章建模对象与物理图景；2-3 句中文。）

== 2. 控制方程
$ ... $ <eq:key>
$ ... $

== 3. 数值离散与求解
（针对 1D/2D 方程组的离散格式、求解算法、稳定性条件。）

== 4. 参数与单位
#figure(
  table(
    columns: (auto, auto, 1fr, auto),
    align: (left, left, left, left),
    [符号], [含义], [说明], [单位],
    ...,
  ),
  caption: [本章关键变量],
)

== 5. 代码实现
- 关键函数：#raw[src/Module.jl:Lstart-Lend]（功能）
- 测试断言：#raw[test/test_xxx.jl:Lx]
- 参考：#link("../wiki/julia/Module-xxx.md")[Module wiki 页面]

== 6. 局限与已知问题
（数值方法的限制、参数敏感性、与 Fortran 实现的差异等。）
```

## 6. 章节模板（手册篇）

```typst
= 章标题

== 1. 目标
（读完本章能做什么）

== 2. 前置条件
（依赖、环境、权限）

== 3. 操作步骤
+ 第一步
+ 第二步

== 4. 示例
（最小可复现的命令序列与预期输出）

== 5. 常见问题
（FAQ；与第 3 节互引）

== 6. 参考
- 源码：#raw[src/Module.jl]
- 脚本：#raw[example/xxx.jl]
```

## 7. 与 wiki/ 的边界

| 维度     | wiki/                          | docs/                          |
| -------- | ------------------------------ | ------------------------------ |
| 受众     | 模型维护者                     | 模型开发者 + 用户              |
| 视角     | 代码视角（行级、函数级）       | 理论视角（章节、公式、推导）   |
| 语言     | 中文                           | 中文 + 英文术语                |
| 数学     | 极简（突出关键方程）           | 完整推导                       |
| 文献引用 | 无                             | 有（references.bib）           |
| 同步     | 源码改 → 改 wiki              | 公式改 → 改 docs + 检查 wiki  |

**禁止**：wiki 中超过 30 行的源代码（CLAUDE.md §6）；docs 中重复 wiki 的逐行描述。

## 8. 编译与预览

```bash
# 全文档编译
typst compile docs/main.typ

# 开发期实时预览（自动重编译）
typst watch docs/main.typ

# 单独编译讲稿式文件（需提供 template 函数）
typst compile --input preamble=docs/preamble.typ docs/土壤水运动_ASAP.typ
```

## 9. 修改纪律

- ❌ 不要修改 `preamble.typ` 中的 `template` 函数签名（现有讲稿式 .typ 依赖它）。
- ❌ 不要删除 `references.bib` 中已声明的引用键（可能导致 cross-ref 失效）。
- ✅ 修改章节后运行 `typst compile docs/main.typ` 验证编译通过。
- ✅ 改动涉及 src/ 时，**双向同步**：先改 docs，再回头校验 wiki。