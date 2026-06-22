// docs/manual/02_快速开始.typ
// ==========================
// 章：快速开始
// 状态：骨架（待摄取）

#import "../preamble.typ": *

= 快速开始

== 1. 目标

用最少的代码运行 ASAP-model 一步时间步。

== 2. 前置条件

- 已完成 `installation.typ`（即本文档前章）的安装步骤。
- 已阅读 `theory/01_模型总览.typ` 了解主循环结构。

== 3. 操作步骤

+ 单步时间步：`example/complete_example.jl::rootdepth_main(...)`
+ 多日区域应用：`example/regional_example.jl::--mock 1 --duration 168`

== 4. 示例

最小可复现示例：

```julia
using ASAP
rootdepth_main(...)
```

== 5. 常见问题

待补。

== 6. 参考

- `example/complete_example.jl`
- `example/regional_example.jl`