// docs/manual/01_安装.typ
// ======================
// 章：安装
// 状态：骨架（待摄取）

#import "../preamble.typ": *

= 安装

== 1. 目标

完成 ASAP-model 的本地安装与依赖配置，能在 REPL 中执行 `using ASAP`。

== 2. 前置条件

- Julia 1.6+（推荐 1.10 LTS）
- NCDatasets.jl 0.12-0.14（NetCDF P0 子集）
- 操作系统：Linux / macOS / Windows（理论支持，未完整测试 Windows）

== 3. 操作步骤

+ 克隆仓库
  ```bash
  git clone https://github.com/CUG-hydro/ASAP-model.git
  cd ASAP-model
  ```

+ 启动 Julia 项目环境
  ```bash
  julia --project
  ```

+ 安装依赖
  ```julia
  using Pkg
  Pkg.instantiate()
  ```

+ 验证安装
  ```julia
  using ASAP
  ```

== 4. 示例

== 5. 常见问题

- *Q：报错 `NCDatasets not found`？*
  A：运行 `Pkg.add("NCDatasets")`。

== 6. 参考

- 源码：`src/ASAP.jl`
- 依赖：`Project.toml`