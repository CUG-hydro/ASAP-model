# CLAUDE.md — ASAP-model LLM 操作指南

本文件是项目级 LLM 操作守则，Claude Code / Claude API / 其他 LLM agent 在处理本项目时应优先阅读本文件。

---

## 1. 项目背景

ASAP = **A**gricultural **S**ystems **A**nalysis and **P**rediction，一个全球陆面水文模型。
项目用 Julia 重写了原始 Fortran 实现（`fortran/`），将单一 Fortran 模块拆分为多个可独立测试的 Julia 文件。

- **主入口**：`src/ASAP.jl`（`module ASAP`），通过 `include` 加载所有子模块
- **物理过程**：13 种 USDA 土壤类型、3 种 PET（P-T / P-M / S-W）、Richards 1D 土壤水、8 方向侧向流、运动波/扩散波河流路由、18O 同位素追踪（**当前阶段不串联**，见 §6.1）
- **运行方式**：`julia --project -e 'using Pkg; Pkg.test()'` 触发完整测试套件
- **完整示例**：`example/complete_example.jl` 调用 `rootdepth_main(...)` 一次时间步

---

## 2. 仓库结构（速查）

```
ASAP-model/
├── src/
│   ├── ASAP.jl              # 主入口 module
│   ├── SoilParameters.jl    # 13 种土壤类型 + Campbell 导水率
│   ├── SoilInitialization.jl# CLM/固定厚度两种土壤分层
│   ├── Evapotranspiration.jl# 三种 PET 方法（P-T/P-M/S-W）
│   ├── Interception.jl      # 植被截留（最小降水率阈值）
│   ├── extraction.jl        # 根系吸水（Feddes-style 水分胁迫）
│   ├── SoilFluxes.jl        # Richards 三对角求解 + tridag!
│   ├── updatewtd_shallow.jl # 浅层地下水位更新
│   ├── RootDepth.jl         # 主算法：rootdepth_main(...)
│   ├── helper.jl            # find_jwt, flowdir
│   ├── modules/             # 二级模块（水位/河流子系统）
│   └── backup/              # [已废弃] 禁止修改与启用
├── test/                    # 单元测试（runtests.jl 链式 include）
│   └── wtable/              # 水位/河流子系统测试
├── fortran/                 # 原始 Fortran 实现（参考用）
├── docs/                    # Typst 文档（数学公式推导）
├── example/                 # 可运行示例
├── wiki/                    # LLM 维护的中文 Wiki（详见 §3）
└── Project.toml             # Julia 依赖
```

---

## 3. Wiki 结构与约定

`wiki/` 是本项目的 LLM 维护知识库，覆盖 17 个 Julia 文件 + 11 个 Fortran 文件 + 测试 + 文档。

### 3.1 目录树

```
wiki/
├── README.md            # Wiki 首页（项目背景、目录树、快速导航）
├── index.md             # 页面索引（按类别组织）
├── log.md               # 时间线（每次 ingest/update/lint 追加条目）
├── conventions.md       # 命名约定、单位、文档模板
├── _meta/
│   ├── status.md        # 摄取状态表（页面 + 已知源码问题）
│   └── cross-refs.md    # 跨页面引用登记表
├── julia/               # 16 个 Julia 模块页面（1 主入口 + 1 主算法 + 14 子模块）
├── fortran/             # 8 个 Fortran 模块页面（1 README + 7 个 module）
└── mapping/
    ├── julia-fortran-对照.md   # 双向对照表
    └── algorithm-索引.md        # 物理过程 → 实现函数
```

### 3.2 页面模板

每个 `julia/`、`fortran/` 下的页面遵循以下 7 段模板（中文）：

```markdown
# {中文标题}

> 源文件：`src/{file}.jl:{行号段}`
> Fortran 来源：`fortran/{module}.f90:{子程序}`（如适用）
> 测试：`test/{file}.jl`
> 状态：已摄取 / 待复核

## 1. 功能概述
（2-3 句中文）

## 2. 函数签名
julia
function name(args::T) :: ReturnType
julia

## 3. 算法 / 公式
（用 LaTeX 写出关键方程；如来自 docs/*.typ 则注明）

## 4. 关键变量与单位
| 符号 | 含义 | 单位 |
| ---- | ---- | ---- |

## 5. 与 Fortran 对应
| Fortran 子程序 | Julia 函数 | 差异 |
| -------------- | ---------- | ---- |

## 6. 引用
- 行号：`src/{file}.jl:L10-L25` …
- 测试断言：`test/{file}.jl:L40` 的 @test ...
- 文档：`docs/{file}.typ`

## 7. 已知问题与备注
（如有悬空 export、废弃参数、未使用 import 等）
```

---

## 4. 核心操作

### 4.1 读 wiki（Query）

回答用户问题时：

1. 先读 `wiki/index.md` 定位相关页面
2. 读取 `wiki/julia/<模块>.md` 获取实现细节
3. 通过 `wiki/mapping/julia-fortran-对照.md` 反查 Fortran 原版
4. 通过 `wiki/mapping/algorithm-索引.md` 反查物理过程
5. 引用代码时给出 `file:line` 定位

### 4.2 更新 wiki（Ingest 源码变更）

当 `src/`、`fortran/`、`test/`、`docs/` 发生变更时：

1. 用 `sqz_read_file` 读取变更后的源文件
2. 用页面模板（§3.2）撰写更新
3. 更新 `wiki/_meta/status.md` 中对应行
4. 追加 `wiki/log.md` 条目（格式：`## [YYYY-MM-DD] update | 文件:修改摘要`）
5. 检查 `wiki/_meta/cross-refs.md` 中相关引用是否需要更新
6. 运行 lint（§4.3）

### 4.3 Lint（健康检查）

周期性健康检查：

1. 用 ripgrep 检查所有 markdown 内链目标文件是否存在
2. 检查孤立页面（无入链）
3. 对照最新源代码核对事实陈述
4. 在 `wiki/log.md` 追加 `[YYYY-MM-DD] lint | 摘要` 条目
5. 修复发现的断裂链接或孤立页面

### 4.4 关键工具

- **sqz_read_file / sqz_list_dir**：读取源码与 wiki（首选，节省 token）
- **sqz_grep**：在 wiki / 源码中搜索（首选）
- **ripgrep (rg)**：搜索文件内容（用于 lint）
- **Project.toml**：依赖管理
- **test/runtests.jl**：测试入口（`julia --project -e 'using Pkg; Pkg.test()'`）
- **rtk**：命令代理（git/rg/test 等）

---

## 5. 源码组织约定

- `src/*.jl`：顶层模块（按物理过程划分，无嵌套 module）
- `src/modules/*.jl`：二级模块（水位/河流子系统）
- `src/modules/Tracing/*.jl`：三级模块（同位素追踪）
- `src/backup/`：**已废弃**，禁止修改与启用
- `test/*.jl` 与 `test/wtable/*.jl`：两级测试套件
- `fortran/`：原始 Fortran 实现（参考用，对应关系见 `wiki/mapping/julia-fortran-对照.md`）

### 5.1 命名约定

| 符号                        | 含义         |
| --------------------------- | ------------ |
| `θ` / `theta`               | 体积含水量   |
| `ψ` / `psi`                 | 基质势（m）  |
| `κ` / `kappa`               | 导水率       |
| `ρ` / `rho`                 | 密度相关参数 |
| `Δt`                        | 时间步长     |
| `α` / `β` / `γ` / `δ` / `λ` | 其他参数     |

### 5.2 单位约定

- 长度：m
- 时间：s
- 流量：m³/s（转 mm 需 ×1000）
- 质量：kg
- 温度：K（热力学）
- 数组索引：遵循 Julia 1-based 约定

---

## 6. 禁止事项

- ❌ **启用 `src/backup/` 中的任何代码**（已废弃，被新模块替代）
- ❌ **删除或修改** `wiki/_meta/status.md` 中已记录的已知源码问题清单（待后续 PR 修复）
- ❌ **在 wiki 中复制超过 30 行的源代码**（保持 wiki 是「地图」而非源码副本）
- ❌ **修改 `fortran/` 目录**（参考用）
- ❌ **修改 `docs/*.typ` 文档**（学术写作风格规范）
- ❌ **不修改 `Project.toml` 依赖版本约束**（`compat` 段）
- ❌ 不删除有意义的注释或有意义的注释掉的代码（等未来可恢复的代码，或有意义的注释）

### 6.1 现阶段不实现的范围（2026-06-22 起生效）

用户明确以下三类功能**当前阶段均不需要实现**，禁止在本阶段任务中引入：

- ❌ **NetCDF I/O 扩展**：仅 `src/io/NetCDF.jl::NetCDFIO` 已落地的 P0 子集（`read_initial` / `read_wtdnc`）是当前终点。**禁止**实现写出（`WRITEOUTPUTNC*` / `WRITEHISTORYNC*`）、其它读取子程序（`READLATLON` / `READVEG` / `READHVEG` / `READSMOIEQ` / `READFLOWDIRECTION` / `READRIVERPARAMETERS` / `READHISTORYNC*`）或新增 NetCDF 子模块。
- ❌ **MPI 并行化**：Julia 端是单进程设计。**禁止**引入 `MPI.jl` 或任何并行域分解 / halo 通信代码（`fortran/module_parallel.f90` 保持只读）。
- ❌ **¹⁸O 同位素串联**：`src/modules/Tracing/IsotopeTracing.jl` 中已实现的 `lateral_isotope!` / `lateralflow_with_isotope!` / `updatedeepwtable!` / `updatewtd_simple` 维持现状（保留单测）。**禁止**启用 `src/SoilFluxes.jl` 末段被注释的 ¹⁸O 三对角组装、**禁止**在 `src/RootDepth.jl` 中调用 `lateral_isotope!` 或 `updatedeepwtable!`、**禁止**实现 `fortran/module_rootdepth.f90` 的 Majoube 平衡分馏 α 公式。

后续如需启用上述任一功能，必须由用户明确重新授权，并在 `wiki/log.md` 追加 `[YYYY-MM-DD] enable | 范围` 条目。

---

## 7. 已知源码问题（截至 2026-06-22）

10 个问题已识别并修复（参见 commit 摘要）；剩余清理项见 `wiki/_meta/status.md` §5：

| #   | 问题                                             | 状态       | 修复位置                                    |
| --- | ------------------------------------------------ | ---------- | ------------------------------------------- |
| 1   | `SoilParameters.init_soil_param.fieldcp` 占位 0  | ✅ 已修复   | `src/SoilParameters.jl:73-L78`              |
| 2   | `SoilFluxes.jl` 氧 18 注释段 + docstring 不一致  | ✅ 已清理   | `src/SoilFluxes.jl`、`src/RootDepth.jl`     |
| 3   | `Modules.jl` 悬空 export `wtable!`/`updatewtd!`  | ✅ 已删除   | `src/modules/Modules.jl`                    |
| 4   | `ASAP.jl` 悬空 export `updatewtd_qlat`           | ✅ 已删除   | `src/ASAP.jl`                               |
| 5   | `extraction.kroot = k - 1` 注释误导              | ✅ 已澄清   | `src/extraction.jl`                         |
| 6   | `updatewtd_shallow.flag` 未初始化隐患            | ✅ 注释加固 | `src/updatewtd_shallow.jl`                  |
| 7   | `Interception.minpprate` docstring 错误          | ✅ 已修正   | `src/Interception.jl`                       |
| 8   | `potevap_shutteworth_wallace` 拼写错误           | ✅ 新增别名 | `src/Evapotranspiration.jl`                 |
| 9   | `SoilInitialization.DataFrames` 悬空 import      | ✅ 已删除   | `src/SoilInitialization.jl`、`Project.toml` |
| 10  | `rivers_*` 中形参 `length` 与 `Base.length` 同名 | ✅ 已重命名 | `rivers_kw_flood.jl`、`rivers_dw_flood.jl`  |

---

## 8. 测试入口

```bash
# 完整测试套件
julia --project -e 'using Pkg; Pkg.test()'

# 单文件测试
julia --project test/test_soil_parameters.jl

# 完整示例
julia --project example/complete_example.jl
```

---

## 9. 跨语言映射速查

| Julia 函数                    | Fortran 子程序                | 差异                                               |
| ----------------------------- | ----------------------------- | -------------------------------------------------- |
| `rootdepth_main`              | `ROOTDEPTH`                   | 接口对齐；新增泛型签名                             |
| `potevap_shutteworth_wallace` | `POTEVAP_Shutteworth_Wallace` | 拼写相同；Julia 返回元组                           |
| `soilfluxes`                  | `SOILFLUXES`                  | 同位素段 Julia 已禁用（**现阶段不串联**，见 §6.1） |
| `extraction`                  | `EXTRACTION`                  | 1-based 一致                                       |
| `interception`                | `INTERCEPTION`                | 接口一致                                           |
| `lateral_flow!`               | `LATERALFLOW4`                | 8 邻居（D8）                                       |
| `gw2river!`                   | `GW2RIVER`                    | 三种河床交换                                       |
| `rivers_kw_flood!`            | `RIVERS_KW_FLOOD`             | 形参 `length`→`river_length`                       |
| `rivers_dw_flood!`            | `RIVERS_DW_FLOOD`             | 形参 `length`→`river_length`                       |
| `flooding!`                   | `FLOODING`                    | 8 邻居漫流                                         |
| `updatewtd_shallow`           | `UPDATESHALLOWWTD`            | Julia 简化（去除 qlat）                            |
| `find_jwt`                    | 无对应                        | Julia 工具                                         |

完整对照表见 `wiki/mapping/julia-fortran-对照.md`。

---

## 10. 收尾约定

- 完成任何实质性修改后，**运行** `julia --project -e 'using Pkg; Pkg.test()'`
- 完成 wiki 变更后，**追加 log.md 条目**
- 修改源码后，**更新对应 wiki 页面 + status.md**
- 检测到新问题时，**追加到 `_meta/status.md` §5**

如有疑问，先读 `wiki/README.md` 与 `wiki/conventions.md`。
