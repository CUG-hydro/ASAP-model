# 命名、单位与文档模板

> 全 Wiki 共用的命名约定、单位约定与页面模板。新建页面请遵循此规范。

## 1. 命名约定

### 1.1 物理量符号（Unicode 优先）

| 符号 | 含义 | 备注 |
|---|---|---|
| `θ` | 体积含水量 (m³/m³) | `theta`；范围 `[θ_wilt, θ_sat]` |
| `ψ` | 基质势 (m) | `psi`；饱和时为负；POTLEAF = POTWILT = −153.0 m |
| `κ` | 导水率 (m/s 或 m²/s) | `kappa`；`Ksat` 为饱和值，`κlat` 为侧向因子 |
| `ρ` | 密度 (kg/m³) | `rho`；仅在 `gw2river` 等模块间接使用 |
| `Δt` | 时间步长 (s) | `Δt` 写作 Unicode；常量子模块使用 `Δt` |
| `α` `β` `γ` `δ` `λ` | 其他参数 | 见各模块定义（如 PT 系数 α≈1.26、SW 系数 γ 摩阻线性化） |

### 1.2 命名风格

- **函数名**：snake_case 动词或名词短语（`get_soil_params`、`cal_K`、`updatewtd_shallow`、`potevap_shutteworth_wallace`）。
- **就地修改函数**：以 `!` 结尾（`lateral_flow!`、`gw2river!`、`rivers_kw_flood!`）。
- **类型/Struct**：PascalCase（`SoilType`）。
- **模块**：PascalCase（`module ASAP`、`module Evapotranspiration`）。
- **常量**：UPPER_SNAKE（`NVTYP`、`NSTYP`、`POTWILT_LOCAL`、`BIOPARMS`）。
- **形参缩写**：`is/ie/js/je` 网格范围、`nzg` 土壤层数、`dt` 时间步长、`Δt` Unicode 步长、`slz` 层中心深度、`dz` 层厚度、`z₋ₕ` 层下边界深度（地表为 0、向下为负）。
- **下标变量**：`θ_wtd` 水位处含水量、`θ_eq` 平衡含水量、`θ_sat` 饱和含水量、`θ_cp` 凋萎/残余含水量。

### 1.3 维数与缩写

- 矩阵 `M<:Matrix{T}`、三维数组 `A3<:Array{T,3}`、向量 `V<:Vector{T}`、`T<:Real`（多数实现为 `Float64`）。
- `landmask::Matrix{Int}`（0 = 海洋，1 = 陆地）。
- `fdepth` 根系深度因子（m）；同时用于 Ksat/θsat/ψsat 的 `exp(z/fdepth)` 缩放。
- `freedrain::Bool`（Julia 风格）或 `freedrain::Int`（Fortran 兼容）；`updatewtd_shallow` 仍接收 `Int`。

## 2. 单位约定

| 量 | 单位 | 说明 |
|---|---|---|
| 长度 | **m（米）** | 深度、地形、水位、流速距离均为 m |
| 时间 | **s（秒）** | 时间步长、累计时间均为 s |
| 面积 | m² | 单元面积、河宽 × 河长 |
| 流量 | **m³/s** | 河道出口流量 `qnew` |
| 降水 | mm | 单步降水量 `precip`、穿透降水 `ppdrip` |
| 蒸散发 | mm | `et_s`、`et_i`、`et_c`；由 m × 1e3 转换 |
| 体积含水量 | m³/m³ | 无量纲比值 |
| 导水率 | m/s 或 m²/s | 视上下文：饱和值 Ksat 为 m/s；侧向流量 κlat 为 m/s；达西通量时 m²/s |
| 温度 | K | 开尔文；地表温度、参考气温 |
| 压力 | Pa | 大气压；模块内 `presshp` 为 hPa 时显式标注 |
| 同位素 ¹⁸O | ‰ 或 m³/m³ | 见 `IsotopeTracing` 约定 |

**m → mm 换算**：所有从 `soilfluxes`、`rivers_*`、`updatewtd_shallow` 等返回的米制通量在 `RootDepth-主算法` 中通过 `* 1.0e3` 转为 mm 再累加到状态变量。**反向 mm → m 用 `* 1.0e-3`**。

## 3. 数组索引约定

- 全部数组遵循 **Julia 1-based** 约定（第 1 个元素索引为 1）。
- 网格范围 `is:ie`（i 方向）、`js:je`（j 方向）；**主循环跳过 `is+1:ie-1` 与 `js+1:je-1` 的边界 halo**。
- 土壤层索引 `k = 1` 对应最深层，`k = nzg` 对应表层；`nzg+1` 层节点 `z₋ₕ[nzg+1] = 0`（地表）。
- 地下水位所在层 `jwt = find_jwt(wtd, z₋ₕ)`，`jwt == 1` 表示地下水位低于最深层。
- `Array{Int,3}` 中 `[2, is:ie, js:je]` 第 1 维为土壤层序号（1/2），第 2/3 维为网格。

## 4. 文档页面模板

每个模块/函数页面统一采用以下 7 节结构（与阶段 A 已建页面保持一致）：

```markdown
# <模块名> — <中文副标题>

> 源文件：src/<file>.jl:L1-L<max>
> Fortran 来源：fortran/<file>.f90:<subroutine>
> 测试：test/<test_file>.jl
> 状态：已摄取 | 已摄取（含悬空 export 备注）| 待复核

## 1. 功能概述
（≤ 5 行的功能定位，包含模块在整体流程中的位置）

## 2. 函数签名
（Julia 函数签名块，可附返回值说明）

## 3. 算法 / 公式
（LaTeX 公式 + 算法分节描述；引用 docs/*.typ 中的推导）

## 4. 关键变量与单位
（Markdown 表格：符号 / 含义 / 单位）

## 5. 与 Fortran 对照
（表格：Fortran 子程序 / Julia 函数 / 差异）

## 6. 引用
（行号、测试断言、文档参考）

## 7. 已知问题与备注
（悬空 import、悬空 export、未导出函数、拼写错误等）
```

## 5. 状态标记

- **已摄取**：源码、签名、公式、引用四节齐备。
- **已摄取（含悬空 export 备注）**：基本齐备但有遗留问题已登记在第 7 节。
- **待复核**：源码已读但存在歧义（如 `updatewtd_shallow.flag` 未初始化、`extraction.kroot = k - 1` 行为待与 Fortran 比对），需下一轮核实。

## 6. 引用风格

- 行号引用：`src/<file>.jl:L<start>-L<end>`，例如 `src/SoilFluxes.jl:L82-L134`。
- 测试引用：`test/<file>.jl:L<n>`，写明断言内容。
- 文档引用：`docs/<file>.typ §<n>`（Typst 文档），引用其中的公式编号或章节。