# SoilInitialization - 土壤分层

> 源文件：src/SoilInitialization.jl:L1-L81
> Fortran 来源：fortran/module_initial.f90（CLM 指数分层、预设 DZ2 分层）
> 测试：test/test_soil_initialization.jl
> 状态：已摄取

## 1. 功能概述
该模块负责 ASAP 模型土壤柱的垂向离散化，提供两套土壤分层方案：基于 CLM（Community Land Model）指数增长节点深度的 `initializesoildepth_clm`，以及按 40 层预定义厚度（`DZ2`）直接查表的 `initializesoildepth`。函数返回层中心深度 `slz`（或节点深度 `z₋ₕ`）与层厚度 `dz`，统一约定地表以下为负、地表以上为正。

## 2. 函数签名
```julia
function initializesoildepth_clm(nzg::Int)::NamedTuple{(:slz, :dz), NTuple{2,Vector{Float64}}}
function initializesoildepth(nzg::Int)::Tuple{Vector{Float64}, Vector{Float64}}
```

## 3. 算法 / 公式

### CLM 指数分层
节点下边界采用指数展开：
$$z_k = 0.025 \left(e^{0.5(k-0.5)} - 1\right), \quad k = 1,\dots,n_{zg}$$

层厚度取相邻节点中点差分：
$$dz_k = \tfrac{1}{2}(z_{k+1} - z_{k-1}),\quad k=2,\dots,n_{zg}-1$$
端点：$dz_1 = \tfrac{1}{2}(z_1 + z_2)$，$dz_{n} = z_n - z_{n-1}$。

层中心：$slz_k = \tfrac{1}{2}(z_{k-1} + z_k)$，$k=2,\dots,n$；$slz_1 = \tfrac{1}{2} z_1$。

### 固定层厚分层
直接以 `DZ2` 数组倒数第 `n` 个元素为顶层（地表），逐层向下累加：
$$z_{k} = z_{k+1} - dz_k, \quad dz_k = DZ2[n_{zg}-k+1]$$

## 4. 关键变量与单位
| 符号 | 含义 | 单位 |
|---|---|---|
| nzg | 土壤层数 | — |
| slz / z₋ₕ | 节点深度（地表为 0，地下为负） | m |
| dz | 层厚度 | m |
| DZ2 | 40 层预设厚度查找表（0.1–540 m） | m |
| vctr4 | CLM 节点下边界中间量 | m |

## 5. 与 Fortran 对应
| Fortran 子程序 | Julia 函数 | 差异 |
|---|---|---|
| CLM `setsoildepth`/`init_clm_soil` 指数节点 | `initializesoildepth_clm` | Julia 显式生成 `dz` 与 `slz`，Fortran 内联在初始化块 |
| 预设 DZ2 查表初始化 | `initializesoildepth` | 公式一致；Julia 显式校验 `nzg ≤ length(DZ2)` |
| 输出 `slz, dz` 反转顺序 | 同名 | Julia 使用 `reverse` 一次性反转向量；Fortran 多采用 `do` 循环逆序填入 |

## 6. 引用
- 常量 DZ2：src/SoilInitialization.jl:L4-L8
- 函数 initializesoildepth_clm：src/SoilInitialization.jl:L26-L52
- 函数 initializesoildepth：src/SoilInitialization.jl:L65-L80
- 测试断言：test/test_soil_initialization.jl:L6 `length(z₋ₕ) == nzg+1`，L7 `length(dz) == nzg`，L14 `z₋ₕ[1] ≈ -1000`

## 7. 已知问题与备注
- 注释中明确指出：**该算法存在问题，slz 计算到了 N+1 层**（src/SoilInitialization.jl:L23-L24, L48）。`slz` 长度为 `nzg+1`，但 `dz` 仅有 `nzg` 个元素，最后一项 `slz[nzg+1] = vctr4[nzg] + 0.5*dz[nzg]` 使用了未与 `dz` 一一对应的下界值，建议后续对齐两数组长度或补充 `dz[nzg+1]` 哨兵。
- `using DataFrames`（L1）在当前实现中未被使用，属于悬空 import，调用者无需引入 DataFrames 即可使用本模块。
- `initializesoildepth_clm` 返回 `NamedTuple` `(slz=..., dz=...)`，而 `initializesoildepth` 返回普通 `Tuple`，API 不一致；下游调用需注意解构方式。
- 反转后 `dz` 顺序与 `DZ2` 表底→表顶顺序一致（最厚层位于最深处），与地表向上为正、向下为负的符号约定一致。
- 文档：docs/土壤水运动_ASAP.typ 给出分层原则与 ASAP 默认 `nzg` 取值。
