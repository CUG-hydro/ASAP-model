# rivers_dw_flood! — 扩散波河道洪水路由

> 源文件：`src/modules/rivers_dw_flood.jl`:L1-L100
> Fortran 来源：`fortran/module_wtable.f90:L996-L1222`（`subroutine RIVERS_DW_FLOOD`）
> 测试：`test/wtable/test_rivers_dw_flood.jl`
> 状态：已摄取

## 1. 功能概述
对河漫滩洪水采用扩散波（Diffusion Wave）方程进行路由：洪水态 (`floodheight > 0.05` 或下游超载) 启用扩散波隐式更新；非洪水态回退至运动波。曼宁糙率固定 `n = 0.03`，流速使用 `clamp(speed, 0.01, length/δt)` 防止 CFL 失稳。

## 2. 函数签名
```julia
const g0 = 9.81  # 重力加速度 (m/s²)

function rivers_dw_flood!(
  imax::Int, js::Int, je::Int, Δt::T, δt::T,
  fd::Matrix{Int}, bfd::Matrix{Int}, qnew::M,
  qs::M, qrf::M, delsfcwat::M,
  slope::M, depth::M, width::M,
  length::M, maxdepth::M, area::M,
  riverarea::M, floodarea::M, riverchannel::M,
  qmean::M, floodheight::M, topo::M
) where {T<:AbstractFloat, M<:Matrix{T}}
```

## 3. 算法 / 公式

**扩散波**（保留压力项、忽略对流惯性，对应 `docs/汇流_扩散波.typ` §3 推导）：

$$
q^{n+1} = \frac{q^n - g\, y^n \Delta t\, S_h^n}{1 + g\, \Delta t\, n^2\, q^n / \left(R^{4/3}\, y^n\right)}
$$

水力半径 `R = A / (2y + W)`；水面比降 `S_h = (z_0 - z_1)/L`。冻结系数线性化（frozen-coefficient linearization）将非线性摩阻项近似为 `γⁿ·qⁿ⁺¹`。

**运动波回退**（非洪水态）：

$$
V = \frac{1}{n} R^{2/3}\sqrt{S_0},\qquad Q = V \cdot A
$$

**水力半径与 `width = 0` 兜底**：

$$
R = \begin{cases}
A/(2y + W), & W > 0\\
y, & W = 0
\end{cases},\qquad
\mathrm{flowwidth} = \begin{cases}
W, & W \ne 0\\
\sqrt{A}, & W = 0
\end{cases}
$$

## 4. 关键变量与单位
| 符号 | 含义 | 单位 |
|---|---|---|
| qnew | 网格出口流量 | m³ s⁻¹ |
| depth | 河道水深 `y` | m |
| floodheight | 洪泛区水深 | m |
| slope | 河床比降 `S₀` | — |
| g0 | 重力加速度 | m s⁻² |
| n | 曼宁糙率（固定 0.03） | s m⁻¹ᐟ³ |
| γ | 摩阻线性化系数 | m⁻¹ |
| qmean | 时间累加平均流量 | m³ |

## 5. 与 Fortran 对应
| Fortran 子程序 | Julia 函数 | 差异 |
|---|---|---|
| `RIVERS_DW_FLOOD` (L996) | `rivers_dw_flood!` | 接口完全对齐；`g0` 提升为模块级 `const` |
| 内联硬编码分流 | `apply_specific_diversions` | 共享 `gw2river.jl` 中的实现 |

## 6. 引用
- 行号：`src/modules/rivers_dw_flood.jl`:L1 `g0`；L19-L26 入流汇入；L29-L46 扩散波/运动波分支；L48 `qmean` 累加
- 测试断言：`test/wtable/test_rivers_dw_flood.jl` 中 `@test_nowarn rivers_dw_flood!` 洪水/非洪水两套用例；`@test qnew != qnew_before`；`@test qmean ≈ qmean_before .+ qnew .* dtlr atol=1e-8`
- 文档：`docs/汇流_扩散波.typ` §3（含 boxed 公式 `q^{n+1}` 完整推导）

## 7. 已知问题与备注
- 曼宁糙率 `n = 0.03` 硬编码，未作为参数暴露，无法按河段/植被差异化
- `width == 0` 时使用 `√area` 作为等效流宽，属于经验兜底，缺乏物理依据
- `speed` 上下限 `clamp(..., 0.01, length/δt)` 隐含 CFL 限制，但阈值 0.01 m/s 在山区可能过低
- 共享 `apply_specific_diversions` 来自 `gw2river.jl`，模块间存在隐式耦合
- 与 `rivers_kw_flood!` 存在功能重叠：同一时间步应仅调用其中一个（Fortran 中由外部开关控制）
