# rivers_kw_flood! — 运动波河道洪水路由

> 源文件：`src/modules/rivers_kw_flood.jl`:L1-L101
> Fortran 来源：`fortran/module_wtable.f90:L800-L993`（`subroutine RIVERS_KW_FLOOD`）
> 测试：`test/wtable/test_river_routing.jl`（"Kinematic Wave Flood Routing"）
> 状态：已摄取

## 1. 功能概述
采用一维运动波（Kinematic Wave）方程对河道-洪泛区进行洪水演算：使用 Manning 公式估算流速，按 `fd` 流向在每一步向河道汇入侧向径流（`qrf`、`qs`、`delsfcwat`），并基于 `riverchannel` 阈值将水深在河道与洪泛区之间重分配，最后累加 `qmean`。

## 2. 函数签名
```julia
function rivers_kw_flood!(
  imax::Int, jmax::Int, is::Int, ie::Int, js::Int, je::Int,
  Δt::T, δt::T,
  fd::Matrix{Int}, bfd::Matrix{Int},
  qnew::M, qs::M, qrf::M,
  delsfcwat::M, depth::M,
  riverarea::M, floodarea::M, floodheight::M,
  width::M, length::M, maxdepth::M, slope::M, area::M, topo::M, riverchannel::M,
  qmean::M
) where {T<:AbstractFloat, M<:Matrix{T}}
```

## 3. 算法 / 公式

运动波忽略惯性项与附加比降（参见 `docs/汇流_扩散波.typ` §1 推导），由曼宁公式直接给出流速：

$$
V = \frac{1}{n} R^{2/3}\sqrt{S_0},\qquad
Q = V \cdot A
$$

$$
R = \frac{A}{2y + W},\qquad A = W \cdot y
$$

水面高程比较（驱动水深更新）：

$$
z_0 = \mathrm{topo}_{i,j} - \mathrm{maxdepth}_{i,j} + y_{i,j},\quad
z_1 = \mathrm{topo}_{i_1,j_1} - \mathrm{maxdepth}_{i_1,j_1} + \max(y_{i_1,j_1},0)
$$

河道-洪泛区重分配：当 `riverarea·maxdepth`（即 `riverchannel`）阈值跨越，则把超出 `maxdepth` 的部分从河道水深转入洪泛水深 `floodheight`。

## 4. 关键变量与单位
| 符号 | 含义 | 单位 |
|---|---|---|
| qnew | 网格出口流量 | m³ s⁻¹ |
| depth | 河道水深 | m |
| floodheight | 洪泛区水深 | m |
| riverarea / floodarea | 河道/洪泛面积 | m² |
| riverchannel | 阈值 `maxdepth·riverarea` | m³ |
| slope | 河床比降 `S₀` | — |
| width / length | 河宽 / 河长 | m / m |
| n | 曼宁糙率 | s m⁻¹/³ |
| qmean | 时间累加平均流量 | m³ |

## 5. 与 Fortran 对应
| Fortran 子程序 | Julia 函数 | 差异 |
|---|---|---|
| `RIVERS_KW_FLOOD` (L800) | `rivers_kw_flood!` | 接口完全对齐；将 `deltat`/`dtlr` 拆为 `Δt`/`δt` 参数语义化 |
| 内联硬编码分流 | `apply_specific_diversions` | 抽成可复用的独立函数 |
| 局部 `inflow` 数组 | `qin = zeros(size(q))` | 名称 `qin` 取代原 Fortran 变量 |

## 6. 引用
- 行号：`src/modules/rivers_kw_flood.jl`:L13-L23 计算外部入流；L26-L32 流向入流；L34-L46 Manning 公式主体
- 测试断言：`test/wtable/test_river_routing.jl` 中 `@test_nowarn rivers_kw_flood!` 与 `@test all(qnew .>= 0.0)`、`@test all(depth .>= 0.0)`
- 文档：`docs/汇流_扩散波.typ` §1（运动波 vs 扩散波对比）§3（圣维南方程组线性化推导）

## 7. 已知问题与备注
- 运动波假设在缓坡河漫滩失效，应让 `rivers_dw_flood!` 在洪水态接管（见 `rivers_dw_flood-扩散波路由.md`）
- `length` 参数与 `Base.length` 同名，存在命名冲突隐患
- `riverchannel` 在 Fortran 中是 `maxdepth·riverarea`（体积阈值），在 Julia 中保留同样语义
- `n` 曼宁糙率取默认 0.03，未在参数列表暴露
