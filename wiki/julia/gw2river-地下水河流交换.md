# gw2river! — 地下水-河流交换

> 源文件：`src/modules/gw2river.jl`:L1-L125
> Fortran 来源：`fortran/module_wtable.f90:L744-L796`（`subroutine gw2river`）、`L1325-L1382`（`subroutine moveqrf`）
> 测试：`test/wtable/test-river.jl`、`test/wtable/test_river_routing.jl`
> 状态：已摄取

## 1. 功能概述
本模块计算地下水位 (wtd) 与河流水面之间的双向通量 `qrf`，并通过 `moveqrf!` 将小河流网格的通量重新分配到下游河道，最后由 `apply_specific_diversions` 应用流域内若干硬编码的人工分流（Taquari 河等）。

## 2. 函数签名
```julia
function gw2river!(
  imax::Int, jmax::Int, is::Int, ie::Int, js::Int, je::Int, nzg::Int,
  slz::V, Δt::T, soiltxt::Array{Int,3},
  landmask::Matrix{Int}, wtd::M, maxdepth::M,
  riverdepth::M, width::M, length::M,
  area::M, fdepth::M, qrf::M
) where {T<:AbstractFloat,V<:Vector{T},M<:Matrix{T}}

function moveqrf!(imax::Int, js::Int, je::Int, fd::Matrix{Int},
                  qrf::M, area::M, width::M) where {T<:AbstractFloat,M<:Matrix{T}}

function apply_specific_diversions(i::Int, j::Int, q::M, dsnew::T)
  where {T<:AbstractFloat,M<:Matrix{T}}
```

## 3. 算法 / 公式

河床导水率沿深度按指数衰减（Miguez-Macho 等, 2007, Eq. 2d）：

$$
K_{rb} = K_{sat} \cdot \mathrm{clamp}\!\left(\exp\!\left(\frac{-\mathrm{maxdepth}+1.5}{\mathrm{fdepth}}\right),\ 0.1,\ 1.0\right)
$$

$$
q_{rf} =
\begin{cases}
K_{rb}\, W L\,(wtd - z_{riv})\,\dfrac{\Delta t}{A}, & wtd > z_{riv} \quad (\text{GW} \to \text{River}) \\[6pt]
-\min\!\left(K_{rb}\, W L\,(z_{riv}-wtd)\,\dfrac{\Delta t}{A},\ r_d\right)\min\!\left(\dfrac{WL}{A},1\right), & -d_{max} < wtd \le z_{riv} \\[6pt]
-K_{sat}\,\Delta t, & wtd \le -d_{max} \quad (\text{仅渗透})
\end{cases}
$$

每日排水上限 50 mm（即 `min(..., Δt·0.05/86400)`），避免突变。

## 4. 关键变量与单位
| 符号 | 含义 | 单位 |
|---|---|---|
| wtd | 地下水位（地表为 0，向上为正） | m |
| maxdepth | 河道最大深度 | m |
| riverdepth | 河道当前水深 | m |
| fdepth | 河床导水率衰减深度 | m |
| width/length/area | 河宽/河长/单元面积 | m / m / m² |
| qrf | 河流-地下水通量 | m |
| Ksat | 饱和水力传导度 | m s⁻¹ |

## 5. 与 Fortran 对应
| Fortran 子程序 | Julia 函数 | 差异 |
|---|---|---|
| `gw2river` (L744) | `gw2river!` | 同名移植，类型参数化；保留 `K_rb` 指数衰减与 50 mm/day 上限 |
| `moveqrf` (L1325) | `moveqrf!` | 直接移植小河流 → 下游重分配逻辑 |
| `RIVERS_KW_FLOOD` 中的硬编码分流 | `apply_specific_diversions` | 抽成独立函数，按 `(i,j)` 索引加减下游格点流量 |

## 6. 引用
- 行号：`src/modules/gw2river.jl`:L1-L125（`gw2river!`、L26-L53 主体；L54-L80 `moveqrf!`；L82-L125 `apply_specific_diversions`）
- 测试断言：`test/wtable/test-river.jl` 中 `@test all(isfinite.(qrf))` 与 `@test any(qrf .> 0.0)`
- 测试断言：`test/wtable/test_river_routing.jl` 中 `@test_nowarn moveqrf!` 与 `@test qrf[i,j] == 0.0`（小河流）
- 文档：`docs/汇流_扩散波.typ` 提供圣维南方程组动量守恒背景

## 7. 已知问题与备注
- `K_rb` 公式注释中指出"未做积分，认为 dl = dh"——属于强简化版本，仅在 losing river 严格成立
- 50 mm/day 的排水上限可能对洪泛平原快速排泄过强
- `apply_specific_diversions` 中硬编码网格坐标（Taquari 河 4498/4535 等）仅适用于特定流域，移植到其他流域需清空
- `width` 在 `gw2river!` 内既作为 `Matrix{T}` 又在循环中与 `length`/`area` 同名局部读取，需注意 `Base.length` 的命名冲突
