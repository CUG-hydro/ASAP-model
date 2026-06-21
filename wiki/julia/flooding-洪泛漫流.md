# flooding! — 洪泛区漫流

> 源文件：`src/modules/flooding.jl`:L1-L70
> Fortran 来源：`fortran/module_wtable.f90:L1226-L1322`（`subroutine FLOODING`）
> 测试：`test/wtable/test_river_routing.jl`（"Flooding Calculation"）
> 状态：已摄取

## 1. 功能概述
对洪泛平原 (floodplain) 上的积水按 8 邻域最低高程方向进行漫流扩散：寻找最低邻格、计算坡降、并按对角线距离 (×1/√2) 调整；主河道方向优先，且当 `delsfcwat < 0` 时受地表水变化约束。所有更新汇集到 `delsfcwat`，时间分割 `ntsplit = 1`。

## 2. 函数签名
```julia
function flooding!(
  imax::Int, jmax::Int, is::Int, ie::Int, js::Int, je::Int, Δt::T,
  fd::Matrix{Int}, bfd::Matrix{Int}, topo::M, area::M,
  riverwidth::M, riverlength::M, riverdepth::M,
  floodheight::M, delsfcwat::M
) where {T<:AbstractFloat, M<:Matrix{T}}
```

## 3. 算法 / 公式

**邻格高差**（含对角线距离调整）：

$$
dh = \bigl(\mathrm{fh}_{i,j} + \mathrm{topo}_{i,j}\bigr) - \bigl(\mathrm{fh}_{ii,jj} + \mathrm{topo}_{ii,jj}\bigr)
$$

$$
dh \mathrel{/}= \sqrt{2} \quad \text{当}\ ii\ne i\ \text{且}\ jj\ne j\ \text{（对角线邻）}
$$

**漫流通量**：

$$
d_{ij} = \max\!\left(\mathrm{fh}_{i,j} - \max\!\left(\tfrac{1}{2}\bigl(\mathrm{topo}_{ilow,jlow} - \mathrm{topo}_{i,j} + d_{total}\bigr),\ 0\right),\ 0\right)
$$

**主河道扣除**（仅当最低邻即下游流向 `(i1,j1)`）：

$$
d_{ij} \mathrel{-}= \frac{W_r \cdot \mathrm{fh}_{i,j} \cdot L_r}{A}
$$

**地表水亏空约束**：

$$
d_{ij} = \max\!\bigl(\min(d_{ij},\ \mathrm{fh}_{i,j} + \mathrm{delsfcwat}_{i,j}),\ 0\bigr) \quad \text{当}\ \mathrm{delsfcwat}_{i,j} < 0
$$

## 4. 关键变量与单位
| 符号 | 含义 | 单位 |
|---|---|---|
| floodheight | 洪泛区水深 | m |
| topo | 地形高程 | m |
| delsfcwat | 地表水增量 Δ | m |
| dflood | 局部漫流累积 | m |
| riverwidth / riverlength | 河道宽 / 长 | m / m |
| area | 单元面积 | m² |
| ntsplit | 时间分割数 | — |

## 5. 与 Fortran 对应
| Fortran 子程序 | Julia 函数 | 差异 |
|---|---|---|
| `FLOODING` (L1226) | `flooding!` | 主体移植；`ntsplit` 显式赋 1（Fortran 中可外部传入） |
| 内嵌循环 8 邻域搜索 | 8 邻域方向数组 | 保留原 D8 邻格遍历顺序 |

## 6. 引用
- 行号：`src/modules/flooding.jl`:L12 初始化 `dflood`；L18-L24 邻格高差与对角线因子；L26-L32 主河道扣除；L34-L36 地表水约束
- 测试断言：`test/wtable/test_river_routing.jl` 中 `@test_nowarn flooding!`；`@test all(isfinite.(floodheight))`；`@test abs(sum(delsfcwat)) < 1e-10`（水量守恒）
- 文档：`docs/汇流_扩散波.typ`（圣维南方程组给出洪泛区漫流的连续性方程基础）

## 7. 已知问题与备注
- 时间步 `ntsplit = 1` 固定，无法对陡坡洪泛区做次步长细分
- 仅按单步 D8 漫流，无显式扩散项，平原区扩散偏慢
- 漫流量采用 0.5 系数做相邻单元平均（`0.5*(topo[ilow,jlow] - topo[i,j] + dtotal)`），缺乏严格推导
- `delsfcwat < 0` 截断只约束源端，未在目标端作相应扣除，可能在小幅破坏水量守恒
- 函数未使用 `bfd` 反向流向参数，存在悬空参数
