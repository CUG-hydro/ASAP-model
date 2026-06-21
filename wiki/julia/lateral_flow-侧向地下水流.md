# lateral_flow! — 侧向地下水流（D8 八邻居差分格式）

> 源文件：src/modules/lateral_flow.jl:L1-L60
> 测试：test/wtable/test_watertable.jl
> 状态：已摄取

## 1. 功能概述

基于达西定律和 D8（八方向）三角形有限差分格式，计算每个陆面网格与 8 个邻居之间的侧向地下水交换通量。深度依赖的有效导水率随地下水位埋深呈指数衰减，并按对角线方向除以 $\sqrt{2}$ 修正。

## 2. 函数签名

```julia
function lateral_flow!(
    imax::Int, jmax::Int, is::Int, ie::Int, js::Int, je::Int,
    wtd::M, qlat::M, fdepth::M,
    topo::M, landmask::Matrix{Int}, Δt::T,
    area::M, κlat::M
) where {T<:AbstractFloat, M<:Matrix{T}}
```

就地更新 `qlat`（净侧向流量，m/s），并写回通量结果。

## 3. 算法 / 公式

**达西通量**（沿 $\xi$ 方向单位网格边长）：

$$
q_\xi = \frac{1}{2}\bigl(\kappa_\mathrm{cell}^\mathrm{self} + \kappa_\mathrm{cell}^\mathrm{neighbor}\bigr)\cdot\frac{h_\mathrm{neighbor} - h_\mathrm{self}}{d_\xi}
$$

对角线邻居：$d_\xi = \sqrt{2}$，并乘以角度因子

$$
f_\mathrm{angle} = \frac{\sqrt{\tan(4\pi/32)}}{2\sqrt{2}}
$$

**有效导水率**（按水位埋深分支）：

$$
\kappa_\mathrm{cell} =
\begin{cases}
0, & f_\mathrm{depth} < 10^{-6} \\
f_\mathrm{depth}\,\kappa_\mathrm{lat}\,\exp\!\bigl((wtd+1.5)/f_\mathrm{depth}\bigr), & wtd < -1.5 \\
\kappa_\mathrm{lat}\,(wtd + 1.5 + f_\mathrm{depth}), & \text{otherwise}
\end{cases}
$$

## 4. 关键变量与单位

| 符号 | 含义 | 单位 |
|---|---|---|
| `wtd` | 地下水位深度（负值） | m |
| `topo` | 地形高程 | m |
| `fdepth` | 包气带深度因子 | m |
| `κlat` | 侧向饱和水力传导度 | m/s |
| `κcell` | 有效导水率 | m²/s |
| `head = topo − wtd` | 地下水水头 | m |
| `qlat` | 净侧向流量（输出） | m/s |

## 5. 与 Fortran 对应

| Fortran 子程序 | Julia 函数 | 差异 |
|---|---|---|
| `LateralFlow::compute` | `lateral_flow!` | 接受类型参数 `M<:Matrix{T}`；`qlat` 改为就地修改 |
| — | `flowdir` | 新增 D8 流向辅助函数，位于 test 引用处 |

## 6. 引用

- 行号：src/modules/lateral_flow.jl:L1-L60
- 角度因子：L7 `fangle = sqrt(tan(4π/32)) / (2√2)`
- 有效导水率：L21-L27（深度分支判断）
- 八方向累加：L29-L46（4 对角 + 4 邻接）
- 测试断言：test/wtable/test_watertable.jl:L7-L18（`flowdir` 八方向）；L23-L40（`qlat` 有限性）

## 7. 已知问题与备注

- 文件顶部包含 docstring 但与 `Modules.jl` 中的 `export` 重复声明，无需额外处理。
- 陆面掩膜 `landmask` 仅在主循环判断，未对边缘网格（`is/ie/js/je` 边界）做 halo 保护；外层调用方需保证 padding。
- 暂未实现时间步积分（`Δt` 形参未参与计算），需后续与 `qsprings` 协调。
