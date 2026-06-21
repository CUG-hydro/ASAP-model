# IsotopeTracing — 同位素（¹⁸O）侧向流与深层水位追踪

> 源文件：src/modules/Tracing/IsotopeTracing.jl:L1-L209
> 测试：test/wtable/test_isotope_tracing.jl
> 状态：已摄取

## 1. 功能概述

本模块在侧向地下水流的基础上耦合稳定同位素 ¹⁸O 的对流传输。它将 8 邻居差分简化为带纬度修正的 4 主方向格式（北/南/西/东），并配套提供深层地下水位更新与简化版 `updatewtd_simple`，用于大规模耦合模拟中的深层含水层估算。

## 2. 函数签名

```julia
function lateral_isotope!(
    imax, jmax, is, ie, js, je, nzg,
    soiltxt::Array{Int,3}, wtd::M, qlat::M, fdepth::M, topo::M, landmask,
    Δt::T, area::M, lats::M, dxy::T,
    slz::Vector{T}, o18::Array{T,3}, smoi::Array{T,3},
    qlato18::M, qlatin::M, qlatout::M,
    qlatino18::M, qlatouto18::M,
    qlatinsum::M, qlatoutsum::M,
    qlatino18sum::M, qlatouto18sum::M
) where {T<:AbstractFloat, M<:Matrix{T}}

function updatedeepwtable!(
    imax, jmax, js, je, nzg,
    slz::Vector{T}, dz::Vector{T}, soiltxt::Array{Int,3},
    wtd::M, bottomflux::M, rech::M,
    qslat::M, qlat::M, landmask::Matrix{Int},
    Δt::T, smoi::Array{T,3}, smoieq::Array{T,3},
    smoiwtd::M, qsprings::M
)
```

## 3. 算法 / 公式

**水位层同位素浓度归一化**（按含水量除算为液态浓度）：

$$
c_\mathrm{wtd}^{(i,j)} = \frac{o_{18}^{(k_\mathrm{wtd}, i, j)}}{\mathrm{smoi}^{(k_\mathrm{wtd}, i, j)}}
$$

**四方向通量（含纬度 cos 修正）**：

$$
q_N = \frac{1}{2}(\kappa_i + \kappa_{j+1})\,(h_{j+1} - h_j)\cos\!\bigl(\pi_2^\mathrm{rad}(l + \tfrac{\Delta_y}{2})\bigr)
$$

$$
q_W = \frac{1}{2}(\kappa_i + \kappa_{i-1})\,(h_{i-1} - h_i)\,/\cos(\pi_2^\mathrm{rad}\,l)
$$

其中 $\pi_2^\mathrm{rad} = \pi/180$，$\Delta_y$ 为格距。

**深层水量平衡**：

$$
\mathrm{totwater} = q_\mathrm{lat} - q_\mathrm{slat} - \mathrm{deeprech}
$$

`rech` 在浅层（$wtd < \mathrm{slz}[1] - \mathrm{dz}[1]$）由排导率与底通量差分得到。

## 4. 关键变量与单位

| 符号 | 含义 | 单位 |
|---|---|---|
| `o18` | ¹⁸O 含量（三维：层×i×j） | ‰ 或 m³/m³（取决于约定） |
| `o18wtd` | 水位处液态同位素浓度 | ‰ |
| `lats` | 网格中心纬度 | ° |
| `dxy` | 网格间距 | ° |
| `qlatinsum` 等 | 累积入/出流量 | mm |
| `totwater` | 深层格点净水量 | m |

## 5. 与 Fortran 对应

| Fortran 子程序 | Julia 函数 | 差异 |
|---|---|---|
| `IsotopeFlow::lateral_iso` | `lateral_isotope!` | 显式参数列表；4 方向而非 8 方向 |
| `WTD::deep_update` | `updatedeepwtable!` | 引入 `bottomflux` 显式回写 |
| — | `updatewtd_simple` | 简化版本，原 Fortran 未单列子程序 |

## 6. 引用

- 入口函数：src/modules/Tracing/IsotopeTracing.jl:L13-L57 `lateral_isotope!`
- 四方向计算：L101-L120（北/南/西/东）
- 深层更新：L122-L209 `updatedeepwtable!`
- 简化版：L186-L260 `updatewtd_simple`
- 测试断言：test/wtable/test_isotope_tracing.jl:L86-L116（**Isotope Conservation** 质量守恒）、L46-L73（梯度传输）

## 7. 已知问题与备注

- `π2r` 定义为模块级常量，仅供本文件内部使用，未导出。
- 4 方向格式不包含对角线（区别于 `lateral_flow!` 的 8 邻居），需保证调用方不依赖对角线贡献。
- 累积数组 `*sum` 在每次调用末尾乘以 `1e3` 转换 m→mm，跨调用会持续累加。
- `updatewtd_simple` 中部分 `soil.θ_sat` / `soil.θ_cp` 字段依赖外部 `soil` 全局，建议替换为参数注入。
