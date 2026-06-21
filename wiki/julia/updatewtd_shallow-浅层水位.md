# updatewtd_shallow — 浅层水位更新

> 源文件：src/updatewtd_shallow.jl:L1-L67
> 测试：test/test_updatewtd_shallow.jl
> 状态：已摄取

## 1. 功能概述

本模块实现浅层地下水位（WTD，Water Table Depth）的逐时间步更新。在自由排水（freedrain）或浅埋条件下，基于各土壤层含水量与平衡含水量之差判断水位升降方向，并依据质量守恒计算水位的离散跳跃与释放/吸收的补给量 `rech`。

## 2. 函数签名

```julia
function cal_factor(z::T, fdepth::T) where {T<:Real}
function updatewtd_shallow(
    nzg::Int, freedrain::Int, z₋ₕ::Vector{Float64},
    dz::Vector{Float64}, soiltxt::Int, θ_eq::Vector{Float64},
    θ_wtd::Float64, θ::Vector{Float64}, wtd::Float64, fdepth::Float64
)
```

返回值：`Tuple{Float64, Float64}`，依次为更新后的地下水位深度 `wtd` 与补给量 `rech`。

## 3. 算法 / 公式

饱和含水量的深度衰减因子：

$$
\mathrm{cal\_factor}(z, f_\mathrm{depth}) = \mathrm{clamp}\!\left(\exp\!\left(\frac{z + 1.5}{f_\mathrm{depth}}\right),\ 0.1,\ 1.0\right)
$$

水位上升的补给释放量：

$$
\mathrm{rech} = (w_\mathrm{td}^{\mathrm{old}} - w_\mathrm{td}) \cdot (\theta_\mathrm{sat}(z_\mathrm{jwt}) - \theta_\mathrm{eq}(k_\mathrm{wt}))
$$

通过 `find_jwt(wtd, z₋ₕ)` 定位水位所在层索引，再依据 $\theta_{kwt} > \theta_{eq,kwt}$ 判定上升，若饱和则整层抬升一层。

## 4. 关键变量与单位

| 符号 | 含义 | 单位 |
|---|---|---|
| `wtd` | 地下水位深度（负值到地表为 0） | m |
| `fdepth` | 根系/包气带深度因子 | m |
| `θ_eq` | 各层平衡含水量 | m³/m³ |
| `θ_sat` | 饱和含水量（深度修正后） | m³/m³ |
| `z₋ₕ` | 层边界深度数组 [nzg+1] | m |
| `rech` | 释放/吸收补给量 | m |

## 5. 与 Fortran 对应

| Fortran 子程序 | Julia 函数 | 差异 |
|---|---|---|
| `WTD::shallow_update` | `updatewtd_shallow` | 命名采用小写下划线；内联 `find_jwt`；输出改为多返回值 |
| — | `cal_factor` | 工具函数，无 Fortran 直接对应 |

## 6. 引用

- 行号：src/updatewtd_shallow.jl:L1-L67（核心循环 L33-L60）
- 测试断言：test/test_updatewtd_shallow.jl:L6-L26（`find_jwt` 边界）、L52-L62（`cal_factor` 指数衰减）
- `cal_factor(-1.5, fdepth) ≈ 1.0` 验证 $z=-1.5$ 时因子归一
- `cal_factor(-100, 2) ≈ 0.1` 验证深度截断下限

## 7. 已知问题与备注

- `freedrain` 参数已声明但当前未在主循环内使用。
- 函数体中存在 `if flag == 0` 比较，但 `flag` 在浅层饱和检测中未必被赋值，存在局部变量未初始化的隐患。
- 未导出 `cal_factor`，但测试通过 `ASAP.cal_factor` 访问，依赖模块的隐式可见性。
