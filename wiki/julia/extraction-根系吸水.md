# extraction — 根系吸水与冠层蒸腾

> 源文件：src/extraction.jl:L45-L186
> Fortran 来源：fortran/extraction.f90:extraction（即将映射，待阶段2补充）
> 测试：test/test_extraction.jl
> 状态：已摄取

## 1. 功能概述

`extraction` 是 ASAP 模型冠层-土壤耦合的核心函数，依据 Penman-Monteith 方程计算潜在冠层蒸腾（`pet_c`）与土壤蒸发（`pet_s`），并按根系活性权重从各土壤层提取水分，返回各层含水量变化 `dθ`、水分亏缺 `watdef` 及深层地下水变化 `dθ_deep`。函数同步维护每层的非活跃天数 `inactivedays`，实现对冰冻层、极小活性层与连续干旱层的自适应抑制。

## 2. 函数签名

```julia
function extraction(
  nzg::Int, z₋ₕ::Vector{Float64}, dz::Vector{Float64},
  Δt::Float64, soiltxt::Int,
  wtd::Float64, θ::Vector{Float64}, θ_wtd::Float64,
  Δ::Float64, γ::Float64, λ::Float64,
  lai::Float64, ra_a::Float64, ra_c::Float64, rs_c_factor::Float64,
  R_a::Float64, R_s::Float64, petfactor_s::Float64, petfactor_c::Float64,
  inactivedays::Vector{Int}, maxinactivedays::Int,
  hhveg::Float64, fdepth::Float64, icefac::Vector{Int8}
) :: Tuple{Float64,Float64,Float64,Vector{Float64},Float64}
```

返回 `(pet_s, pet_c, watdef, dθ, dθ_deep)`。

## 3. 算法 / 公式

**水分提取便利性**（src/extraction.jl:L82-L86，源自 Fortran 原型）：

$$
\text{easy}_k = \max\!\left(-\frac{(\mathrm{POTLEAF}-\psi_k)\,\mathrm{soilfactor}_k}{h_{\mathrm{veg}}-z_k},\;0\right)
$$

其中 $\psi_k = \psi_{\mathrm{sat},k}(\theta_{\mathrm{sat},k}/\theta_k)^b$，$\mathrm{soilfactor}_k = \mathbb{1}_{\text{icefac}_k=0}$。

**根区土壤水胁迫因子**（L141）：

$$
f_{\mathrm{swp}} = \begin{cases}0, & \text{rootfc}=0\\ \mathrm{clamp}(\theta_{\mathrm{root}}/\mathrm{rootfc},\,0,\,1), & \text{其他}\end{cases}
$$

**冠层气孔阻力**（L144）：

$$
r_{s,c} = \begin{cases}5000, & f_{\mathrm{swp}}=0\\ \min(r_{s,c}^{\mathrm{factor}}/f_{\mathrm{swp}},\;5000), & \text{其他}\end{cases}
$$

**Penman-Monteith 蒸腾**（L155-L159）：

$$
\mathrm{pet}_c = \max\!\left(\frac{C_c\,\mathrm{petfactor}_c}{(\Delta+\gamma\,\bigl(1+r_{s,c}/(r_{a,a}+r_{a,c})\bigr))}\,\frac{\Delta t}{\lambda},\;0\right)
$$

其中分配系数 $C_c = 1/(1+R_a R_c/(R_s(R_c+R_a)))$，$R_c=(\Delta+\gamma)r_{a,c}+\gamma r_{s,c}$。

## 4. 关键变量与单位

| 符号 | 含义 | 单位 |
|---|---|---|
| `POTLEAF` | 叶片水势阈值 | m（−153.0） |
| `POTWILT` | 凋萎点水势 | m（−153.0） |
| `POTFC` | 田间持水量水势 | m（−3.366） |
| `easy[k]` | 水分提取便利性 | 无量纲 |
| `ψ_k` | 土壤基质势 | m |
| `θ_sat`, `θ_fc`, `θ_min` | 饱和/田间/凋萎含水量 | m³/m³ |
| `fdepth` | 根系衰减深度因子 | m |
| `icefac[k]` | 冰冻抑制因子（0/1） | flag |
| `fswp` | 水分胁迫因子 | [0,1] |
| `rs_c`, `rs_s` | 冠层/土壤表面阻力 | s/m |
| `pet_c`, `pet_s` | 冠层蒸腾/土壤蒸发 | m |
| `watdef` | 未满足水分亏缺 | m |
| `inactivedays` | 连续非活跃天数 | d |

## 5. 与 Fortran 对应

| Fortran 子程序 | Julia 函数 | 差异 |
|---|---|---|
| `extraction`（fortran/extraction.f90） | `extraction`（src/extraction.jl） | Julia 版引入 `fdepth` 深度因子调整 $\theta_{\mathrm{sat}}$、$\psi_{\mathrm{sat}}$（L76-L78）；新增 `maxeasy×0.001` 阈值剪枝（L99）与非活跃层额外抑制（L101-L106）；kroot 返回 `k-1` 而非 `k`（注释 L71 提示待复核） |

## 6. 引用

- 行号：src/extraction.jl:L1-L4（常量）；L76-L86（便利性）；L122-L141（fswp）；L144（rs_c）；L155-L186（PET 与提取）
- 测试断言：test/test_extraction.jl:L37 `@test dθ_deep == 0.0`；L48 `@test pet_c == 0.0`（零 LAI）；L60 `@test all(dθ .== 0.0)`（全冰冻）；L82 `@test inactivedays`（非活跃日更新）；L98-L102 时间步长线性比例；L118 提取守恒（`sum(dθ) ≈ pet_c*1e-3`）
- 文档：docs/extraction.typ（待摄取）

## 7. 已知问题与备注

1. `kroot = k - 1`（L67-L72 注释）：源码注释质疑"为何返回 k-1 而非 k"，可能导致根区最浅层被跳过，需与 Fortran 行为对照。
2. 当前实现中根系不能直接从地下水层（`kroot < jwt`）抽水，注释 L73 明确指出。
3. `dθ_deep` 始终返回 0.0，地下水直接耦合尚未启用；测试 L46 显式断言此点。
4. `POTLEAF == POTWILT = -153.0` 为常量，未暴露为参数；如需敏感性分析应参数化。
5. `rs_s`（L148）依赖表层 `θ[nzg]`，而非全根区平均；浅层优先假设隐含其中。
6. 导出表：`export extraction`（L1）已注册；模块内无悬空 export。