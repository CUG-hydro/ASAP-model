# SoilFluxes — 土壤水运动

> 源文件：src/SoilFluxes.jl:L1-L401
> Fortran 来源：fortran/SoilWater.f90:SoilWater（迁移自 CLM4.5 / CoLM）
> 测试：test/test_soil_fluxes.jl
> 状态：已摄取

## 1. 功能概述

本模块实现一维 Richards 方程的 Crank-Nicolson 有限差分离散，求解土壤剖面（`nzg` 层，含地下水位 `jwt`）在时间步长 `dt` 内的水分通量与含水量更新。其核心任务包括：基于 Campbell 关系的导水率 `K(θ)` 与水力扩散率 `D(θ)` 计算、三对角系统的 Thomas 算法求解、顶部入渗能力约束、含水量超饱和 / 凋萎含水量修正，以及自由排水与受限（`freedrain`）两种底部边界条件。

## 2. 函数签名

```julia
function soilfluxes(
  nzg::Int, dt::Float64,
  z₋ₕ::Vector{Float64}, Δz::Vector{Float64}, soiltxt::Int,
  θ_wtd::Float64, transp::Vector{Float64}, transpdeep::Float64,
  θ::Vector{Float64}, wtd::Float64, precip::Float64, pet_s::Float64,
  fdepth::Float64, qlat::Float64, qrf::Float64, flood::Float64,
  icefactor::Vector{Int8}; freedrain::Bool=true
) :: Tuple

function tridag!(
  a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64},
  r::Vector{Float64}, u::Vector{Float64}
) :: Nothing
```

## 3. 算法 / 公式

### 3.1 Richards 方程与通量分解（来自 docs/土壤水运动_ASAP.typ:第 1-2 式）

$$
Q_{k+1/2} = -K_{k+1/2}(\theta) - D_{k+1/2}(\theta)\,\frac{\theta_{k+1}^{n+1} - \theta_k^{n+1}}{\Delta z_{k+1/2}}
$$

其中第一项为**重力通量**，第二项为**毛管通量**。

### 3.2 隐式离散三对角方程

合并后整理为 $a\theta_{k-1}^{n+1} + b\theta_k^{n+1} + c\theta_{k+1}^{n+1} = r$ 形式（docs/土壤水运动_ASAP.typ 中第 3 ~ N-1 层代码）：

```julia
aa[k] =  D_mid[k]   / Δz₊ₕ[k]
cc[k] =  D_mid[k+1] / Δz₊ₕ[k+1]
bb[k] = -(aa[k] + cc[k] + Δz[k] / dt)
rr[k] = -θ[k] * Δz[k] / dt - K_mid[k+1] + K_mid[k] + transp[k] / dt
```

### 3.3 Campbell 关系

$$
K(\theta) = K_\text{sat}\left(\frac{\theta}{\theta_\text{sat}}\right)^{2b+3},\quad
D(\theta) = -K_\text{sat}\,\psi_\text{sat}\,\frac{b}{\theta_\text{sat}}\left(\frac{\theta}{\theta_\text{sat}}\right)^{b+2}
$$

为模拟根系吸水对水力参数的影响，对 $K_\text{sat}$、$\theta_\text{sat}$、$\psi_\text{sat}$ 沿深度按 `fdepth` 作 `exp(z/fdepth)` 缩放（src/SoilFluxes.jl 中 `_Ksat`、`_θsat`、`_ψsat` 计算）。

### 3.4 入渗能力约束

$$
I_\text{max} = K_\text{sat}\,\Delta t,\quad \text{runoff} = \max(0,\,|Q_{N+1/2}| - I_\text{max})
$$

### 3.5 Green-Ampt 改进（来自 docs/下渗_黄土高原.typ:王文焰 2003）

对于黄土高原干土入渗，采用**椭圆面积修正**：

$$
\frac{\mathrm{d}V}{\mathrm{d}t} = \frac{16\,K(\theta_s)}{(\pi + 4)(\theta_s - \theta_l)}\,\frac{H + 0.5 L - S_m}{L}
$$

系数 $16/(\pi+4)$ 为湿润锋椭圆浸润面积的等效修正（原 Green-Ampt 取 $1/2$），$S_m$ 为负的基质势，$H$ 为积水深度，$L$ 为入渗锋面深度。

## 4. 关键变量与单位

| 符号                | 含义                                 | 单位           |
| ------------------- | ------------------------------------ | -------------- |
| `nzg`               | 土壤层数                             | -              |
| `jwt`               | 地下水位所在层号                     | -              |
| `dt`                | 时间步长                             | s              |
| `z₋ₕ`               | 层界面深度                           | m              |
| `Δz`                | 层厚度                               | m              |
| `θ`, `θ_wtd`        | 体积含水量 / 地下水位处含水量        | m³/m³          |
| `θ_sat`, `θ_cp`     | 饱和含水量 / 凋萎含水量              | m³/m³          |
| `Ksat`, `b`, `ψsat` | 饱和导水率 / Campbell b / 饱和基质势 | m s⁻¹ / - / m  |
| `K_mid`, `D_mid`    | 层间 K(θ) 与 D(θ)                    | m s⁻¹ / m² s⁻¹ |
| `Q`, `flux`         | 层间水分通量数组（向下为负）         | m s⁻¹          |
| `Imax`              | 入渗能力                             | m              |
| `runoff`            | 超渗产流量                           | m              |
| `rech`              | 地下水补给量（= Q[1]）               | m              |
| `transp`            | 根系吸水提取量                       | m              |
| `icefactor`         | 冰冻抑制因子 (0~1)                   | -              |
| `freedrain`         | 自由排水标志                         | Bool           |

## 5. 与 Fortran 对应

| Fortran 子程序     | Julia 函数             | 差异                                                                                                                                                  |
| ------------------ | ---------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------- |
| `SoilWater()`      | `soilfluxes(...)`      | Julia 用 kwargs（`freedrain::Bool`）替代 Fortran 整型标志；返回元组替代 COMMON 块；删除了氧18同位素过程（src/SoilFluxes.jl 末尾 `o18ratio` 段被注释） |
| `tridi()`          | `tridag!(...)`         | Julia 原地修改右端向量 `r` 与解 `u`，省略 Fortran 的 `bet` 工作数组（用 `gam` 替代）                                                                  |
| `theta_limits()`   | 内联于 `soilfluxes` 中 | Julia 在主循环末尾统一校正超饱和 / 过干                                                                                                               |
| `infil_capacity()` | `Imax = Ksat*dt` 内联  | 不再单独子程序化                                                                                                                                      |

## 6. 引用

- 行号：src/SoilFluxes.jl:L26-L42 顶部边界与入渗能力；L45-L55 K(θ)/D(θ) 计算；L82-L134 三对角组装；L180-L198 底部边界；L300-L328 含水量修正；L350-L401 tridag!
- 公式文档：docs/土壤水运动_ASAP.typ 第 1、2 式与 jwt 周围分层规则（第 1~`jwt-3` 层 Q=0；`jwt-2` 层只算下边界 Q；`jwt-1` 层全计算并约束方向）
- 测试断言：test/test_soil_fluxes.jl L40 自由排水底层 Q[1]<0；L100 非自由排水 `flux[1]==0.0`；L160 `flux_nonfree==0.0 && rech_nonfree==0.0`；tridag! 解析解测试以 `[2,1;1,2,1;1,2] x=[1,2,3]` 验证

## 7. 已知问题与备注

1. **docstring 清理（2026-06-21）**：函数 docstring 已移除氧 18 相关形参（`o18`、`precipo18`、`tempsfc`、`qlato18`、`transpo18`）的条目，返回值改为 `(et_s, runoff, rech, flux, qrfcorrect, updated_θ)`；`RootDepth.jl` 中 `# transpo18[i, j] += ...` 与 `# o18[:, i, j] .= updated_o18` 等注释段已清理。
2. **遗留死代码（部分清理）**：`o18ratio = o18 ./ θ` 等 9 行氧 18 同位素三对角求解代码仍保留在函数尾部（src/SoilFluxes.jl 末段）作为参考，待 `IsotopeTracing` 模块与主循环完全串联后恢复。
3. **`smoiwtd` 未回写**：注释 `# smoiwtd = smoiwtd - Q[1]` 表明自由排水时地下水位含水量的更新被推迟到调用方（避免双重计算）。
4. **`fdepth` 含义模糊**：Fortran 中 `fd` 为根系深度因子，本模块既用于 `Ksat/θsat/ψsat` 缩放，又作为 `pet_s_actual` 判断的隐含阈值（`θ_cp` 比较），需保持调用一致性。
5. **`jwt` 变量未在签名中**：`jwt` 通过 `wtd` 与 `z₋ₕ` 在函数内部反算（搜索 `for k in max(jwt-1,2)` 等），是耦合 Fortran 习惯的遗留。
6. **`transpdeep` 与 `qlat`、`qrf`、`flood`**：形参已声明但当前实现未在可见代码中消费，需查阅 `git log` 确认是否在最近重构中被搁置。
7. **冰冻因子 `f_ice`**：src/SoilFluxes.jl L52 处出现裸引用 `f_ice`，应来自父模块作用域（`icefactor` 经 `prod(icefactor)` 或类似映射传入），建议在最终版统一为 `icefactor[k]`。
