# soilfluxes.f90

> 源文件：`fortran/soilfluxes.f90`
> 行数：609
> Module 声明：(无 module 包裹，独立子程序文件)
> Julia 来源：`src/SoilFluxes.jl::soilfluxes` + `src/updatewtd_shallow.jl::updatewtd_shallow`
> 测试：Fortran 侧无独立单元测试；逻辑等价于 `test/test_soil_fluxes.jl`
> 状态：已摄取

## 1. 功能概述

本文件实现 ASAP 模型的 **1D Richards 方程求解器** `SOILFLUXES`：对单一像元 `(i, j)` 的 `nzg` 层土壤柱，在一个时间步 `dtll` 内，按 Crank-Nicolson 隐式格式组装三对角系统，调用 `tridag` 求解新含水量 `smoi`，并分别处理顶部（入渗能力截断 + 土壤蒸发）与底部（自由排水 / 受限无流）边界条件。求解后依次进行：含水量超饱和 → 凋萎含水量的双向限幅、底层重力排水累积到 `rech`、调用 `UPDATESHALLOWWTD` 更新地下水位 `wtd`、同位素 ¹⁸O 三对角求解、侧向地下水流 `qlat` 在水位以上的逐层分配。整个流程是 ASAP 物理过程的核心能量/质量守恒环节。

## 2. 函数签名

```fortran
subroutine SOILFLUXES(i, j, nzg, freedrain, dtll, slz, dz, soiltxt, smoiwtd, &
                      transp, transpdeep, smoi, wtd, rech, deeprech, precip, &
                      pet_s, et_s, runoff, flux, fdepth, qlat, qlatflux, qrf, &
                      qrfcorrect, flood, icefactor, smoieq, o18, precipo18, &
                      tempsfc, qlato18, transpo18)
```

形参类型与含义：

| 形参 | 类型 | 含义 | 单位 |
|---|---|---|---|
| `i, j` | `integer` | 像元坐标（用于诊断输出） | - |
| `nzg` | `integer` | 土壤层数 | - |
| `freedrain` | `integer` | 自由排水标志 (0=受限, 1=自由) | - |
| `dtll` | `real` | 主时间步长 | s |
| `slz` | `real(nzg+1)` | 层中心深度（向下为负） | m |
| `dz` | `real(nzg)` | 层厚度 | m |
| `soiltxt` | `integer(2)` | 土壤类型（1=深/2=浅） | - |
| `smoiwtd` | `real` | 地下水位处含水量 | m³/m³ |
| `transp` | `real(nzg)` | 根系吸水（mm 步内累计） | m |
| `transpdeep` | `real` | 深层根系吸水（未消费） | m |
| `smoi` | `real(nzg)` | 当前含水量（输入+输出） | m³/m³ |
| `wtd` | `real` | 地下水位（负值） | m |
| `rech` | `real` | 重力排水（输出） | m |
| `deeprech` | `real` | 深层 rech（未消费） | m |
| `precip` | `real` | 降水 | mm |
| `pet_s, et_s` | `real` | 土壤 PET / 实际 ET（输出） | mm |
| `runoff` | `real` | 累计产流（输出） | m |
| `flux` | `real(nzg+1)` | 累计层间通量（m），向下为负 | m |
| `fdepth` | `real` | 根系深度因子 | m |
| `qlat` | `real` | 侧向流总通量（m 步内） | m |
| `qlatflux` | `real(0:nzg+1)` | 各层侧向流分配（输出） | m |
| `qrf` | `real` | 地下水补给河道 | m |
| `qrfcorrect` | `real` | qrf 修正量（输出） | m |
| `flood` | `real` | 漫流层输入 | m |
| `icefactor` | `integer*1(nzg)` | 冰冻因子 (0=冻, 1=液态) | - |
| `smoieq` | `real(nzg)` | 平衡含水量 | m³/m³ |
| `o18` | `real(nzg)` | ¹⁸O 浓度（输入+输出） | ‰·m³/m³ |
| `precipo18` | `real` | 降水 ¹⁸O 比 | ‰ |
| `tempsfc` | `real` | 地表温度（用于平衡分馏） | K |
| `qlato18` | `real` | 侧向流 ¹⁸O 总通量 | ‰·m |
| `transpo18` | `real` | 累计蒸腾 ¹⁸O 通量（输出） | ‰·m |

## 3. 算法 / 公式

### 3.1 几何辅助向量

```fortran
vctr2(k) = 1.0 / dz(k)                      ! 层厚倒数
vctr4(k) = 0.5 * (slz(k) + slz(k+1))        ! 层中心 (k 与 k+1 节点中点)
vctr5(k) = vctr4(k) - vctr4(k-1)            ! 中心间距
vctr6(k) = 1.0 / vctr5(k)                   ! 间距倒数（仅 k≥2）
```

`vctr4` 与 `vctr6` 用于层间通量的距离归一化。

### 3.2 土壤水力参数 Campbell 关系

按深度分支取土类（`slz(k) < -0.30` 用深土 `soiltxt(1)`，否则用 `soiltxt(2)`），再以 `fdepth` 做垂向衰减：

$$
K_\text{sat,eff} = \text{slcons}(n_\text{soil}) \cdot \mathrm{clip}\!\left(\exp\!\left(\frac{z+1.5}{f_\text{depth}}\right),\, [0.1,\, 1.0]\right)
$$

$$
\theta_\text{sat,eff} = \text{slmsts}(n_\text{soil}) \cdot \mathrm{clip}\!\left(\exp\!\left(\frac{z+1.5}{f_\text{depth}}\right),\, [0.1,\, 1.0]\right)
$$

$$
\psi_\text{sat,eff} = \text{slpots}(n_\text{soil}) \cdot \mathrm{clip}\!\left(\exp\!\left(-\frac{z+1.5}{f_\text{depth}}\right),\, [1.0,\, 10.0]\right)
$$

`wgpmid` 是基于线性外推的层间含水量，被 `min(., θ_sat,eff)` 截断。

水力传导与扩散：

$$
K_\text{mid}(k) = i_\text{ice} \cdot K_\text{sat,eff} \cdot \left(\frac{w_{gp,\text{mid}}}{\theta_\text{sat,eff}}\right)^{2b+3}
$$

$$
D_\text{mid}(k) = -i_\text{ice} \cdot \frac{K_\text{sat,eff} \cdot \psi_\text{sat,eff} \cdot b}{\theta_\text{sat,eff}} \cdot \left(\frac{w_{gp,\text{mid}}}{\theta_\text{sat,eff}}\right)^{b+2}
$$

其中 `icefac = (icefactor(k)==0) ? 0 : 1`，即冰冻层完全抑制水力通量。

### 3.3 三对角系统组装（Crank-Nicolson）

对 $k = \max(\text{iwtd}, 3) \mathrel{..} n_\text{zg}$：

$$
a_k = D_\text{mid}(k) \cdot v_6(k),\quad c_k = D_\text{mid}(k+1) \cdot v_6(k+1)
$$

$$
b_k = -(a_k + c_k + \Delta z(k)/\Delta t)
$$

$$
r_k = -\theta_k^n \Delta z(k)/\Delta t - K_\text{mid}(k+1) + K_\text{mid}(k) + \text{transp}(k)/\Delta t
$$

其中 `iwtd` 是 `wtd` 所在的层下标（`freedrain==0` 时计算；`freedrain==1` 时 `iwtd=0`，全部层参与计算）。

### 3.4 顶部边界条件

$$
Q_{n_\text{zg}+1/2} = (-\text{precip} + \text{pet}_s)\cdot 10^{-3} - \text{flood} \quad [\text{m}]
$$

入渗能力截断（避免超过 `K_sat * dt`）：

$$
\text{runoff} = \max(0,\, -Q_{n+1/2} - \text{slcons}(n_\text{soil}) \cdot \Delta t)
$$

$$
Q_{n+1/2} = -\min\!\left(\text{slcons}(n_\text{soil}) \cdot \Delta t,\; |Q_{n+1/2}|\right)
$$

随后装配 `rr(nzg)`、`bb(nzg)`，并按 `iwtd-1 == nzg` 与否分支处理（极端情形水位在地表时 `aa=cc=0`）。

### 3.5 底部边界条件

**`freedrain == 0`（受限排水）**：

按水位深度分 3 段子情形：

| `iwtd` 区间 | 处理 |
|---|---|
| `iwtd ≤ 2`（水位极深） | `aa(1)=0`，`cc(1)=D(2)v_6(2)`，`bb(1)=-(cc(1)+dz(1)/dt)`，并显式组装 `k=2` |
| `1 < iwtd ≤ 3` | `k=1..iwtd-3` 设 `aa=cc=0, bb=1, rr=smoi(k)`（强制不变），`k=iwtd-1` 设 `aa=0, cc=D(k+1)v_6(k+1)`，加 `min(kfmid(k) + D(k)v_6(k)(smoi(k)-smoi(k-1)), 0)` 防止向上渗漏；`k=iwtd-2` 设 `aa=cc=0, bb=-dz(k)/dt, rr=...+max(-kfmid(k+1) - D(k+1)v_6(k+1)(smoi(k+1)-smoi(k)), 0)` |
| `iwtd > 3` | 水位以上层冻结 (`aa=cc=0, bb=1, rr=smoi`)；水位层 `iwtd-1` 仅计算下边界；`iwtd-2` 仅做下方通量 |

**`freedrain == 1`（自由排水）**：

$$
K_\text{mid}(1) = K_\text{sat,eff} \cdot (\theta_1/\theta_\text{sat,eff})^{2b+3}
$$

$$
aa(1)=0,\quad cc(1)=D(2)v_6(2),\quad bb(1)=-(cc(1) + \Delta z(1)/\Delta t)
$$

$$
rr(1) = -\theta_1 \Delta z(1)/\Delta t - K_\text{mid}(2) + K_\text{mid}(1) + \text{transp}(1)/\Delta t
$$

### 3.6 求解与通量重算

调用 `tridag(aa, bb, cc, rr, smoi, nzg)` 更新 `smoi`。随后按层间通量分解：

$$
\text{gravflux}(k) = -K_\text{mid}(k) \Delta t
$$

$$
\text{capflux}(k) = -a_k(\theta_k - \theta_{k-1}) \Delta t
$$

$$
\text{vt3di}(k) = \text{capflux}(k) + \text{gravflux}(k)
$$

水位层 `iwtd-1` 做方向约束：若 `capflux > -gravflux` 则保留；若 `capflux > 0` 则将 `gravflux = -capflux`，`vt3di = 0`；否则两侧清零。

新含水量：

$$
\theta_k^\text{old} \mathrel{+}= \frac{\text{vt3di}(k) - \text{vt3di}(k+1) - \text{transp}(k)}{\Delta z(k)}
$$

### 3.7 含水量限幅

**超饱和（自上而下）**：

$$
\Delta\theta = \max((\theta_k^\text{old} - \theta_\text{sat,eff}) \cdot \Delta z(k),\; 0)
$$

若 $k < n_\text{zg}$，则将超额水分推到下层 `smoiold(k+1)`，并把 `vt3di(k+1) += Δθ`；若 $k = n_\text{zg}$，计入 `runoff`；并相应把 `gravflux(k+1)` 调整为 `max(0, gravflux + capflux)`。

**过干（自下而上）**：

对 $k = n_\text{zg} \mathrel{..} 1$，若 $\theta_k^\text{old} < \theta_\text{cp,eff}$：

$$
\Delta\theta = \max((\theta_\text{cp,eff} - \theta_k^\text{old}) \cdot \Delta z(k),\; 0)
$$

顶部层优先扣减 PET（`et_s = max(0, pet_s - Δθ*1e3)`）；下方层从 `smoiold(k-1)` 借水：`smoiold(k-1) -= Δθ/dz(k-1)`，`vt3di(k) += Δθ`，并把 `smoiold(k) = θ_cp,eff`。

### 3.8 自由排水补给累积

```fortran
if (freedrain == 1) then
   rech = vt3di(1)                       ! [m] 底层重力排水
   smoiwtd = smoiwtd - vt3di(1)          ! 水位桶扣减
end if
```

### 3.9 ¹⁸O 同位素求解（独立三对角）

`o18ratio(k) = o18(k) / smoiold(k)`，`o18dz = o18`（备份）。

对每个通量 `vt3di(k)` 与 `capflux`、`gravflux`，先做"方向修正"：

- 若 `capflux(k) < 0` → `gravflux(k) += capflux(k); capflux(k) = 0`
- 计算 `fluxdiff = vt3di(k) - (gravflux(k)+capflux(k))`，正部分并入 `capflux`，负部分并入 `gravflux`

随后按 Crank-Nicolson 风格的同位素方程组装：

$$
b_k = 1 + \frac{\Delta z^{-1}(k)}{2}\!\left[\frac{\text{transp}(k)}{\theta_k^\text{old}} - \frac{\text{gravflux}(k)}{\theta_k^\text{old}} + \frac{\text{capflux}(k+1)}{\theta_k^\text{old}}\right]
$$

$$
r_k = o_{18,k} - \frac{\Delta z^{-1}(k)}{2}\!\left[\text{transp}(k)\,o_{18r,k} - \text{gravflux}(k)\,o_{18r,k} - \text{capflux}(k)\,o_{18r,k-1} - \text{gravflux}(k+1)\,o_{18r,k+1} - \text{capflux}(k+1)\,o_{18r,k}\right]
$$

$$
a_k = -\frac{\Delta z^{-1}(k)}{2}\frac{\text{capflux}(k)}{\theta_{k-1}^\text{old}},\quad c_k = \frac{\Delta z^{-1}(k)}{2}\frac{\text{gravflux}(k+1)}{\theta_{k+1}^\text{old}}
$$

**顶部边界**（含 Majoube 平衡分馏）：

$$
\alpha = \frac{1}{\exp(1137/T_\text{sfc}^2 - 0.4156/T_\text{sfc} - 0.0020667)}
$$

$$
b_{n_\text{zg}} \mathrel{+}= \frac{\Delta z^{-1}(n_\text{zg})}{2} \cdot \frac{\alpha \cdot \text{et}_s \cdot 10^{-3}}{\theta_{n_\text{zg}}^\text{old}}
$$

$$
r_{n_\text{zg}} \mathrel{-}= \frac{\Delta z^{-1}(n_\text{zg})}{2} \cdot (\alpha \cdot \text{et}_s \cdot 10^{-3}) \cdot o_{18r,n_\text{zg}}
$$

并按 `vt3di(n+1) - et_s*1e-3` 正负分支，决定降水同位素 `precipo18` 与表层水分混合的符号。

求解后 `transpo18 += 0.5*(o18(k)/smoiold(k) + o18ratio(k))*transp(k)`，并对负值或 `o18 > smoiold` 触发 `write(6,*)` 诊断。

### 3.10 地下水位更新与侧向流分配

调用 `UPDATESHALLOWWTD(i, j, nzg, freedrain, slz, dz, soiltxt, smoieq, smoiwtd, smoiold, wtd, rech, fdepth)` 更新 `wtd`，再按 `wtd` 重算 `iwtd`，并令 `kwtd = max(iwtd-1, 1)`。

**侧向流分配（`qgw = qlat - qrf`）**：

- 若 `qgw > 0`：自 `kwtd-1` 起向上逐层填充至 `θ_sat,eff`：
  - 若 `qlatlayer ≤ (θ_sat - smoiold)·dz`：一次性加入该层，`o18(k) += (qlatlayer/qgw) * qgwo18 / dz`
  - 若处于 `kwtd-1` 层：将余量累计并继续向上
  - 否则加入该层后强制截断到 `θ_sat,eff`，剩余 `qlatlayer` 与 `qgwo18` 继续上传；到达 `nzg` 时计入 `runoff`

- 若 `qgw < 0`：自 `kwtd` 起向下逐层抽取，约束到 `smoieq(k)`（最深层用 `soilcp`），当 ¹⁸O 不足时同时从下一层 `o18(k-1)` 按比例扣减；`k=1` 时记入 `qrfcorrect -= qlatlayer`

最终 `o18 = max(o18, 0)`，`smoi = smoiold`。

## 4. 关键变量与单位

| 符号 | 含义 | 单位 |
|---|---|---|
| `vctr2, vctr4, vctr5, vctr6` | 几何辅助向量（厚度倒数/层中心/间距/间距倒数） | 1/m 或 m |
| `kfmid, diffmid` | 层间 `K(θ)` 与 `D(θ)` | m/s 与 m²/s |
| `aa, bb, cc, rr` | 三对角系数与右端向量 | m/s, m³/m³ |
| `vt3di, gravflux, capflux` | 节点通量（重力 + 毛细分量） | m/s 或 m |
| `flux` | 累计层间通量（步内增量被累加） | m |
| `iwtd, kwtd` | `wtd` 所在层 / 计算层 | - |
| `alpha` | Majoube 平衡分馏系数 | - |
| `dsmoi` | 超饱和 / 过干差额（乘 dz 后为水量） | m³/m² |
| `qgw, qgwo18, qlatlayer, o18frac` | 侧向流与同位素分配中间量 | m, ‰, -, - |
| `tridag` | Thomas 三对角算法（外部子程序） | - |

## 5. 与 Julia 对应

| Fortran 子程序 / 段 | Julia 函数 / 段 | 差异 |
|---|---|---|
| `SOILFLUXES` 主循环 | `soilfluxes(...)`（`src/SoilFluxes.jl:35-L376`） | Julia 用 kwargs (`freedrain::Bool`) 替代 Fortran 整型标志；返回元组替代 COMMON 块；删除了氧18同位素过程（末段被注释） |
| 三对角求解 | `tridag!`（`src/SoilFluxes.jl:377-L401`） | Julia 原地修改右端向量 `r` 与解 `u`；省略 Fortran 的 `bet` 工作数组 |
| Campbell `K(θ)`/`D(θ)` 计算 | `cal_K` + 内联 D 公式 | Julia 显式 `Float64` 类型；Fortran 通过 `(smoi, nsoil)` |
| `fdepth` 深度衰减 | 同上 | 公式一致 |
| 顶部入渗能力截断 | 内联 `Imax = Ksat*dt` | 不再单独子程序化 |
| 底部自由排水 vs 受限 | `if !freedrain` / `else` 分支 | 接口对齐 |
| 含水量限幅（超饱和 / 过干） | `soilfluxes` 末段循环 | Julia 在主循环末尾统一校正；Fortran 同段中分两段独立 `do` 循环 |
| 同位素段（`alpha`、`o18ratio`、`tridag` for o18） | **当前未启用**（`src/SoilFluxes.jl` 末段被注释；`wiki/julia/IsotopeTracing-同位素追踪.md` 是另一路径） | Julia 中 `RootDepth.jl` `# transpo18[i, j] += ...` 等注释段已清理；Fortran 仍为活动代码 |
| `UPDATESHALLOWWTD` | `updatewtd_shallow`（`src/updatewtd_shallow.jl:26-L102`） | 接口对齐；Julia 内联 `find_jwt`；输出改为多返回值 |
| 侧向流分配 `qgw` 逐层填充 | 内联于 `SoilFluxes.jl` 末段（`qgw > 0` 分支） | Julia 同样按 `kwtd-1` 起向上填至 θ_sat,eff；`qgw < 0` 分支已被注释 |

## 6. 引用

- 行号：
  - `fortran/soilfluxes.f90:L1-L29` 形参与几何辅助
  - `fortran/soilfluxes.f90:L33-L36` K(θ)/D(θ) 公式
  - `fortran/soilfluxes.f90:L18-L25` 顶部边界 + 入渗能力截断
  - `fortran/soilfluxes.f90:L63-L100` 受限排水（iwtd 分支）
  - `fortran/soilfluxes.f90:L101-L116` 自由排水底部
  - `fortran/soilfluxes.f90:L120-L180` 三对角求解后通量重算
  - `fortran/soilfluxes.f90:L182-L235` 超饱和校正
  - `fortran/soilfluxes.f90:L237-L271` 过干校正（含 PET 截断）
  - `fortran/soilfluxes.f90:L292-L297` 重力排水累积到 rech
  - `fortran/soilfluxes.f90:L302-L394` ¹⁸O 三对角系统组装与求解
  - `fortran/soilfluxes.f90:L398-L488` 侧向流分配 + ¹⁸O 同步
- 关联：`module_rootdepth.f90:994 SOILFLUXES`（实际入口），`module_wtable.f90:1226 FLOODING`（后续调用）
- 文档：`docs/土壤水运动_ASAP.typ`（对应 Julia 侧；Fortran 公式与之一致）
- 公式来源：Campbell (1974) 导水率关系；Majoube (1971) 平衡分馏公式
- 对照：`wiki/julia/SoilFluxes-土壤水运动.md`、`wiki/mapping/julia-fortran-对照.md` §2（Richards 章节）

## 7. 已知问题与备注

1. **无独立 module 声明**：与 `main.f90` 类似，`SOILFLUXES` 没有 `MODULE` 包裹，依赖调用方提供 `slcons`、`slmsts`、`slpots`、`slbs`、`soilcp`、`tridag`、`UPDATESHALLOWWTD` 等模块符号；编译时需在调用上下文中 `USE module_rootdepth`。
2. **`transpdeep`、`deeprech`、`qlatflux(0)` 与 `qlatflux(nzg+1)` 等形参**：在当前 Fortran 实现中已被声明，但部分边界场景（`qgw > 0` 与 `qgw < 0` 分支）未消费全部元素（`qgw < 0` 时 `qlatflux(0)` 未累加）。Julia 侧对应形参亦未完全消费，已在 `wiki/julia/SoilFluxes-土壤水运动.md §7` 登记。
3. **`freedrain` 整型标志**：`0`/`1` 双值；Julia 改为 `Bool`，避免整数歧义。
4. **同位素段诊断输出**：`write(6,*) 'O18 less than zero!!!'` 等语句直接写 stdout 6 号单元，并触发于任何 `o18 < 0` 或 `o18 > smoiold` 的像元；运行时会在 `rootdepth_main` 频繁触发，仅做调试用。
5. **`pppendepth` 形参已注释**：早期降水穿透深度计算（`if(vt3di(k).lt.-1.e-6)then ... endif` 段）整段被注释，对应 Julia 中已删除的 `pppendepthold` 维护逻辑。
6. **`alpha` 平衡分馏仅在地表使用**：深层 `o18` 通量计算未应用平衡分馏系数，符合水-汽交换主要发生在表层的物理图像。
7. **`qgw < 0` 时最深层 `k=1` 的 `o18` 截断**：当 `o18(k)*dz(k) < o18out` 且 `k>1` 时按 `o18(k)*dz(k) / (o18(k)*dz(k) + o18(k-1)*dz(k-1))` 比例从相邻层借水；这是简单近似，未严格守恒。
8. **`qsprings`、`qslat`、`transptop` 形参未列入**：本子程序签名不包含这些变量；它们在 `module_rootdepth.f90::ROOTDEPTH` 主调度中处理。
9. **`icefactor` 与 `icefac`**：`icefactor(k)` 为 `integer*1`（0/1），`icefac` 为 `real`（0./1.）；当 `icefactor(k)==0` 时 `icefac=0`，完全抑制 K 与 D；非冻结态 `icefac=1`。