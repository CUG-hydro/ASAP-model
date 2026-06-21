# ASAP-model 物理过程 → 实现函数 正向索引

> 用途：以**物理过程**为主键，给出一行数学描述、Julia 实现、Fortran 实现、测试位置与文档参考。
> 状态：基于 `wiki/julia/*.md`、`wiki/fortran/*.md`、`src/` 与 `fortran/` 源码直接核对，覆盖至少 18 类物理过程。
> 路径约定：所有路径均为绝对路径；行号以源码实测为准。

## 1. 土壤水力参数化（13 类 USDA 分类）

| 物理过程 | 数学描述 | Julia 实现 | Fortran 实现 | 测试位置 | 文档参考 |
|---|---|---|---|---|---|
| 13 类土壤饱和含水量 `θ_sat` | 查表：`θSAT[1..13]` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/SoilParameters.jl:25` `const θSAT` | `fortran/module_rootdepth.f90` 模块级 `slmsts(13)` | `test/test_soil_parameters.jl:L6-L22` | `docs/土壤水运动_ASAP.typ` §土壤参数 |
| 13 类凋萎含水量 `θ_cp` | 查表 `SOILCP[1..13]` | `src/SoilParameters.jl:26` `const SOILCP` | `module_rootdepth.f90` `soilcp(13)` | `test/test_soil_parameters.jl:L6-L22` | 同上 |
| Clapp-Hornberger `b` 指数 | 查表 `SLBS[1..13]` | `src/SoilParameters.jl:27` `const SLBS` | `module_rootdepth.f90` `slbs(13)` | 同上 | 同上 |
| 饱和导水率 `K_sat` | 查表 `KSAT[1..13]` (m/s) | `src/SoilParameters.jl:28` `const KSAT` | `module_rootdepth.f90` `slcons(13)` | 同上 | 同上 |
| 饱和基质势 `ψ_sat` | 查表 `ΨSAT[1..13]` (m) | `src/SoilParameters.jl:30` `const ΨSAT` | `module_rootdepth.f90` `slpots(13)` | 同上 | 同上 |
| 侧向流比例因子 `K_lat/K_sat` | 查表 `KLATFACTOR[1..13]` | `src/SoilParameters.jl:31` `const KLATFACTOR` | `module_rootdepth.f90` `klatfactor(13)` | 同上 | 同上 |
| `SoilType` 类型化封装 | `struct SoilType` | `src/SoilParameters.jl:13-L21` | 无对应（Fortran 用裸数组） | 同上 | — |
| `get_soil_params(soil_type)` | 返回 `SoilType` | `src/SoilParameters.jl:43-L57` | 直接数组访问 `slmsts(nsoil)` | 同上 | — |
| `init_soil_param(nzg)` | 凋萎点初始化 | `src/SoilParameters.jl:68-L79` | `module_rootdepth.f90:1891 init_soil_param` | `test/test_soil_parameters.jl:L26` | — |
| Campbell 导水率 | $K = K_{sat}(\theta/\theta_{sat})^{2b+3}$ | `src/SoilParameters.jl:94-L96` `cal_K` | `module_rootdepth.f90:1906 khyd` | `test/test_soil_parameters.jl:L51` | `docs/土壤水运动_ASAP.typ` §Campbell |
| 深度依赖衰减 | `slmsts(z) = slmsts·clip(exp((z+1.5)/fdepth), 0.1, 1)` | `src/SoilFluxes.jl` 内联（`_Ksat/_θsat/_ψsat`） | `module_rootdepth.f90` 同公式 | `test/test_soil_fluxes.jl` | `docs/土壤水运动_ASAP.typ` §fdepth |

## 2. 土壤层初始化（CLM / 固定厚度）

| 物理过程 | 数学描述 | Julia 实现 | Fortran 实现 | 测试位置 | 文档参考 |
|---|---|---|---|---|---|
| CLM 指数节点深度 | $z_k = 0.025(e^{0.5(k-0.5)}-1)$ | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/SoilInitialization.jl:26-L52` `initializesoildepth_clm` | `fortran/module_rootdepth.f90:937 INITIALIZESOILDEPTHCLM` | `test/test_soil_initialization.jl:L6` | `docs/土壤水运动_ASAP.typ` §分层 |
| 固定 40 层厚度（`DZ2`） | 查表 `DZ2[end-n+1..end]` 倒序 | `src/SoilInitialization.jl:65-L80` `initializesoildepth` | `fortran/module_rootdepth.f90:973 INITIALIZESOILDEPTH` | `test/test_soil_initialization.jl:L14` | 同上 |
| 平衡含水量（理论法） | 每层 Brent 求根使稳态通量为零 | **无对应** | `fortran/module_initial.f90:141 EQSOILMOISTUREtheor` | — | `module_initial.f90` |
| 平衡含水量（迭代法） | Newton 迭代至稳态 | **无对应** | `fortran/module_initial.f90:321 EQSOILMOISTUREiter` | — | 同上 |
| Brent 求根 | 反二次插值 + 二分法 | **无对应** | `fortran/module_initial.f90:638 zbrent` + `func1/2/3` | — | `module_initial.f90` |
| `find_jwt(wtd, z₋ₕ)` | 二分查找水位所在层 | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/helper.jl:4-L29` | 内联于 `SOILFLUXES`、`UPDATESHALLOWWTD` 等多处 | `test/test_updatewtd_shallow.jl:L6` | — |

## 3. 三种潜在蒸散（PET）

| 物理过程 | 数学描述 | Julia 实现 | Fortran 实现 | 测试位置 | 文档参考 |
|---|---|---|---|---|---|
| Priestley-Taylor (PT) | $\lambda ET = \alpha \cdot \Delta R_n / (\Delta+\gamma)$ | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/Evapotranspiration.jl:70-L101` `potevap_priestly_taylor` | `fortran/module_rootdepth.f90:546 POTEVAP_Priestly_Taylor` | `test/test_evapotranspiration.jl:L9` | `docs/土壤水运动_ASAP.typ` §PET |
| Penman-Monteith (P-M, FAO) | $\lambda ET = \frac{\Delta(R_n-G)+\rho_a c_p VPD/r_a}{\Delta+\gamma(1+r_c/r_a)}$ | `src/Evapotranspiration.jl:103-L166` `potevap_penman_monteith` | `fortran/module_rootdepth.f90:569 POTEVAP_Penman_Monteith` | `test/test_evapotranspiration.jl:L25` | 同上 |
| Shuttleworth-Wallace (S-W, 双源) | $ET = C_c PM_c + C_s PM_s$；返回 12 元组 | `src/Evapotranspiration.jl:168-L288` `potevap_shutteworth_wallace` | `fortran/module_rootdepth.f90:653 POTEVAP_Shutteworth_Wallace` | `test/test_evapotranspiration.jl:L48-L60` | `wiki/julia/Evapotranspiration-蒸散发.md` |
| 31 类 USGS 生物物理参数 | `[z0m, 叶宽]` 表 | `src/Evapotranspiration.jl:15-L43` `BIOPARMS/RL/ZDIS` | `module_rootdepth.f90` 同类数组 | `test/test_evapotranspiration.jl` | — |
| 饱和水汽压斜率 `Δ` | Clausius-Clapeyron 数值导数 | `src/Evapotranspiration.jl` 内联 | `module_rootdepth.f90` 内联 | — | — |
| 湿度计常数 `γ` | `γ = c_p·P / (0.622·λ)` | `src/Evapotranspiration.jl:11-L13` | `module_rootdepth.f90` 同 | — | — |

## 4. 截留与截留蒸发

| 物理过程 | 数学描述 | Julia 实现 | Fortran 实现 | 测试位置 | 文档参考 |
|---|---|---|---|---|---|
| 截留容量 | $S_{max} = 0.2 \cdot LAI$ | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/Interception.jl:15` `intercepmax` | `fortran/module_rootdepth.f90:277 INTERCEPTION` 内联 | `test/test_interception.jl:L42` | `docs/土壤水运动_ASAP.typ` §截留 |
| 穿透降水 + 截留蒸发 | 水量平衡：$P_d = P - D$；$E_i = \min(S_{max}, PET_i)$ | `src/Interception.jl:18-L54` `interception` | `fortran/module_rootdepth.f90:277 INTERCEPTION` | `test/test_interception.jl:L20` | 同上 |
| 最小降水率阈值 | `precip < minpprate` 触发蒸发 | `src/Interception.jl:3` 常量（**接口未读**） | `module_rootdepth.f90:277` 全局参数 | — | — |
| 截留蒸发（SW 分量） | `pet_i`（取气孔阻力为 0 的极限） | `src/Evapotranspiration.jl:168-L288` `potevap_shutteworth_wallace` | `module_rootdepth.f90:653 POTEVAP_Shutteworth_Wallace` | 同 §3 | — |

## 5. 根系吸水与水分胁迫

| 物理过程 | 数学描述 | Julia 实现 | Fortran 实现 | 测试位置 | 文档参考 |
|---|---|---|---|---|---|
| 水分提取便利性 | $\text{easy}_k = \max\!\left(-\frac{(\text{POTLEAF}-\psi_k)\,\text{soilfactor}_k}{h_{veg}-z_k},\;0\right)$ | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/extraction.jl:82-L86` | `fortran/module_rootdepth.f90:307 EXTRACTION` | `test/test_extraction.jl:L37` | `docs/extraction.typ` |
| 土壤基质势 `ψ` | $\psi = \psi_{sat}(\theta_{sat}/\theta)^b$ | `src/extraction.jl:80` | `module_rootdepth.f90:307` 内联 | — | — |
| 水分胁迫因子 `f_swp` | `clamp(θ_root/rootfc, 0, 1)` | `src/extraction.jl:122-L141` | `module_rootdepth.f90:307` 内联 | `test/test_extraction.jl:L60` | — |
| 冠层气孔阻力 `rs_c` | $\min(r_{s,c}^{factor}/f_{swp}, 5000)$ s/m | `src/extraction.jl:144` | `module_rootdepth.f90:307` 内联 | — | — |
| Feddes 减少函数 | 整合 `fswp` × `easy` | `src/extraction.jl:155-L186` | `module_rootdepth.f90:307 EXTRACTION` | `test/test_extraction.jl:L118` | — |
| 非活跃层计数 | `inactivedays` 累积抑制 | `src/extraction.jl:101-L106` | `module_rootdepth.f90` 共享 `inactivedays` | `test/test_extraction.jl:L82` | — |
| 主算法时步入口 | 串联 PET → 截留 → 提取 → 土壤水 → 水位 | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/RootDepth.jl:90-L278` `rootdepth_main` | `fortran/module_rootdepth.f90:32 ROOTDEPTH` | `test/test_rootdepth.jl` | `docs/土壤水运动_ASAP.typ` §RootDepth |
| 冰冻层抑制 | `icefactor[k]=0` 跳过 | `src/SoilFluxes.jl:52`（裸引用 `f_ice`） | `module_rootdepth.f90:994 SOILFLUXES` 内部 | `test/test_soil_fluxes.jl` | — |

## 6. Richards 方程与 Thomas 算法

| 物理过程 | 数学描述 | Julia 实现 | Fortran 实现 | 测试位置 | 文档参考 |
|---|---|---|---|---|---|
| 1D Richards 方程 | $\partial\theta/\partial t = \partial/\partial z[D(\theta)\partial\theta/\partial z] + \partial K/\partial z$ | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/SoilFluxes.jl:35-L376` `soilfluxes` | `fortran/soilfluxes.f90:2 SOILFLUXES` / `fortran/module_rootdepth.f90:994` | `test/test_soil_fluxes.jl` | `docs/土壤水运动_ASAP.typ` §Richards |
| Crank-Nicolson 离散 | 隐式时间积分，`(1-α)·F^n + α·F^{n+1}` | `src/SoilFluxes.jl:82-L134` | `fortran/soilfluxes.f90` 中段 | 同上 | 同上 |
| 三对角组装 | $a\theta_{k-1}^{n+1} + b\theta_k^{n+1} + c\theta_{k+1}^{n+1} = r$ | `src/SoilFluxes.jl:82-L134` | `fortran/soilfluxes.f90` | 同上 | 同上 |
| Thomas 算法 | 顺消 + 回代 | `src/SoilFluxes.jl:377-L401` `tridag!` | `fortran/module_rootdepth.f90:1868 tridag` | `test/test_soil_fluxes.jl:L120+`（解析解 `[2,1;1,2,1;1,2]x=[1,2,3]`） | — |
| 含水量限幅 | `clamp(θ, θ_cp, θ_sat)` | `src/SoilFluxes.jl:300-L328` 内联 | `fortran/soilfluxes.f90` 内联 | `test/test_soil_fluxes.jl` | — |
| 层间通量 `Q` | 达西定律 + 离散 | `src/SoilFluxes.jl:45-L55` | `fortran/soilfluxes.f90` | 同上 | — |

## 7. 顶部 / 底部边界

| 物理过程 | 数学描述 | Julia 实现 | Fortran 实现 | 测试位置 | 文档参考 |
|---|---|---|---|---|---|
| 自由排水底部边界 | `flux[1] = 0`（水位处通量为 0） | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/SoilFluxes.jl:180-L198` `freedrain=true` | `fortran/soilfluxes.f90` `freedrain=1` | `test/test_soil_fluxes.jl:L40` | `docs/土壤水运动_ASAP.typ` §边界 |
| 非自由排水底部边界 | `Q[1] = 0`（无深层出流） | `src/SoilFluxes.jl:180-L198` `freedrain=false` | `fortran/soilfluxes.f90` `freedrain=0` | `test/test_soil_fluxes.jl:L100 / L160` | 同上 |
| 入渗能力 | $I_{max} = K_{sat} \cdot \Delta t$；`runoff = max(0, |Q| − I_max)` | `src/SoilFluxes.jl:26-L42` 内联 | `fortran/soilfluxes.f90` 内联 | `test/test_soil_fluxes.jl` | — |
| Green-Ampt 改进 | $dV/dt = \frac{16 K(\theta_s)}{(\pi+4)(\theta_s-\theta_l)} \cdot \frac{H+0.5L-S_m}{L}$（王文焰 2003） | `src/SoilFluxes.jl` 内联（**未完整实现**，`docs/下渗_黄土高原.typ`） | **无对应** | — | `docs/下渗_黄土高原.typ` |

## 8. 浅层 / 深层地下水位更新

| 物理过程 | 数学描述 | Julia 实现 | Fortran 实现 | 测试位置 | 文档参考 |
|---|---|---|---|---|---|
| 浅层水位更新 | 离散跳跃（按饱和/过干判定升降） | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/updatewtd_shallow.jl:26-L102` `updatewtd_shallow` | `fortran/module_rootdepth.f90:1603 UPDATESHALLOWWTD` | `test/test_updatewtd_shallow.jl` | `docs/土壤水运动_ASAP.typ` §WTD |
| 深度因子 `cal_factor` | $\text{clamp}(\exp((z+1.5)/f_{depth}), 0.1, 1.0)$ | `src/updatewtd_shallow.jl:2-L25` `cal_factor` | `module_rootdepth.f90` 内联 | `test/test_updatewtd_shallow.jl:L52` | — |
| 深层水位更新 | `totwater = qlat + qspring + rech` 驱动 | `src/modules/Tracing/IsotopeTracing.jl:197-L279` `updatedeepwtable!` | `fortran/module_wtable.f90:185 UPDATEDEEPWTABLE` / `:468 UPDATEWTD` | `test/wtable/test_isotope_tracing.jl` | `wiki/julia/IsotopeTracing-同位素追踪.md` |
| 简化深层更新 | 累积 `bottomflux - rech` 差 | `src/modules/Tracing/IsotopeTracing.jl:280-L337` `updatewtd_simple` | **无对应** | — | 同上 |
| 地下水补给驱动 | `rech = Q[1]`（最底层通量） | `src/SoilFluxes.jl` 返回 `rech` | `fortran/soilfluxes.f90` 写 `rech` | `test/test_soil_fluxes.jl` | — |

## 9. 8 方向侧向地下水流（D8 达西）

| 物理过程 | 数学描述 | Julia 实现 | Fortran 实现 | 测试位置 | 文档参考 |
|---|---|---|---|---|---|
| D8 8 邻居差分 | 4 主向 + 4 对角（×1/√2 距离修正） | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/modules/lateral_flow.jl:8-L58` `lateral_flow!` | `fortran/module_wtable.f90:138 LATERAL` / `:269 LATERALFLOW` | `test/wtable/test_watertable.jl:L23` | `docs/侧向地下水流.typ` |
| 达西通量 | $q = -K \cdot dH/d\xi$ | `src/modules/lateral_flow.jl:21-L27` | `module_wtable.f90:269` | 同上 | 同上 |
| 有效导水率深度依赖 | $\kappa = f_{depth}\cdot\kappa_{lat}\cdot\exp((wtd+1.5)/f_{depth})$ | `src/modules/lateral_flow.jl:21-L27` 内联 | `module_wtable.f90:269 LATERALFLOW` | 同上 | — |
| 角度因子 | $\sqrt{\tan(4\pi/32)}/(2\sqrt{2})$ | `src/modules/lateral_flow.jl:7` `fangle` | `module_wtable.f90:138 LATERAL` | — | — |
| 4 主方向 + O18 | 北/南/西/东 + 纬度 cos 修正 | `src/modules/Tracing/IsotopeTracing.jl:20-L72` `lateral_isotope!` | `fortran/module_wtable.f90:344 LATERALFLOW4` | `test/wtable/test_isotope_tracing.jl:L46` | — |
| 流向编码（D8 → 下游像元） | `(fd[i,j] → (ii,jj))` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/helper.jl:32-L55` `flowdir` | `fortran/module_wtable.f90:1384 FLOWDIR` | `test/wtable/test_watertable.jl:L7` | — |

## 10. 地下水-河流交换（Miguez-Macho 2007）

| 物理过程 | 数学描述 | Julia 实现 | Fortran 实现 | 测试位置 | 文档参考 |
|---|---|---|---|---|---|
| 河床导水率衰减 | $K_{rb} = K_{sat} \cdot \exp(-d/f_{depth})$（沿深度指数） | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/modules/gw2river.jl:19-L72` `gw2river!` | `fortran/module_wtable.f90:744 GW2RIVER` | `test/wtable/test-river.jl` | `docs/汇流_扩散波.typ` |
| 每日排水上限 | `min(..., Δt·0.05/86400)` = 50 mm/day | `src/modules/gw2river.jl:42` | `module_wtable.f90:744` 同 | `test/wtable/test_river_routing.jl` | Miguez-Macho et al. 2007 |
| `qrf` 下游重分配 | 小河流 `qrf` 按 `fd` 推到下游 | `src/modules/gw2river.jl:73-L97` `moveqrf!` | `fortran/module_wtable.f90:1325 MOVEQRF` | `test/wtable/test_river_routing.jl` | — |
| 硬编码人工分流 | 按 `(i,j)` 索引加减下游流量 | `src/modules/gw2river.jl:98-L125` `apply_specific_diversions` | `module_wtable.f90` 内联（Taquari 河等） | `test/wtable/test_river_routing.jl` | — |
| 河床连通性判断 | `wtd > -maxdepth` → 河-含水层连通 | `src/modules/gw2river.jl:30-L40` 内联 | `module_wtable.f90:744` 内联 | `test/wtable/test-river.jl` | — |

## 11. 运动波河流路由（Manning）

| 物理过程 | 数学描述 | Julia 实现 | Fortran 实现 | 测试位置 | 文档参考 |
|---|---|---|---|---|---|
| Manning 流速 | $v = R^{2/3} S_0^{1/2} / n$ | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/modules/rivers_kw_flood.jl:34-L46` `rivers_kw_flood!` | `fortran/module_wtable.f90:800 RIVERS_KW_FLOOD` | `test/wtable/test_river_routing.jl` | `docs/汇流_扩散波.typ` §1 |
| 河道-洪泛区重分配 | 跨越 `riverchannel = maxdepth·riverarea` 阈值 | `src/modules/rivers_kw_flood.jl:30-L33` | `module_wtable.f90:800-L993` | 同上 | 同上 |
| 子步长 `dtlr` | 5 min 显式 KW 步长 | `src/modules/rivers_kw_flood.jl` 参数 `δt` | `module_wtable.f90:800` | 同上 | — |
| 流量累加 | `qmean += qnew * Δt` | `src/modules/rivers_kw_flood.jl:48` | `module_wtable.f90` | 同上 | — |

## 12. 扩散波河流路由

| 物理过程 | 数学描述 | Julia 实现 | Fortran 实现 | 测试位置 | 文档参考 |
|---|---|---|---|---|---|
| 扩散波方程 | 保留压力项、忽略惯性：$q^{n+1} = \gamma^n q^n + S_h$ | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/modules/rivers_dw_flood.jl:8-L100` `rivers_dw_flood!` | `fortran/module_wtable.f90:996 RIVERS_DW_FLOOD` | `test/wtable/test_rivers_dw_flood.jl` | `docs/汇流_扩散波.typ` §3 |
| 冻结系数线性化 | 摩阻项 $γ^n·q^{n+1}$（隐式） | `src/modules/rivers_dw_flood.jl:36-L46` | `module_wtable.f90:996-L1222` | 同上 | 同上 |
| 洪水态/非洪水态切换 | `floodheight > 0.05` 或下游超载 | `src/modules/rivers_dw_flood.jl:29-L32` | `module_wtable.f90:996` | 同上 | — |
| CFL 限制 | `clamp(speed, 0.01, length/δt)` | `src/modules/rivers_dw_flood.jl:46` | `module_wtable.f90` | — | — |

## 13. 洪泛区漫流

| 物理过程 | 数学描述 | Julia 实现 | Fortran 实现 | 测试位置 | 文档参考 |
|---|---|---|---|---|---|
| D8 漫流通量 | $dh = dh / \sqrt{2}$（对角线）+ 坡降驱动 | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/modules/flooding.jl:6-L70` `flooding!` | `fortran/module_wtable.f90:1226 FLOODING` | `test/wtable/test_river_routing.jl` | `docs/汇流_扩散波.typ` |
| 主河道扣除 | 最低邻即下游流向 `(i1,j1)` 时扣除河道容量 | `src/modules/flooding.jl:26-L32` | `module_wtable.f90:1226` | 同上 | — |
| 地表水变化约束 | `delsfcwat < 0` 时截断源端 | `src/modules/flooding.jl:34-L36` | `module_wtable.f90:1226` | 同上 | — |

## 14. ¹⁸O 同位素追踪

| 物理过程 | 数学描述 | Julia 实现 | Fortran 实现 | 测试位置 | 文档参考 |
|---|---|---|---|---|---|
| 侧向流同位素对流 | 4 方向通量乘以 `(o18wtd + o18wtd_neighbor)/2` | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/modules/Tracing/IsotopeTracing.jl:20-L72` `lateral_isotope!` | `fortran/module_wtable.f90:344 LATERALFLOW4`（参数 `o18wtd, qlato18`） | `test/wtable/test_isotope_tracing.jl:L86` | `wiki/julia/IsotopeTracing-同位素追踪.md` |
| 水位处液态浓度归一 | `o18wtd = o18wtd / smoiwtd` | `src/modules/Tracing/IsotopeTracing.jl:25-L30` | `module_wtable.f90:344` | 同上 | — |
| 累积 `*sum` | `qlatinsum += qlat*Δt * 1e3`（m→mm） | `src/modules/Tracing/IsotopeTracing.jl:60-L70` | `module_wtable.f90:344` | 同上 | — |
| 土壤层 ¹⁸O 对流 | `o18ratio = o18 ./ θ` | `src/SoilFluxes.jl` 末 30 行（**已注释**） | `fortran/module_rootdepth.f90:994 SOILFLUXES`（`transpo18` 同步累加） | — | `module_rootdepth.f90` |
| 平衡分馏 `α` | $\alpha = 1/\exp(1137/T^2 - 0.4156/T - 0.0020667)$（Majoube） | **无对应** | `fortran/module_rootdepth.f90` 内联 | — | Majoube 1971 |
| O18 降水气候态读取 | 月线性插值到日 | **无对应** | `fortran/module_forcings.f90:1174 READO18CLIM` | — | — |

## 15. 冰冻因子

| 物理过程 | 数学描述 | Julia 实现 | Fortran 实现 | 测试位置 | 文档参考 |
|---|---|---|---|---|---|
| 冰冻指示符 | `icefactor[k] ∈ {0,1}`，由 ERA5 4 层土壤温度判定 | 通过父模块作用域注入 | `fortran/module_forcings.f90:692 READFORCINGSSOILT` | — | — |
| 冰冻层水力抑制 | `K(θ) *= (液态水比)^{2b+3}` | `src/SoilFluxes.jl:52`（裸引用 `f_ice`） | `fortran/module_rootdepth.f90:994 SOILFLUXES` | — | — |
| 根系吸水跳过冰冻层 | `icefactor[k]==0` 时 `transp[k]=0` | `src/extraction.jl` 内联（`soilfactor_k`） | `module_rootdepth.f90:307 EXTRACTION` | — | — |

## 16. 主入口 / 时步循环

| 物理过程 | 数学描述 | Julia 实现 | Fortran 实现 | 测试位置 | 文档参考 |
|---|---|---|---|---|---|
| ASAP 主模块 | `module ASAP` 顶层 include | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/ASAP.jl:5-L27` | `fortran/module_driver.f90:1 program driver` | — | — |
| 时步主循环 | `while year != 2024`：LATERAL → GW2RIVER → ROOTDEPTH → FLOODING → RIVERS_KW_FLOOD | `src/RootDepth.jl:90-L278` `rootdepth_main` | `fortran/module_driver.f90` 时步段 | `test/test_rootdepth.jl` | `module_driver.f90` |
| 模块聚合 | 11 个 `!` 函数 export | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/modules/Modules.jl:1-L14` | 无对应（Fortran 子程序均 PUBLIC） | `test/wtable/test_watertable.jl` | — |
| `updatewtd_qlat` 历史实现 | 含侧向流的水位更新 | `/mnt/z/GitHub/jl-pkgs/ASAP-model/src/backup/updatewtd_qlat.jl`（已迁移） | `fortran/module_rootdepth.f90:1725 UPDATEWTDQLAT` | — | — |
| `wtable!` 主调度 | 侧向流 + 深层 + 河流全流程 | `src/modules/Modules.jl:11`（**悬空 export**） | `fortran/module_wtable.f90:13 WTABLE` | — | — |

## 17. NetCDF I/O / 强迫 / 并行（Fortran 独占）

| 物理过程 | 数学描述 | Julia 实现 | Fortran 实现 | 测试位置 | 文档参考 |
|---|---|---|---|---|---|
| 静态场读取 | `READINITIAL/READWTD/READVEG/READFLOWDIRECTION` 等 12 子程序 | **无对应** | `fortran/module_io.f90` | — | `wiki/fortran/module_io.md` |
| Restart 读取 | `READHISTORYNC` + `READHISTORYVAR3DNC` | **无对应** | `fortran/module_io.f90:823 / 1173` | — | 同上 |
| 日 / 月输出 | `WRITEOUTPUTNC_par` + `WRITEOUTPUTNC_INFIL_par` | **无对应** | `fortran/module_io.f90:2115 / 3311` | — | 同上 |
| 小时 ERA5 强迫 | 11 变量（风/温/压/露点/降水/辐射/4 层土温） | **无对应** | `fortran/module_forcings.f90:14 READFORCINGS` | — | `wiki/fortran/module_forcings.md` |
| 累积辐射 / 降水 | `READFORCINGSACC` 差分 | **无对应** | `fortran/module_forcings.f90:337` | — | 同上 |
| 融雪强迫 | `READFORCINGSSNOW`（ERA5-Land `smlt`） | **无对应** | `fortran/module_forcings.f90:557` | — | 同上 |
| MPI 域分解 | `DIVIDEDOMAIN` 按地形把陆地块分给 worker | **无对应** | `fortran/module_parallel.f90:252` | — | `wiki/fortran/module_parallel.md` |
| Halo 通信 | 上下后左右 `SENDBORDERS4` | **无对应** | `fortran/module_parallel.f90:558` | — | 同上 |
| 自定义 MPI 类型 | `domblock` / `arraysection` / `updownhalo` 等 11 个 | **无对应** | `fortran/module_parallel.f90`（CONTAINS 前） | — | 同上 |
| 网格插值 | `gdtost2` / `gdtost3` / `TRNCL1/2` / `HTINT2` 等 11 个 | **无对应** | `fortran/interp_lib.f90` | — | — |

## 18. 测试覆盖索引

| 测试文件 | 覆盖过程 |
|---|---|
| `/mnt/z/GitHub/jl-pkgs/ASAP-model/test/test_soil_parameters.jl` | §1 土壤参数化全部 |
| `/mnt/z/GitHub/jl-pkgs/ASAP-model/test/test_soil_initialization.jl` | §2 土壤层初始化 |
| `/mnt/z/GitHub/jl-pkgs/ASAP-model/test/test_evapotranspiration.jl` | §3 三种 PET |
| `/mnt/z/GitHub/jl-pkgs/ASAP-model/test/test_interception.jl` | §4 截留 |
| `/mnt/z/GitHub/jl-pkgs/ASAP-model/test/test_extraction.jl` | §5 根系吸水与水分胁迫 |
| `/mnt/z/GitHub/jl-pkgs/ASAP-model/test/test_soil_fluxes.jl` | §6 Richards + §7 边界 + §8 补给 |
| `/mnt/z/GitHub/jl-pkgs/ASAP-model/test/test_updatewtd_shallow.jl` | §8 浅层水位 |
| `/mnt/z/GitHub/jl-pkgs/ASAP-model/test/wtable/test_watertable.jl` | §9 8 方向侧向流 + D8 |
| `/mnt/z/GitHub/jl-pkgs/ASAP-model/test/wtable/test_isotope_tracing.jl` | §9 O18 侧向 + §8 深层 + §14 同位素 |
| `/mnt/z/GitHub/jl-pkgs/ASAP-model/test/wtable/test-river.jl` | §10 地下水-河流交换 |
| `/mnt/z/GitHub/jl-pkgs/ASAP-model/test/wtable/test_river_routing.jl` | §10 `moveqrf!` + §11 运动波 + §13 漫流 |
| `/mnt/z/GitHub/jl-pkgs/ASAP-model/test/wtable/test_rivers_dw_flood.jl` | §12 扩散波 |

## 19. 文档参考

- `docs/土壤水运动_ASAP.typ` — 主方程（Richards、PET、根系吸水、截留、WTD）
- `docs/下渗_黄土高原.typ` — Green-Ampt 椭圆面积修正（王文焰 2003）
- `docs/汇流_扩散波.typ` — 圣维南方程组、运动波与扩散波推导
- `docs/extraction.typ` — Feddes 根系吸水
- `docs/侧向地下水流.typ` — D8 差分细节
- `docs/cetz/Figure_Soil_ASAP.typ`、`Figure_Soil_Kong.typ` — 土壤参数图件
- `wiki/julia/` — 14 篇过程级 wiki 页（参见 `julia-fortran-对照.md` §0）
- `wiki/fortran/` — 7 篇模块级 wiki 页

## 20. 引用

- Julia 源：`/mnt/z/GitHub/jl-pkgs/ASAP-model/src/*.jl`、`src/modules/*.jl`、`src/modules/Tracing/*.jl`
- Fortran 源：`/mnt/z/GitHub/jl-pkgs/ASAP-model/fortran/*.f90`
- 测试根目录：`/mnt/z/GitHub/jl-pkgs/ASAP-model/test/`
- 既有 wiki：`/mnt/z/GitHub/jl-pkgs/ASAP-model/wiki/julia/`、`/mnt/z/GitHub/jl-pkgs/ASAP-model/wiki/fortran/`