# 跨页面引用登记表

> 记录至少 10 条符号级别的跨页面引用。每条说明**源页面 → 目标符号 → 目标页面**，便于读者从一个模块的公式或函数名直接跳转。

## 1. 数据流引用（主循环依赖）

| # | 源页面 | 源符号 / 调用点 | 目标符号 | 目标页面 |
|---|---|---|---|---|
| 1 | `RootDepth-主算法.md` 步骤 ① | `potevap_shutteworth_wallace(...)` | Shuttleworth-Wallace 12 元组（Δ/γ/λ/ra_a/ra_c/rs_c/R_a/R_s/petfactor_s/petfactor_c/petstep_w/petstep_i） | [`Evapotranspiration-蒸散发.md`](../julia/Evapotranspiration-蒸散发.md) §2、§3 |
| 2 | `RootDepth-主算法.md` 步骤 ② | `interception(...)` | `ppdrip, etstep_i, new_intercepstore` | [`Interception-截留.md`](../julia/Interception-截留.md) §2、§3 |
| 3 | `RootDepth-主算法.md` 步骤 ④.1 | `extraction(...)` | `pet_s, pet_c, watdef, dθ, dsmoideep` | [`extraction-根系吸水.md`](../julia/extraction-根系吸水.md) §2、§3 |
| 4 | `RootDepth-主算法.md` 步骤 ④.2 | `soilfluxes(...)` | `et_s_step, runoff, rechstep, flux_step, qrfcorrect, updated_smoi` | [`SoilFluxes-土壤水运动.md`](../julia/SoilFluxes-土壤水运动.md) §2、§3 |
| 5 | `RootDepth-主算法.md` 步骤 ④.4 | `updatewtd_shallow(...)` | `(wtd, rech_additional)` | [`updatewtd_shallow-浅层水位.md`](../julia/updatewtd_shallow-浅层水位.md) §2、§3 |
| 6 | `RootDepth-主算法.md` 步骤 ④.1（内部 `find_jwt`） | `find_jwt(wtd, z₋ₕ)` | 水位所在层索引 | [`ASAP-主入口.md`](../julia/ASAP-主入口.md) §5 + `src/helper.jl:L4-L14` |

## 2. 公式 / 参数级引用（跨页面共享）

| # | 源页面 | 源公式 / 参数 | 目标符号 | 目标页面 |
|---|---|---|---|---|
| 7 | `SoilFluxes-土壤水运动.md` §3.3 | Campbell 导水率 `K(θ) = Ksat (θ/θsat)^(2b+3)` | `cal_K(θ, θ_sat, Ksat, b)` | [`SoilParameters-土壤参数.md`](../julia/SoilParameters-土壤参数.md) §2、§3 |
| 8 | `SoilFluxes-土壤水运动.md` §3.5 | Green-Ampt 椭圆面积修正 `16/(π+4)` 系数 | 引用 `docs/下渗_黄土高原.typ`（王文焰 2003） | 同上 §3.5 注释 |
| 9 | `extraction-根系吸水.md` §3 | `POTLEAF = POTWILT = -153.0 m`、`POTFC = -3.366 m` | 凋萎/田间持水量水势阈值 | [`SoilParameters-土壤参数.md`](../julia/SoilParameters-土壤参数.md) §4（POTWILT_LOCAL） |
| 10 | `SoilFluxes-土壤水运动.md` §3.3 | Ksat/θsat/ψsat 深度 `exp(z/fdepth)` 缩放 | `fdepth` 双重语义 | [`extraction-根系吸水.md`](../julia/extraction-根系吸水.md) §4 + [`updatewtd_shallow-浅层水位.md`](../julia/updatewtd_shallow-浅层水位.md) §3 |

## 3. 模块聚合引用（`Modules.jl` 出口）

| # | 源页面 | 源符号 | 目标符号 | 目标页面 |
|---|---|---|---|---|
| 11 | `modules-水位模块聚合.md` §2 | `lateral_flow!` | D8 八邻居达西通量 | [`lateral_flow-侧向地下水流.md`](../julia/lateral_flow-侧向地下水流.md) §2、§3 |
| 12 | `modules-水位模块聚合.md` §2 | `gw2river!` | 地下水-河流双向通量 | [`gw2river-地下水河流交换.md`](../julia/gw2river-地下水河流交换.md) §2、§3 |
| 13 | `modules-水位模块聚合.md` §2 | `rivers_kw_flood!` / `rivers_dw_flood!` | 运动波 / 扩散波路由 | [`rivers_kw_flood-运动波路由.md`](../julia/rivers_kw_flood-运动波路由.md) + [`rivers_dw_flood-扩散波路由.md`](../julia/rivers_dw_flood-扩散波路由.md) |
| 14 | `modules-水位模块聚合.md` §2 | `flooding!` | 洪泛区漫流 | [`flooding-洪泛漫流.md`](../julia/flooding-洪泛漫流.md) §2、§3 |
| 15 | `modules-水位模块聚合.md` §2 | `lateral_isotope!` / `updatedeepwtable!` | ¹⁸O 侧向流 + 深层水位 | [`IsotopeTracing-同位素追踪.md`](../julia/IsotopeTracing-同位素追踪.md) §2、§3 |
| 16 | `modules-水位模块聚合.md` §7 | `wtable!` / `updatewtd!` | 悬空 export | [`_meta/status.md`](./status.md) §5（全局问题汇总） |

## 4. 河流链引用

| # | 源页面 | 源符号 | 目标符号 | 目标页面 |
|---|---|---|---|---|
| 17 | `lateral_flow-侧向地下水流.md` §1 | D8 侧向流 → `qlat` | `gw2river!` 中读取 `qlat`/传导度 | [`gw2river-地下水河流交换.md`](../julia/gw2river-地下水河流交换.md) §3 |
| 18 | `gw2river-地下水河流交换.md` §3 | `qrf` 双向通量 | `moveqrf!` 下游重分配 | 同 §2 `moveqrf!` 签名 |
| 19 | `gw2river-地下水河流交换.md` §3 | `qrf` 输入到 `rivers_kw_flood!` | KW 河道汇入 | [`rivers_kw_flood-运动波路由.md`](../julia/rivers_kw_flood-运动波路由.md) §2 |
| 20 | `rivers_kw_flood-运动波路由.md` §3 | Manning 公式 `V = 1/n R^(2/3) √S0` | 同样公式在 DW 中回退使用 | [`rivers_dw_flood-扩散波路由.md`](../julia/rivers_dw_flood-扩散波路由.md) §3 |
| 21 | `rivers_kw_flood-运动波路由.md` §3 + `rivers_dw_flood-扩散波路由.md` §3 | `floodheight` 跨函数传递 | 漫流通量约束 | [`flooding-洪泛漫流.md`](../julia/flooding-洪泛漫流.md) §3 |

## 5. 同位素链引用

| # | 源页面 | 源符号 | 目标符号 | 目标页面 |
|---|---|---|---|---|
| 22 | `IsotopeTracing-同位素追踪.md` §3 | `lateral_isotope!` 4 方向通量 | 基础 D8 格式 | [`lateral_flow-侧向地下水流.md`](../julia/lateral_flow-侧向地下水流.md) §3（8 邻居对比） |
| 23 | `IsotopeTracing-同位素追踪.md` §3 | `o18` / `o18ratiopp` 形参 | 主循环中未消费 | [`RootDepth-主算法.md`](../julia/RootDepth-主算法.md) §5（同位素遗留问题） |
| 24 | `SoilFluxes-土壤水运动.md` §7（已知问题） | `o18ratio = o18 ./ θ` 注释段 | 同 §5 待恢复 | [`RootDepth-主算法.md`](../julia/RootDepth-主算法.md) §5 |

## 6. 命名 / 单位引用（指向 `conventions.md`）

| # | 源页面 | 符号 | 引用位置 |
|---|---|---|---|
| 25 | `SoilFluxes-土壤水运动.md` §4 | `K_mid, D_mid, Q, flux, Imax, rech` 单位 m/s, m²/s, m | [`conventions.md` §2](../conventions.md) |
| 26 | `gw2river-地下水河流交换.md` §4 | `qrf` 单位 m | [`conventions.md` §2](../conventions.md) |
| 27 | `RootDepth-主算法.md` §4 | `θ/θ_eq/θ_wtd/wtd/flux/infilflux` 单位 m³/m³, m, mm | [`conventions.md` §2](../conventions.md) |
| 28 | `IsotopeTracing-同位素追踪.md` §4 | `o18/o18wtd` 单位 ‰ | [`conventions.md` §2](../conventions.md) |
| 29 | `Forcings-ERA5.md` §3 | `read_initial` / `read_wtdnc` 静态场 + 初始水位读取 | [`io-NetCDF.md`](../julia/io-NetCDF.md) §2、§3 |
| 30 | `Forcings-ERA5.md` §6 | ERA5 路径约定（变量分目录、每日一文件） | [`README.md` §区域应用：数据准备清单](../../README.md#区域应用数据准备清单) |

| # | 源页面 | 源符号 / 调用点 | 目标符号 | 目标页面 |
|---|---|---|---|---|
| 31 | `example-regional.md` §3.5 | `eqsoilmoisturetheor(nzg, nsoil, z₋ₕ, dz, fdepth, wtd)`（2-D 循环） | 平衡含水量逐层 `zbrent` 求根 | [`SoilInitialization-土壤分层.md`](../julia/SoilInitialization-土壤分层.md) §3 |
| 32 | `example-regional.md` §3.6 | `rootdepth_main(...)` 47 形参对齐 | 主算法泛型签名 + 6 步主循环 | [`RootDepth-主算法.md`](../julia/RootDepth-主算法.md) §1、§2 |
| 33 | `example-regional.md` §5 | `read_initial(paths.static)` / `read_wtdnc(paths.wtd)` | NetCDF 读取 + 符号约定 | [`io-NetCDF.md`](../julia/io-NetCDF.md) §2、§3 |
| 34 | `example-regional.md` §3.4 | `era5_paths_for(era5_root, date)` 多日滚动 + `read_mock_hourly_forcings` | 与 `read_hourly_forcings` 的 4 点差异 | [`Forcings-ERA5.md`](../julia/Forcings-ERA5.md) §2、§7 |
| 35 | `example-regional.md` §7 #2 | `icefactor .= 0` 简化路径 | 真实冻结因子 `read_soil_temps` 4 层 STL 映射 | [`Forcings-ERA5.md`](../julia/Forcings-ERA5.md) §3.7 |
| 36 | `example-regional.md` §6（测试引用） | `test/test_regional_example.jl` 4 个 testset（端到端 / 静态场往返 / wtd 符号 / 类型断言） | `Pkg.test()` 入口与断言粒度 | [`CLAUDE.md` §8](../../CLAUDE.md) |
| 37 | `example-regional.md` §1（数据准备） | `static.nc` / `wtd.nc` / ERA5 / LAI 文件约定 | README 数据准备清单 | [`README.md` §区域应用：数据准备清单](../../README.md#区域应用数据准备清单) |

## 引用统计

- **总计登记**：37 条
- **覆盖页面**：16 个 Julia 页面（新增 `example-regional.md`）+ `_meta/status.md` + `_meta/cross-refs.md` + `conventions.md` + `index.md` + `log.md` + `README.md` + `CLAUDE.md`
- **主要链路**：① RootDepth 主循环 → 5 个子模块（条目 1-6）；② 公式参数共享（条目 7-10）；③ Modules.jl 聚合（条目 11-16）；④ 河流链（条目 17-21）；⑤ 同位素链（条目 22-24）；⑥ 单位命名引用（条目 25-28）；⑦ ERA5 强迫 + 数据准备（条目 29-30）；⑧ Regional 区域应用串联 7 条新引用（条目 31-37），覆盖 `eqsoilmoisturetheor` / `rootdepth_main` / `read_initial`+`read_wtdnc` / `read_mock_hourly_forcings` 差异 / `icefactor` 简化 / 测试入口 / 数据准备清单。

如发现新增跨页面引用，请按上述格式追加本表，并更新最后一节统计。