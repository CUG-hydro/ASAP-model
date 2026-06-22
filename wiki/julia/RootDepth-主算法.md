# RootDepth — 主算法入口

> 源文件：`src/RootDepth.jl`（279 行）
> Fortran 来源：`fortran/module_rootdepth.f90`（≈1700 行，移植参考）
> 测试：`test/test_rootdepth.jl`、`test/test_soil_fluxes.jl`（间接）
> 状态：已摄取（含同位素追踪遗留问题）

## 1. 功能概述

`rootdepth_main` 是 ASAP-model 的**核心驱动函数**，按时间步长 `Δt` 对 (is:ie)×(js:je) 网格逐格点执行完整的水量平衡计算：Shuttleworth-Wallace PET、植被截留、根系吸水、一维 Richards 求解、浅层地下水位更新、洪水与侧向流分摊、入渗深度统计。函数**原地修改所有状态变量**（无显式返回值），通过类型参数化支持 `Float64`/`Float32` 双精度运行。

模块不显式以 `module RootDepth ... end` 包裹，是 `module ASAP` 内部第 18 行的 `include("RootDepth.jl")` 注入的顶层脚本。

## 2. 泛型函数签名

```julia
function rootdepth_main(
  freedrain::Int, is::Int, ie::Int, js::Int, je::Int, nzg::Int,
  z₋ₕ::V, dz::V, Δt::Float64,
  landmask::Matrix{Int}, veg::M, hveg::M,
  soiltxt::Array{Int,3}, wind::M, temp::M,
  qair::M, press::M, netrad::M,
  rshort::M, lai::M, precip::M,
  qsrun::M, θ::A3, θ_eq::A3,
  θ_wtd::M, wtd::M, waterdeficit::M,
  watext::A3,
  watextdeep::M, rech::M,
  deeprech::M,
  et_s::M, et_i::M, et_c::M,
  intercepstore::M, ppacum::M,
  pInfiltDepth::M, pInfiltDepthK_old::Matrix{Int8}, qlat::M,
  qlatsum::M, qsprings::M, inactivedays::Array{Int,3},
  maxinactivedays::Int, fdepth::M,
  steps::Float64, floodheight::M, qrf::M,
  delsfcwat::M, icefactor::Array{Int8,3}, wtdflux::M,
  et_s_daily::M, et_c_daily::M, transptop::M,
  infilk::Matrix{Int8}, infilflux::A3, infilfluxday::A3,
  infilcounter::Array{Int16,3}, hour::Int,
  o18::A3, o18ratiopp::M, tempsfc::M,
  qlato18::M, transpo18::M, upflux::A3) where {
  T<:Real,V<:Vector{T},M<:Matrix{T},A3<:Array{T,3}}
```

### 类型参数语义

| 参数             | 含义         | 约定                                            |
| ---------------- | ------------ | ----------------------------------------------- |
| `T<:Real`        | 数值类型     | `Float64` 或 `Float32`，决定整体精度            |
| `V<:Vector{T}`   | 一维数组     | `slz`、`dz`（层深度与厚度）                     |
| `M<:Matrix{T}`   | 二维网格数组 | 大部分气象强迫与状态变量，尺寸 `[ie, je]`       |
| `A3<:Array{T,3}` | 三维土壤柱   | `θ`、`θ_eq`、`watext`、`o18` 等 `[nzg, ie, je]` |

`Array{Int,3}`、`Matrix{Int8}`、`Array{Int16,3}` 等特殊整型数组保持**具体类型而非泛型**，因为它们用于 icefactor、infilk 等离散状态。

## 3. 主循环流程（6 步）

主循环由两层嵌套组成：外层 `for j in (js+1):(je-1); for i in (is+1):(ie-1)`，跳过 halo 边界；内层 `for itime in 1:round(Int, steps)` 进行 `steps` 个子步积分（默认 `steps=1.0`）。每格点按下列顺序执行 6 个步骤：

### 步骤 ① Shuttleworth-Wallace PET（L131-L142）

```julia
Δ, γ, λ, ra_a, ra_c, rs_c, R_a, R_s, petfactor_s, petfactor_c, petstep_w, petstep_i =
  potevap_shutteworth_wallace(Δt, temp[i,j], netrad[i,j], rshort[i,j],
    press[i,j], qair[i,j], wind[i,j], lai[i,j], veg[i,j], hveg[i,j], floodflag)
```

- 累计 `et_s[i,j] += petstep_w`（土壤/水面蒸发）；
- 若 `floodflag == 1 && veg ≤ 1`，扣减 `delsfcwat[i,j] -= petstep_w * 1e-3`（从洪水中蒸发）；
- 水体或裸地（`round(Int, veg) ≤ 1`）`continue`，跳过植被相关步骤。

### 步骤 ② 植被截留（L147-L152）

```julia
ppdrip, etstep_i, new_intercepstore = interception(
  precip[i,j], lai[i,j], intercepstore[i,j], petstep_i)
```

累计 `et_i[i,j]`、更新 `intercepstore`。详见 [`Interception-截留.md`](./Interception-截留.md)。

### 步骤 ③ 分步与通量预计算（L154-L163）

将每步外部强迫按 `steps` 切片：`ppdrip_step`、`floodstep`、`qlatstep`、`qrfstep`、`qlato18step`；初始化 `flux = zeros(Float64, nzg+1)` 与 `qlatflux = zeros(Float64, nzg+2)` 用于累计通量；记录 `wtd_old`、`dθ = zeros(Float64, nzg)` 预初始化（**修复**：避免 `UndefVarError`，见 commit `d8fdeb6`）。

### 步骤 ④ 子步循环（`for itime in 1:round(Int, steps)`，L169-L220）

依次执行：

1. **根系吸水** `extraction(...)`：返回 `(pet_s, pet_c, watdef, dθ, dsmoideep)`，更新 `et_c`、`waterdeficit`、`watext`、`transptop`、`et_c_daily`。详见 [`extraction-根系吸水.md`](./extraction-根系吸水.md)。
2. **土壤水流** `soilfluxes(...)`：返回 `(et_s_step, runoff, rechstep, flux_step, qrfcorrect, updated_smoi)`，内部调用 Thomas 三对角求解。详见 [`SoilFluxes-土壤水运动.md`](./SoilFluxes-土壤水运动.md)。**同位素追踪代码已被注释**（详见下文 §5）。
3. **状态累加**：把 `petstep_w`、`runoff`、`rechstep`、`et_s_step`、`qrfcorrect` 累加到 `delsfcwat`、`qsrun`、`rech`、`et_s`、`qrf`；`θ[:,i,j] .= updated_smoi` 覆盖含水量。
4. **浅层 WTD 更新** `updatewtd_shallow(...)`：返回 `(wtd, rech_additional)`，累加补给；详见 [`updatewtd_shallow-浅层水位.md`](./updatewtd_shallow-浅层水位.md)。
5. **通量累加**：`flux .+= flux_step`，供步骤 ⑤ 入渗深度判定使用。

### 步骤 ⑤ 入渗深度统计（L222-L274）

遍历 `k = nzg:-1:1`：

- `flux[k+1] < 0`（向下入渗）累加到 `infilflux[k,i,j]`；
- 否则累加到 `upflux[k,i,j]`（上升通量）；
- `infilfluxday[k,i,j]` 无条件累加。

`hour == 0`（每日 0 点）触发 `infilcounter` 自增与 `infilfluxday` 清零。

随后判定入渗深度阈值 `-0.333e-5 m/s`（≈ 0.012 cm/h × 步长调整），从 `k = nzg` 向下扫描确定 `_kInfilt` 与 `_pInfiltDepth`，写回 `pInfiltDepthK_old`、`pInfiltDepth`、`infilk`，并按条件累加 `wtdflux`。

### 步骤 ⑥ 全局守恒收尾（L271-L274）

```julia
if z₋ₕ[max(_kInfilt - 1, 1)] <= wtd_old
  wtdflux[i,j] -= flux[_kInfilt] * 1.0e3
end
infilk[i,j] > _kInfilt && (infilk[i,j] = _kInfilt)
```

`rootdepth_main` 返回 `nothing`，所有变量都是**原地更新**。

## 4. 关键变量与单位

| 符号                   | 含义                               | 单位       | 数组形状                     |
| ---------------------- | ---------------------------------- | ---------- | ---------------------------- |
| `freedrain`            | 自由排水标志                       | 0/1        | 标量                         |
| `is/ie/js/je`          | 网格范围（halo）                   | —          | 4×Int                        |
| `nzg`                  | 土壤层数                           | —          | Int（典型 15）               |
| `z₋ₕ`                  | 层下边界深度（地表为 0，地下为负） | m          | `[nzg+1]`                    |
| `dz`                   | 层厚度                             | m          | `[nzg]`                      |
| `Δt`                   | 时间步长（典型 3600 s）            | s          | Float64                      |
| `landmask`             | 陆地掩码                           | 0/1        | `[ie, je]`                   |
| `veg` / `hveg`         | 植被类型 / 高度                    | — / m      | `[ie, je]`                   |
| `soiltxt[1,i,j]`       | 当前格点土壤类型 (1..13)           | —          | `Array{Int,3}` 第 1 维 1/2   |
| `θ` / `θ_eq` / `θ_wtd` | 体积含水量 / 平衡 / 水位处         | m³/m³      | `[nzg, ie, je]` / `[ie, je]` |
| `wtd`                  | 地下水位深度（地表为 0）           | m          | `[ie, je]`                   |
| `petstep_w/s/c/i`      | SW 四分量（水面/土壤/冠层/截留）   | mm         | 局部变量                     |
| `flux`                 | 层间通量（向下为负）               | m          | `[nzg+1]`                    |
| `infilflux` / `upflux` | 入渗 / 上升通量                    | mm         | `[nzg, ie, je]`              |
| `o18`                  | ¹⁸O 含量（三维）                   | ‰ 或 m³/m³ | `[nzg, ie, je]`              |

## 5. 同位素追踪代码被注释的问题

`rootdepth_main` 内部调用 `soilfluxes(...)` 的位置（L186-L195）原本应接收 4 个返回值：

```julia
# 当前（实际）签名（L187-L193）：
et_s_step, runoff, rechstep, flux_step, qrfcorrect, updated_smoi =
  soilfluxes(nzg, Δt/steps, z₋ₕ, dz, soiltxt[1,i,j], θ_wtd[i,j],
    dθ, dsmoideep, θ[:,i,j], wtd[i,j], ppdrip_step, pet_s,
    fdepth[i,j], qlatstep, qrfstep, floodstep, icefac
    ; freedrain=Bool(freedrain))

# 被注释掉的返回值（L194-L195）：
# θ_eq[:, i, j], o18[:, i, j], o18ratiopp[i, j], tempsfc[i, j],
# qlato18step, transpo18step

# 被注释掉的累加（L202、L209）：
# transpo18[i, j] += transpo18step * 1.0e3
# o18[:, i, j] .= updated_o18
```

**问题清单**：

1. **`soilfluxes` 当前实现不再返回 `updated_o18` 与 `transpo18step`**：`src/SoilFluxes.jl` 末 30 行的氧 18 同位素段（`o18ratio = o18 ./ θ` 等）全部被注释，导致 `soilfluxes` 返回值仅 6 个而非 10 个。
2. **`θ_eq` 与 `tempsfc` 也未回写**：Fortran 原型中这两个变量在每次 `soilfluxes` 调用后被原地更新；Julia 版去掉此行为后，调用方需另行初始化。
3. **`transpo18` 与 `o18` 在主循环中未演化**：由于 (1)，输入的 `o18`/`o18ratiopp`/`tempsfc`/`qlato18`/`transpo18` 始终保持初始值，**同位素追踪在 `rootdepth_main` 主算法层面失效**。
4. **`IsotopeTracing.jl` 已建立但未被 `rootdepth_main` 调用**：`src/modules/Tracing/IsotopeTracing.jl` 提供 `lateral_isotope!`、`updatedeepwtable!`，是 `RootDepth.jl` 同位素路径的替代实现，但 `rootdepth_main` 仅在 `modules/Modules.jl` 聚合后供下游耦合调用，主循环未引用。

**恢复路径建议**：

- 短期：在 `RootDepth.jl` 中显式调用 `lateral_isotope!` 处理 `o18` 侧向传输、`updatedeepwtable!` 处理深层 `o18` 更新；
- 中期：恢复 `SoilFluxes.jl` 注释段，把 `o18`/`θ_eq`/`tempsfc` 回写重新纳入 `soilfluxes` 返回元组；
- 长期：在 `conventions.md` §2 同位素单位统一后，建立 `o18` 守恒测试（参考 `test/wtable/test_isotope_tracing.jl:L86-L116` 的 Isotope Conservation 断言）。

## 6. 与 Fortran 对照

| Fortran（`module_rootdepth.f90`）         | Julia（`rootdepth_main`）                                 | 差异                                     |
| ----------------------------------------- | --------------------------------------------------------- | ---------------------------------------- |
| `subroutine rootdepth(...)` 入口          | 同名泛型函数                                              | Julia 拆为 5 个 `include` 文件再 include |
| `COMMON /VEGTYPE/`、`/SOILTYPE/` 等全局块 | 形参列表（40+ 变量）                                      | Julia 显式传参，避免全局副作用           |
| `REAL*8` 固定双精度                       | `where T<:Real`                                           | Float32 路径需自行验证                   |
| `INTEGER jwt` 作为入参                    | 通过 `find_jwt(wtd, z₋ₕ)` 内部反算                        | Julia 调用前需确保 `z₋ₕ` 已初始化        |
| 内联调用 5+ 个子程序                      | 内联调用 PET/截留/extraction/soilfluxes/updatewtd_shallow | 同结构，函数命名 snake_case              |
| 氧 18 计算主循环内联                      | 同位素段注释掉                                            | 见 §5                                    |
| `DO j = js+1, je-1`                       | `for j in (js+1):(je-1)`                                  | 1-based，halo 边界同                     |

## 7. 引用

- 函数定义：`src/RootDepth.jl` L90-L278
- 泛型 where 子句：L113-L114
- 主循环入口：L121-L122
- Shuttleworth-Wallace 调用：L132-L135
- 截留调用：L147-L148
- 子步循环：L169-L220
- 浅层 WTD 更新调用：L212-L213
- 入渗深度判定：L222-L274
- 同位素注释点：L186、L194-L195、L202、L209
- 示例调用：`example/complete_example.jl` L138-L149
- 文档参考：`docs/土壤水运动_ASAP.typ`、`docs/汇流_扩散波.typ`

## 8. 已知问题与备注

1. **`dθ` 预初始化**（L166）：commit `d8fdeb6` 已修复 `for k in nzg:-1:0` 块中 `dθ[k]` 在 `itime == 0` 时访问的 `UndefVarError`，保留作为安全网。
2. **`rech` 双重累加**：浅层 WTD 返回的 `rech_additional * 1.0e3`（L215）与 `soilfluxes` 返回的 `rechstep * 1.0e3`（L200）相加；Fortran 中该累加来自不同子程序，调用方需核对合理性。
3. ✅ **`length` 命名冲突已修复**（2026-06-22 修复）：`src/modules/rivers_kw_flood.jl:18` 与 `src/modules/rivers_dw_flood.jl:13` 的形参均已重命名为 `river_length`；`helper.jl`/`RootDepth.jl` 中使用 `length(z₋ₕ)` 的局部绑定与 `Base.length` 同名但属函数体内局部作用域，未实际触发冲突。
4. ⏸️ **`o18`/`θ_eq`/`tempsfc` 形参未消费**（按 §6.1 暂缓）：`src/SoilFluxes.jl:347-L351` ¹⁸O 三对角组装段被注释；恢复需用户明确重新授权并追加 `wiki/log.md` 的 `[YYYY-MM-DD] enable | 范围` 条目。
5. **`steps = 1.0` 默认**：典型调用为 `steps = 1.0`，若调大需保证 `Δt/steps` 仍满足 CFL 与 `icefactor` 时间分辨率。
6. **`icefactor` 维度**：`Array{Int8,3}` 的最后一维长度为 26-40（注释 L57 提示实际为 26:40 切片 15 个元素），主循环 L126 复制 `icefactor[i,j,26:40]` 到 `icefac`，调用方需保证该切片存在。
7. **测试入口**：`test/test_rootdepth.jl` 应至少包含：① 单格点水分守恒（`sum(θ_new) ≈ sum(θ_old) + precip - et - runoff`）、② 数值精度切换（`Float32` 与 `Float64` 结果一致至 1e-4）、③ 同位素守恒断言（待 §5 修复后启用）。
