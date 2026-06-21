# module_rootdepth.f90

**路径**：`/mnt/z/GitHub/jl-pkgs/ASAP-model/fortran/module_rootdepth.f90`
**行数**：1914
**Module 声明**：`MODULE module_rootdepth`
**依赖**：`use module_parallel`（MPI 通信、全局数组）、`use module_nrtype`（`SP`/`DP`/`I4B`）

本模块是 ASAP 模型的核心，承担整个土壤-植被-大气界面的水热计算：冠层截留、潜在蒸散、1D Richards 求解、根系吸水、地下水埋深更新、O18 同位素示踪。

## 13 种 USDA 土壤参数（参数数组）

所有数组按 `nstyp=13`（USDA 土壤分类）下标访问。深度依赖衰减通过 `fdepth(i,j)` 实现：`slmsts(k)=slmsts(nsoil)*max(min(exp((slz(k)+1.5)/fdepth),1.),0.1)`。

| 数组 | 含义 | 单位 |
|------|------|------|
| `slmsts(nstyp)` | 饱和含水量 θ_sat | m³/m³ |
| `soilcp(nstyp)` | 凋萎/残余含水量 θ_wp | m³/m³ |
| `slbs(nstyp)` | Clapp & Hornberger b 指数 | – |
| `slcons(nstyp)` | 饱和水力传导度 K_sat | m/s |
| `slpots(nstyp)` | 冒泡压力 ψ_sat | m |
| `slwilt(nstyp)` | 凋萎点含水量 | m³/m³ |
| `klatfactor(nstyp)` | 侧向流比例因子 | – |
| `fieldcp(nzg,nstyp)` | 现场实测田持含水量 | m³/m³ |
| `slporos`, `slfc`, `slwp` | 孔隙度、田持、凋萎 | – |

## 关键 Subroutine / Function 签名

### ROOTDEPTH（主过程，第 32 行）

```fortran
SUBROUTINE ROOTDEPTH(freedrain, is, ie, js, je, nzg, slz, dz, deltat,
   landmask, veg, hveg, soiltxt, wind, temp, qair, press, netrad, rshort,
   lai, precip, qsrun, smoi, smoieq, smoiwtd, wtd, waterdeficit,
   watext, watextdeep, rech, deeprech, et_s, et_i, et_c, intercepstore,
   ppacum, pppendepth, pppendepthold, qlat, qslat, qsprings,
   inactivedays, maxinactivedays, fieldcp, fdepth, steps, floodheight,
   qrf, delsfcwat, icefactor, wtdflux, et_s_daily, et_c_daily,
   transptop, infilk, infilflux, infilfluxday, infilcounter, hour,
   o18, o18ratiopp, tempsfc, qlato18, transpo18, upflux)
```

每时步顶层入口，按 `freedrain` 标志切换自由排水/地下水耦合模式。内部顺序调用 `INTERCEPTION` → `EXTRACTION` → `POTEVAP_*` → `SOILFLUXES` → `UPDATESHALLOWWTD` / `UPDATEWTDQLAT`。

### SOILFLUXES（第 994 行；定义于 `soilfluxes.f90`）

```fortran
subroutine SOILFLUXES(i, j, nzg, freedrain, dtll, slz, dz, soiltxt,
   smoiwtd, transp, transpdeep, smoi, wtd, rech, deeprech, precip,
   pet_s, et_s, runoff, flux, fdepth, qlat, qlatflux, qrf, qrfcorrect,
   flood, icefactor, smoieq, o18, precipo18, tempsfc, qlato18, transpo18)
```

1D Richards 方程隐式求解：组装三对角矩阵 `aa/bb/cc/rr`，调用 `tridag`，再校正含水量不超过 `slmsts` 或低于 `soilcp`。同步更新 O18 浓度（`o18` 数组），使用平衡分馏系数 `alpha = 1/exp(1137/T²-0.4156/T-0.0020667)`。

### EXTRACTION（第 307 行）

```fortran
subroutine EXTRACTION(i, j, nzg, slz, dz, deltat, soiltxt, wtd, smoi,
   smoiwtd, veg, hveg, lai, et_c, pet_c, transp, transpdeep, fdepth,
   waterdeficit, watext, watextdeep, inactivedays, maxinactivedays,
   o18, transpo18, o18ratiopp)
```

根系吸水：使用 Feddes 减少函数（`fswp`）+ 易度函数（`rootfc`），对每层土壤求根系吸力 `transp(k)`，累加 `et_c`。遇冰冻层（`icefactor`）跳过。同位素通量 `transpo18` 同步累加。

### INTERCEPTION（第 277 行）

```fortran
subroutine INTERCEPTION(minpprate, precip, lai, intercepstore,
   ppdrip, pet_i, et_i)
```

冠层截留：LAI 决定最大截留量，超过则产生穿透降水 `ppdrip`，冠层蒸发 `et_i`。

### POTEVAP_Priestly_Taylor（第 546 行）

```fortran
subroutine POTEVAP_Priestly_Taylor(i, j, tempk, rad, presshp, pet)
```

基于净辐射 `rad` 与温度的经验 PET。

### POTEVAP_Penman_Monteith（第 569 行）

```fortran
subroutine POTEVAP_Penman_Monteith(i, j, tempk, rad, rshort, press,
   qair, wind, lai, veg, hveg, pet)
```

完整 P-M 公式：气孔阻抗 `rs_c_factor` 与 LAI、植被类型相关。

### POTEVAP_Shutteworth_Wallace（第 653 行）

```fortran
subroutine POTEVAP_Shutteworth_Wallace(i, j, deltat, tempk, rad,
   rshort, press, qair, wind, lai, veg, hhveg, ...)
```

S-W 双源（冠层+土壤）蒸散模型。

### TRIDAG（第 1868 行）

```fortran
SUBROUTINE tridag(a, b, c, r, u, n)
```

Thomas 算法求三对角线性方程组 `A·u=r`，隐式求解 1D Richards 的核心工具。

### UPDATESHALLOWWTD / UPDATEWTDQLAT（第 1603 / 1725 行）

地下水埋深 `wtd` 更新：根据土壤含水量分布、侧向流 `qlat`、重力排水 `rech` 重算水位。

### INITIALIZESOILDEPTH / INITIALIZESOILDEPTHCLM（第 973 / 937 行）

按 CLM 方案生成 40 层 `slz`/`dz` 节点深度数组。

## 全局状态

- `nstyp=13`（参数化宏）
- `nzg=40`（土壤层数，由 `module_driver` 传入）
- 通过 `module_parallel` 共享 MPI 进程信息
- `o18(nzg,is:ie,js:je)`：O18 浓度场
- `intercepstore(is:ie,js:je)`：冠层截留水库
- `inactivedays(0:nzg+1,...)`：每层累积失活天数

## 关键算法细节

1. **隐式时间积分**：土壤水势与含水量用 `tridag` 同时求，1h 步长稳定。
2. **冰冻层处理**：`icefactor` 控制 `hydcon` 与 `diffmid` 的乘子（液态水分数的 `2b+3` 次方）。
3. **O18 平衡分馏**：表层蒸发使用 Majoube 经验公式（`alpha`），与 `transpo18`、`qlato18` 双重守恒。
4. **垂向积分保护**：超饱和的水量向上/下层强制重分配，超出顶层则记入 `runoff`。
