# 蒸散发 (Evapotranspiration)

> 源文件：src/Evapotranspiration.jl:L1-L200+
> Fortran 来源：fortran/soilfluxes.f90（PET 子程序，对应 ASAP7 主程序）
> 测试：test/test_evapotranspiration.jl
> 状态：已摄取

## 1. 功能概述

本模块是 ASAP-model 中**唯一显式使用 `module ... end` 包裹**的源文件，实现了三种潜在蒸散发（PET）的计算方法：最简化的 Priestley-Taylor、单源的 Penman-Monteith（FAO 风格，含气孔/空气动力学阻力）、以及完整的 Shuttleworth-Wallace 双源模型。模块通过常量数组 `BIOPARMS`（31 类 USGS 土地覆盖的 [z0m, 叶宽]）耦合地表粗糙度参数，并返回土壤/冠层蒸发与截留蒸发分量，用于驱动水量平衡。

## 2. 函数签名

```julia
function potevap_priestly_taylor(tempk::Float64, rad::Float64, presshp::Float64) :: Float64

function potevap_penman_monteith(tempk::Float64, rad::Float64, rshort::Float64,
    press::Float64, qair::Float64, wind::Float64, lai::Float64,
    veg::Float64, hveg::Float64) :: Float64

function potevap_shutteworth_wallace(Δt::Float64, tempk::Float64, rad::Float64,
    rshort::Float64, press::Float64, qair::Float64, wind::Float64,
    lai::Float64, veg::Float64, hhveg::Float64, floodflag::Int) ::
    Tuple{Float64, Float64, Float64, Float64, Float64, Float64, Float64,
          Float64, Float64, Float64, Float64, Float64}
```

Shuttleworth-Wallace 返回 12 元组：`(Δ, γ, λ, ra_a, ra_c, rs_c, R_a, R_s, pet_s, pet_c, pet_w, pet_i)`，分别对应斜率/湿度计常数/潜热、空气/冠层边界层与冠层气孔阻力、对应 Penman-Monteith 阻力项，以及**土壤蒸发、冠层蒸腾、水面蒸发、截留蒸发**四个分量。

## 3. 算法 / 公式

**(a) Priestley-Taylor（最简化，仅需净辐射与温度）**

$$PET = \alpha \frac{\Delta}{\Delta + \gamma} \frac{R_n}{\lambda} \cdot \Delta t$$

其中 $\alpha \approx 1.26$，$\Delta$ 为饱和水汽压-温度曲线斜率，$\gamma$ 为湿度计常数，$R_n$ 由 W/m² 换算为 MJ/day/m²。

**(b) Penman-Monteith（FAO-56，含气孔阻力 $r_c$）**

$$\lambda ET = \frac{\Delta (R_n - G) + \rho_a c_p \frac{VPD}{r_a}}{\Delta + \gamma \left(1 + \frac{r_c}{r_a}\right)}$$

冠层阻力 $r_c$ 通过 Beer-Lambert 辐射衰减求和得到：$r_c = r_{l,min}/(0.5 LAI)$，并由短波辐射调节有效叶面比例 $f_{rad}=\min(1,\,(0.004 R_s+0.05)/(0.81(1+0.004 R_s)))$。空气动力学阻力 $r_a$ 用对数风廓线（高度 $z_m=10\,\mathrm{m}$，$z_0=0.1 h_{veg}$）。

**(c) Shuttleworth-Wallace 双源（耦合冠层 + 地表）**

$$\lambda E = C_c\,PM_c + C_s\,PM_s$$

其中

$$C_c = \left[1 + \frac{R_a R_c}{R_s(R_c+R_a)}\right]^{-1}, \quad
C_s = \left[1 + \frac{R_a R_s}{R_c(R_s+R_a)}\right]^{-1}$$

阻力网络 $R_a, R_c, R_s$ 由冠层边界层阻力 $r_{b,c}$、气孔阻力 $r_{s,c}$（叶宽参数化）、地表阻力 $r_{s,s}$（依赖 $LAI$）与 Monin-Obukhov 相似理论串联得到。水体/裸地分支（$veg\le1$）退化为 Penman 开阔水面方程；截留蒸发项 $pet_i$ 取气孔阻力为 0 的极限。

## 4. 关键变量与单位

| 符号 | 含义 | 单位 |
|---|---|---|
| `tempk` | 参考气温 | K |
| `rad`, `rshort` | 净辐射 / 短波辐射 | W/m² |
| `press`, `presshp` | 大气压力 | Pa / hPa |
| `qair` | 比湿 | kg/kg |
| `wind` | 风速 | m/s |
| `lai` | 叶面积指数 | m²/m² |
| `hveg` | 植被高度 | m |
| `veg` | 土地覆盖类型 (USGS 1-30) | - |
| `floodflag` | 洪水标志 (0/1) | - |
| `Δt` | 计算时间步长 | s |
| `pet_s/c/w/i` | 土壤/冠层/水面/截留蒸发 | mm |

## 5. 与 Fortran 对应

| Fortran 子程序 | Julia 函数 | 差异 |
|---|---|---|
| `PET_PriestyTaylor` (soilfluxes.f90) | `potevap_priestly_taylor` | Julia 用局部 `CP_LOCAL` 覆盖模块常量；时间步长耦合在 PT 中缺失 |
| `PET_PenmanMonteith` | `potevap_penman_monteith` | 参数顺序不同：Fortran 先传湿度，Julia 用比湿+压强计算 `e` |
| `PET_ShutteworthWallace` | `potevap_shutteworth_wallace` | Julia 返回 12 元组；Fortran 仅返回 4 个分量（pet_s/pet_c/pet_w/pet_i） |

## 6. 引用

- 函数定义：`src/Evapotranspiration.jl:L41-L56` (PT)、`L60-L100` (PM)、`L102-L200+` (SW)
- 生物物理参数表：`src/Evapotranspiration.jl:L11-L43`
- 测试断言：`test/test_evapotranspiration.jl:L9` `@test pet > 0.0`，`L25` `@test pet < 100.0`，`L48-L60` SW 元组解构
- 示例调用：`example/example_usage.jl:L43-L57`
- 文档参考：docs/土壤水运动_ASAP.typ（PET 与产流耦合章节）

## 7. 已知问题与备注

- ✅ 拼写错误已通过别名修正（2026-06-22 修复）：`src/Evapotranspiration.jl:7-L8` 同时 export 旧名 `potevap_shutteworth_wallace`（双 t）与正确拼写别名 `potevap_shuttleworth_wallace`；新代码推荐使用别名，对应单元测试 `test/test_evapotranspiration.jl` 的别名一致性断言。
- **唯一 module 包裹**：`src/Evapotranspiration.jl:L4-L6` 是 `src/` 中唯一显式 `module Evapotranspiration ... end`；其它文件为 `include` 顶层脚本。
- `BIOPARMS` 索引 18-21 与 3-6 重复（USDA/IGBP 双分类并存），需在 `veg_int = clamp(veg_int,1,31)` 后安全访问。
- Priestley-Taylor 默认步长假设为日，内部转换 `rad*24*3600*1e-6` 仅适合逐日调用。
- 未导出任何符号：`export` 列表为空（除 L7-L8 的 PET 函数外），需以 `ASAP.potevap_xxx` 或 `include` 形式访问。