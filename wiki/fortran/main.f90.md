# main.f90

> 源文件：`fortran/main.f90`
> 行数：58
> Module 声明：**(无 module 包裹的代码片段)**，备用时步调度脚本
> Julia 来源：`src/RootDepth.jl` 的 `rootdepth_main` 内部串联
> 状态：已摄取

## 1. 功能概述

本文件是 ASAP 模型的**备用单步调度脚本**：与 `module_driver.f90::program driver`（并行版主循环）平行存在，按 **LATERAL → GW2RIVER → ROOTDEPTH → FLOODING → RIVERS_KW_FLOOD** 顺序串行调度各子程序一次。文件全部由顶层调用代码构成，未定义 `module`、未 `use` 任何模块，依赖运行时上下文提供全部变量。其本质是把 `module_driver` 主循环的"一个时间步"剥离为可独立编译的诊断性入口，便于离线验证 `ROOTDEPTH` 与 `RIVERS_KW_FLOOD` 的输入输出。

## 2. 函数签名

`main.f90` **不定义任何子程序**，仅为顶层调用序列：

```fortran
! (片段示意)
call LATERAL(n2, n3, is, ie, js, je, nzg, soiltxt, wtd, qlat, fdepth, &
             topo, landmask, deltatwtd, area, lats, dxy, slz, &
             o18, smoi, qlato18, qlatinsum, qlatoutsum, qlatino18sum, &
             qlatouto18sum)
call GW2RIVER(n2, n3, is, ie, js, je, nzg, slz, deltat, soiltxt, &
              landmask, wtd, maxdepth, riverdepth, riverwidth, riverlength, &
              area, fdepth, qrf)
call ROOTDEPTH(freedrain, is, ie, js, je, nzg, slz, dz, deltat, landmask, &
               veg, hveg, soiltxt, wind, temp, qair, press, netrad, rshort, &
               lai, precip, qsrun, smoi, smoieq, smoiwtd, wtd, &
               waterdeficit, watext, watextdeep, rech, deeprech, &
               et_s, et_i, et_c, intercepstore, ppacum, pppendepth, &
               pppendepthold, qlat*deltat/deltatwtd, qslat, qsprings, &
               inactivedays, maxinactivedays, fieldcp, fdepth, steps, &
               floodheight, qrf, delsfcwat, icefactor, &
               wtdflux(is,js,daypast), et_s_daily(is,js,daypast), &
               et_c_daily(is,js,daypast), transptop(is,js,daypast), &
               infilk(is,js,daypast), infilflux, infilfluxday, &
               infilcounter, hour, o18, o18ratiopp, tempsfc, &
               qlato18*deltat/deltatwtd, transpo18, upflux)
call FLOODING(n2, n3, is, ie, js, je, deltat, fd, bfd, topoflood, area, &
              riverwidth, riverlength, riverdepth, floodheight, delsfcwat)
call RIVERS_KW_FLOOD(n2, n3, is, ie, js, je, deltat, dtlr, fd, bfd, &
                     riverflow, qsrun, qrf, delsfcwat, slope, riverdepth, &
                     riverwidth, riverlength, maxdepth, area, riverarea, &
                     floodarea, riverchannel, riverflowmean, floodheight, &
                     topoflood)
```

## 3. 算法 / 公式

无独立公式，仅按时序串联 `module_wtable` 中的水平水循环与 `module_rootdepth` 中的垂向水循环：

1. **LATERAL**：D8 侧向地下水流，返回 `qlat`、`qlato18`。
2. **GW2RIVER**（条件 `riverswitch == 1`）：地下水位超过河床底时计算基流 `qrf`。
3. **ROOTDEPTH**：核心调度，依次执行 PET → INTERCEPTION → EXTRACTION → SOILFLUXES → UPDATESHALLOWWTD，并累积 `et_s` / `et_c` / `transptop` 日值与 `rech` / `deeprech` / `upflux` 垂向通量。
4. **FLOODING**：D8 漫流检查 `topo + floodheight > topoflood`，交互 `delsfcwat`。
5. **RIVERS_KW_FLOOD**：运动波河道演算 + 漫滩交换，输出 `riverflowmean`。

**侧向通量缩放**：调用 ROOTDEPTH 时将 `qlat * deltat / deltatwtd`（即 `Δt_wtd` 步长内的累计通量折算到主时间步）与 `qlato18 * deltat / deltatwtd` 一并传入，保证 `Δt_wtd` ≠ `Δt` 时侧向流在 ROOTDEPTH 内部被正确归一化。

**累积约定**：

$$
\text{qlatsum} \mathrel{+}= \text{qlat} \cdot 10^{3} \quad [\text{mm}]
$$

$$
\text{qrfsum} \mathrel{+}= \text{qrf} \cdot 10^{3}, \quad \text{qsrunsum} \mathrel{+}= \text{qsrun} \cdot 10^{3}
$$

$$
\text{delsfcwatsum} \mathrel{-}= \text{delsfcwat} \cdot 10^{3}
$$

**日采样**：每 5 天一次（`mod(daysfromstart - 1, 5) == 0` 且 `hour == 0`）保存 `wtd_daily(:, :, daypast)` 与表层含水量 `smoi_daily(nzg, :, :, daypast)`，并推进 `daypast += 1`。

## 4. 关键变量与单位

| 符号 | 含义 | 单位 |
|---|---|---|
| `n2, n3` | 网格 I/J 维数（含 halo） | grid pt |
| `is, ie, js, je` | 子域范围 | grid pt |
| `nzg` | 土壤层数 | - |
| `deltat` | 主时间步（小时） | s |
| `deltatwtd` | 侧向流时间步（通常 = `deltat`） | s |
| `dtlr` | 河道子步长（5 min = 300 s） | s |
| `freedrain` | 自由排水标志 (0/1) | - |
| `riverswitch` | 河流开关（0 关闭 / 1 启用） | - |
| `slz, dz` | 层中心深度 / 层厚 | m |
| `topo, fdepth` | 地形 / 根系深度因子 | m |
| `wtd` | 地下水位深度（负值） | m |
| `smoi`, `smoieq`, `smoiwtd` | 当前 / 平衡 / 水位处含水量 | m³/m³ |
| `qlat, qrf, qsrun, qsprings` | 侧向流 / 基流 / 地表径流 / 泉流量 | m/s 或 m³/s |
| `qlat*deltat/deltatwtd` | 折算到主步的累计侧向通量 | m |
| `pet_s, et_s, et_i, et_c` | PET 与各分量 ET | mm 或 m |
| `rech, deeprech, upflux` | 补给 / 深层补给 / 上升流 | m |
| `icefactor` | 冰冻抑制因子 (Int8 数组) | 0/1 |
| `floodheight, delsfcwat` | 漫滩水深 / 地表水变化 | m |
| `riverdepth, riverwidth, riverlength` | 河道几何 | m |
| `slope, topoflood` | 河道坡度 / 漫滩高程 | m/m, m |
| `o18, o18ratiopp, precipo18, qlato18, transpo18` | ¹⁸O 同位素量 | ‰·m 或 m³/m³ |
| `et_s_daily`, `et_c_daily`, `transptop`, `infilk`, `wtd_daily`, `wtdflux` | 日均采样数组（按 `(n2, n3, daypast)` 维度） | mm 或 m |
| `t1, t2, t3` | MPI_WTIME 时戳（用于性能 profiling） | s |

## 5. 与 Julia 对应

| Fortran 调用 | Julia 对应 | 差异 |
|---|---|---|
| `LATERAL` | `lateral_flow!`（`src/modules/lateral_flow.jl:8-L58`） | Julia 8 邻居 + 角度因子；Fortran 旧版 5 邻居 |
| `GW2RIVER` | `gw2river!`（`src/modules/gw2river.jl:19-L72`） | Julia 类型参数化 `M<:Matrix{T}`；Fortran 显式形参 |
| `ROOTDEPTH` | `rootdepth_main`（`src/RootDepth.jl:90-L278`） | Julia 内部串联所有子模块；Fortran 单次调用传递全部形参 |
| `FLOODING` | `flooding!`（`src/modules/flooding.jl:6-L70`） | Julia `ntsplit=1` 固定；Fortran 可外部传入 |
| `RIVERS_KW_FLOOD` | `rivers_kw_flood!`（`src/modules/rivers_kw_flood.jl:11-L101`） | 接口对齐；Julia 用 `δt` Unicode 名 |
| 顶部脚本序列（`if riverswitch==1`…`end do`） | `rootdepth_main` 中按 `for itime ... end` 顺序编排 | Julia 隐藏调度细节，对外只暴露一次调用 |

> **注意**：Julia 没有与本文件逐行对应的入口；本文件实质是 `module_driver.f90::program driver` 主循环片段的独立副本，仅用于离线诊断。

## 6. 引用

- 行号：`fortran/main.f90:L1-L4` LATERAL 调用（带 `qlat*deltat/deltatwtd` 折算）；L7-L11 ROOTDEPTH 调用（最大形参数）；L20-L27 RIVERS_KW_FLOOD 调用
- 累计量维护：`fortran/main.f90:L5-L6, L13-L15`（`qlatsum`、`qrfsum`、`qsrunsum`、`delsfcwatsum`）
- 日采样：`fortran/main.f90:L17-L21`（每 5 天一次 `wtd_daily`、`smoi_daily` 写入）
- 关联调用：被 `module_driver.f90:381-L460` 主循环等价引用；可视为程序主循环的「单步外提」版
- 对照文档：`wiki/mapping/julia-fortran-对照.md §3（侧向/河流章节）` 与 `wiki/julia/RootDepth-主算法.md`

## 7. 已知问题与备注

1. **未声明 `program` 或 `module`**：文件全部为顶层执行语句，依赖调用上下文提供全部变量（隐式接口，无类型检查）。这是诊断用入口的常见做法，但存在 `IMPLICIT NONE` 缺失风险。
2. **`qlat * deltat / deltatwtd` 缩放**：当 `deltatwtd ≠ deltat` 时（极少见），侧向流会按比例放大；公式注释段（`if(freedrain.eq.0.and.hour.eq.nint(deltat/3600.)) then ... endif`）在主源码中被注释。
3. **冗余调用 `WTABLE`**：源代码中含被注释的 `! call WTABLE(...)`，说明历史上曾用 `WTABLE` 子程序完成整套水平流，目前已经解耦为 `LATERAL` + `GW2RIVER` + `FLOODING` + `RIVERS_KW_FLOOD` 四次显式调用。
4. **`freedrain.eq.0` 与 `nint(deltat/3600)` 的双条件**已被整段注释，意味着侧向流无条件每步调用，但 ROOTDEPTH 仍按 `freedrain` 标志切换边界条件。
5. **`riverswitch` 整型标志**：Julia 中没有对应（`rivers_kw_flood!` 与 `gw2river!` 永远启用），仅 Fortran 入口层开关。
6. **`wtd_daily(:, :, daypast) = wtd(:, :)`** 的赋值使用 `(:, :, daypast)` 三维形参展开，依赖外部 `wtd_daily` 预先按 `(n2, n3, 73)` 等维度分配好（每 5 天 1 次，最多 73 槽位覆盖 1 年）。
7. **无独立测试**：本文件未与 Julia 测试套件直接对应；其逻辑等价于 `RootDepth.jl` 主循环的一次时间步迭代，由 `test/runtests.jl` 间接覆盖。