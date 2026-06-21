# module_wtable.f90

**路径**：`/mnt/z/GitHub/jl-pkgs/ASAP-model/fortran/module_wtable.f90`
**行数**：1417
**Module 声明**：`MODULE module_wtable`
**依赖**：`use module_parallel`、`use module_rootdepth`

本模块负责地下水-河流-洪水的水平水循环：地下水侧向流、河道演进（Kinematic Wave / Diffusive Wave）、漫滩（floodplain）蓄水与交换。`pi4=3.1415927*4.` 为模块参数。

## 关键 Subroutine 签名

### WTABLE（第 13 行）

```fortran
subroutine WTABLE(imax, jmax, is, ie, js, je, nzg, slz, dz, area,
   soiltxt, wtd, bottomflux, rech, qslat, fdepth, topo, landmask, deltat,
   smoi, smoieq, smoiwtd, qsprings)
```

旧版地下水位主入口：调用 `SENDBORDERS` 同步 halo → 计算 `klat` → `lateralflow` → 更新深层水位 `UPDATEDEEPWTABLE` → 分配 `qsprings`。

### LATERAL（第 138 行）

```fortran
subroutine LATERAL(imax, jmax, is, ie, js, je, nzg, soiltxt, wtd, qlat,
   fdepth, topo, landmask, deltat, area, lats, dxy, slz, o18, smoi,
   qlato18, qlatinsum, qlatoutsum, qlatino18sum, qlatouto18sum)
```

主程序默认入口：使用 `LATERALFLOW` 子例程（带地形坡度 `topo`、`fdepth` 衰减），并把 O18 通量 `qlato18` 同步分配给各层 `smoi`。

### UPDATEDEEPWTABLE（第 185 行）

```fortran
subroutine UPDATEDEEPWTABLE(imax, jmax, js, je, nzg, slz, dz, soiltxt,
   wtd, bottomflux, rech, ...)
```

更新深层地下水位（自由排水模式下的等效底边界），由 `bottomflux` 与 `rech` 决定。

### LATERALFLOW（第 269 行）

```fortran
subroutine LATERALFLOW(imax, jmax, is, ie, js, je, wtd, qlat, fdepth,
   topo, landmask, deltat, area, klat)
```

核心水平流计算：用 `klat = slcons(nsoil) * klatfactor(nsoil)` 和水头差 `dH = wtd - wtd(neighbor)`，按 `dH * klat * dxy` 累加至 `qlat`。调用 `SENDBORDERS4` 同步邻居格点。

### LATERALFLOW4（被 LATERALFLOW 调用）

四向（上下左右）侧向流显式计算，对每个像元判断地形坡向与水力梯度。

### UPDATEWTD（第 468 行）

```fortran
subroutine UPDATEWTD(nzg, slz, dz, wtd, qspring, totwater, smoi,
   smoieq, soiltextures, smoiwtd)
```

给定总水量变化 `totwater = qlat + qspring + rech`，通过饱和带厚度与 `slmsts` 计算新水位 `wtd`，并修正 `smoi`（保持 `smoieq` 平衡剖面）。

### GW2RIVER（第 744 行）

```fortran
subroutine GW2RIVER(imax, jmax, is, ie, js, je, nzg, slz, deltat,
   soiltxt, landmask, wtd, maxdepth, riverdepth, riverwidth, riverlength,
   area, fdepth, qrf)
```

地下水补给河道：当地下水位超过河道底部 `riverdepth` 时，按 `K_sat * (wtd-riverdepth) * area` 估算 `qrf`（基流）。

### RIVERS_KW_FLOOD（第 800 行）

```fortran
subroutine RIVERS_KW_FLOOD(imax, jmax, is, ie, js, je, deltat, dtlr,
   fd, bfd, qnew, qs, qrf, delsfcwat, slope, riverdepth, riverwidth,
   riverlength, maxdepth, area, riverarea, floodarea, riverchannel,
   riverflowmean, floodheight, topoflood)
```

**Kinematic Wave** 河道演算 + 漫滩：使用流向 `fd` / 反向 `bfd`、河长 `riverlength`、河宽 `riverwidth`、坡度 `slope`，求解显式 KW 方程；与漫滩 `floodheight` 通过临界水深交换。`dtlr` 为河道子步长（5 min）。

### RIVERS_DW_FLOOD（第 996 行）

```fortran
subroutine RIVERS_DW_FLOOD(imax, js, je, deltat, dtlr, fd, bfd, qnew,
   qs, qrf, delsfcwat, ...)
```

**Diffusive Wave** 河道演算（更稳定但更耗资源），保留压力项。

### FLOODING（第 1226 行）

```fortran
subroutine FLOODING(imax, jmax, is, ie, js, je, deltat, fd, bfd,
   topoflood, area, riverwidth, riverlength, riverdepth, floodheight,
   delsfcwat)
```

漫滩淹没判别：当 `topo + floodheight` 超过 `topoflood` 时，交换河水与漫滩水量，更新 `delsfcwat`（地表水变化）。

### MOVEQRF（第 1325 行）

```fortran
subroutine MOVEQRF(imax, js, je, fd, qrf, area, width)
```

将本地产生的 `qrf`（地下水补给）按流向 `fd` 推送到下游格点。

### FLOWDIR（第 1384 行）

```fortran
subroutine FLOWDIR(is, ie, js, je, fd, ii, jj, i, j)
```

流向工具子例程：根据 `fd` 取得下一格点坐标 `(ii,jj)`。

## 全局状态

- 通过 `module_parallel` 共享进程网格 `pid`、`numtasks`、`domblock` 等
- 不含全局变量，全部通过参数传递
- 输入均以 `(is:ie, js:je)` 子域形式操作

## 与其他模块的接口

- 调用 `module_rootdepth::UPDATESHALLOWWTD` 和 `UPDATEWTDQLAT` 更新局部土壤水
- 调用 `module_parallel::SENDBORDERS` / `SENDBORDERS4` 同步 halo
- 被 `module_driver::program driver` 直接调用：`LATERAL`、`GW2RIVER`、`FLOODING`、`RIVERS_KW_FLOOD`、`RIVERS_DW_FLOOD`、`MOVEQRF`
