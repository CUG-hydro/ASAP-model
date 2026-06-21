# module_forcings.f90

**路径**：`/mnt/z/GitHub/jl-pkgs/ASAP-model/fortran/module_forcings.f90`
**行数**：1615
**Module 声明**：`MODULE module_forcings`
**依赖**：`use module_parallel`、`use interp_lib`、`use netcdf`

本模块负责读取并分发 ERA5 / ERA5-Land / MODIS-LAI 强迫场，并将全球网格数据通过 `MPI_Type_CONTIGUOUS` 自定义类型广播到各子域，最后用 `gdtost2` / `gdtost3` 双线性插值到模式子域。

## 全局参数（ERA5 网格）

```fortran
integer, parameter :: dimerax=250, dimeray=290
integer, parameter :: xdim=1440, ydim=721
integer, parameter :: xdimeraland=3600, ydimeraland=1801
integer, parameter :: dimeralandx=620, dimeralandy=720
real, parameter :: eraswlat=-89.46282, eraswlon=0., dxera=0.703125
```

## 主要 Subroutine

### READFORCINGS（核心小时强迫）

```fortran
SUBROUTINE READFORCINGS(is, ie, js, je, filename, irec, hour, topo,
   lats, lons, vels, temp, pres, qair, tempsfc, landmask, topoera5,
   varpack, request, istart)
```

由 `pid=0` 从 `/mnt/lustre/.../ERA5/...` 读取：风（ws10）、气温（t2m）、气压（sp）、露点（d2m）、降水（tp+sf）、短波（ssrd）、净辐射（sr）、4 层土壤温度（stl1..stl4）共 11 个变量，组装到 `varpack(dimerax, dimeray, 0:23, 11)`，再 `MPI_ibcast` 分发。各 worker 进程做高度订正：`temp = tempint - 0.0065*(topo-topoera5)`。

### READFORCINGSACC

```fortran
SUBROUTINE READFORCINGSACC(is, ie, js, je, filename, irec, hour, lats,
   lons, pcpgl, rshort, netrad, landmask, varpack, req, istart)
```

累计式：读取总降水量 `pcpgl`、短波 `rshort`、净辐射 `netrad` 累积值，做差分后插值。

### READFORCINGSRNET

```fortran
SUBROUTINE READFORCINGSRNET(is, ie, js, je, filename, irec, hour,
   lats, lons, rnet, landmask, netrad, req, istart)
```

单独读取净辐射（短波 `ssr` + 长波 `str`），用 `MPI_isend`（tag 558）逐像元发送。

### READFORCINGSSNOW

```fortran
SUBROUTINE READFORCINGSSNOW(is, ie, js, je, filename, irec, hour,
   lats, lons, snow, landmask, snowera, request, istart)
```

读取 ERA5-Land `smlt` 融雪量，0.1° 高分辨率（`dimeralandx=620, dimeralandy=720`）。

### READFORCINGSSOILT

```fortran
SUBROUTINE READFORCINGSSOILT(is, ie, js, je, nzg, filename, irec,
   hour, lats, lons, icefactor, landmask, topo, topoera, smoipack,
   req, istart)
```

读取 4 层 ERA5 土壤温度，构造 `icefactor`（冰冻指示符）。

### READTOPOERA5

```fortran
subroutine READTOPOERA5(n2, n3, is, ie, js, je, landmask, lats, lons,
   topoera5)
```

读取 ERA5 静态地形（用于高度订正），tag 8883。

### READLAI / READLAICLIM / READLAICHINA

```fortran
subroutine READLAI(n2, is, ie, js, je, lats, lons, year, month, day, lai)
subroutine READLAICLIM(n2, n3, is, ie, js, je, month, lai)
subroutine READLAICHINA(n2, n3, is, ie, js, je, filelai, lai)
```

LAI 读取三套方案：MODIS 8 天产品、气候态、北京师范大学 30s 高分辨率（实际使用）。

### READO18CLIM

```fortran
subroutine READO18CLIM(n2, n3, is, ie, js, je, month, fileO18, o18ratio)
```

读取 O18 降水同位素气候态（半月份），按月线性插值到日。

## 工具函数

- `julday(imonth, iday, iyear)`：日期→年内序日
- `daynumber(day, month, year)`：日期→自 2003 起的天数
- `handle_err(statusnc)` / `check(status)`：NetCDF 错误处理
- `rslf(p, t)`：饱和水汽混合比（Tetens 公式）
- `READFORCINGS6h` / `READFORCINGS3h`：历史 ERA-Interim 兼容版

## 插值

通过 `interp_lib::gdtost2`（双线性）或 `gdtost3`（三次样条）从 ERA5 子域（250×290，南美区域）插值到模式子域 `(is:ie, js:je)`。
