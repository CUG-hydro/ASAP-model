# module_driver.f90

**路径**：`/mnt/z/GitHub/jl-pkgs/ASAP-model/fortran/module_driver.f90`
**行数**：1118
**Module 声明**：`program driver`（主程序，**非** module）
**依赖**：`use module_parallel`、`use module_forcings`、`use module_rootdepth`、`use module_io`、`use module_wtable`、`use module_initial`

主程序负责：MPI 初始化、域分解、静态场读取、初始化、时步循环、日/月输出、关闭。

## 全局参数

```fortran
integer, parameter :: n2=7320, n3=8520, nzg=40
integer, parameter :: restart=1, freedrain=0, riverswitch=1
integer, parameter :: nvar_out=15, writepar=1
real, parameter :: deltat=3600., deltatwtd=3600., deltatriver=300.
real, parameter :: dxy=1./120., steps=1.
integer, parameter :: maxinactivedays=365*24*nint(steps)
```

`n2big=7320, n3big=8520`（全球大网格，由 `module_parallel` 共享）。

## 主要执行流程

1. **MPI 初始化**：`call MPI_INIT(ierr)` → `INITIALIZEDOMAIN(n2,n3,nzg,filetopo)`（域分解）。
2. **子域分配**：`is:nini_x(pid)`、`ie:nend_x(pid)`、`js:nini_y(pid)`、`je:nend_y(pid)`。
3. **静态场读取**（来自 `module_io`）：
   - `READINITIAL` → 土壤质地 `soiltxt`、地形 `topo`、`fdepth`、`landmask`
   - `READVEG` / `READHVEG` → 植被类型与高度
   - `READLATLON` → 经纬度与像元面积
   - `READTOPOERA5` → ERA5 高程
   - 河流：`READFLOWDIRECTION` / `READRIVERPARAMETERS` × 6 次
4. **土壤层定义**：`INITIALIZESOILDEPTH(nzg, slz, dz)` 给出 40 层深度。
5. **初始化分支**：
   - `restart=1`：`READHISTORYNC` + `EQSOILMOISTUREtheor` 恢复历史
   - `restart=0`：`READWTDNC` + `EQSOILMOISTUREtheor` + `INITIALIZE` + 读 O18
6. **时步循环** `DO WHILE (year .ne. 2024)`：
   - `READFORCINGSSOILT`（每小时前的土壤温度）
   - 时间推进 `hour/day/month/year`
   - `READFORCINGS`（小时 ERA5）+ `READFORCINGSSNOW`（融雪）
   - `READFORCINGSACC` 累积辐射/降水
   - `LATERAL` → `GW2RIVER` → `ROOTDEPTH` → `FLOODING` → `RIVERS_KW_FLOOD` × `niter_river`
   - 月初（`day.eq.1 .and. hour.eq.0`）：写日输出 + 入渗输出
   - 月初同条件：写 history 文件
7. **关闭**：`call MPI_FINALIZE(ierr)`

## 关键计时

`MPI_WTIME()` 取 5 个时戳 `t1..t5` 报告各段耗时：I/O、侧向流、ROOTDEPTH、河道。通过 `MPI_REDUCE` 取最大耗时与对应 `pid` 输出。

## 输出文件

- `rootdaily_wt_*.nc`：日均输出（15 变量：土壤水、ET、rech、qrf 等）
- `infiloutput_*.nc`：入渗通量与 O18 比例
- `history_wt_*.nc`：月初始场（用于 restart）
- `pathoutput = '/mnt/lustre/.../SAMERICA/outputparera5/'`（可改）

## 主要子例程调用

`MPI_INIT` → `INITIALIZEDOMAIN` → `READINITIAL/READVEG/READHVEG/READLATLON/READTOPOERA5` → `READFLOWDIRECTION/READRIVERPARAMETERS` → `INITIALIZESOILDEPTH` → `READHISTORYNC/EQSOILMOISTUREtheor` 或 `READWTDNC/INITIALIZE` → `READFORCINGS/READFORCINGSSNOW/READLAICHINA/READO18CLIM` → 循环 `LATERAL/GW2RIVER/ROOTDEPTH/FLOODING/RIVERS_KW_FLOOD` → `writeoutputnc_par/writehistorync_par`。
