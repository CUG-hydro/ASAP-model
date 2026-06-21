# module_io.f90

**路径**：`/mnt/z/GitHub/jl-pkgs/ASAP-model/fortran/module_io.f90`
**行数**：3599
**Module 声明**：`MODULE module_io`
**依赖**：`use netcdf`（NetCDF Fortran 库）；不显式 `use` 其他 ASAP 模块

本模块是 ASAP 模型的全部 NetCDF / 二进制 I/O 入口。包含 24 个 subroutine，按用途分为：**静态场读取**、**初始/历史读取**、**日/月输出**、**history 写入**、**工具**。

## 静态场读取

| Subroutine | 行号 | 作用 |
|------------|------|------|
| `READINITIAL` | 10 | 读土壤质地（二进制直读）、地形、`fdepth`、`landmask` |
| `READLATLON` | 212 | 读纬度/经度/像元面积（`dxy`、`swlat`、`swlon`） |
| `READWTD` / `READWTDNC` | 255 / 309 | 读初始地下水位（二进制 / NetCDF） |
| `READVEG` / `READHVEG` | 378 / 447 | 读 MODIS 植被类型与高度 |
| `READSMOIEQ` | 517 | 读平衡土壤含水量 |
| `READFLOWDIRECTION` | 591 | 读河流流向 `fd` 与反向 `bfd` |
| `READRIVERPARAMETERS` | 693 | 通用读取河网属性（按 `irec` 选择字段） |
| `handle_err` | 747 | NetCDF 错误处理 |

## 历史 / 初始读取

| Subroutine | 行号 | 作用 |
|------------|------|------|
| `READHISTORY` | 758 | 二进制历史读取 |
| `READHISTORYNC` | 823 | NetCDF 历史读取（smoi、wtd、intercepstore、inactivedays） |
| `READHISTORYVARNC` | 1050 | 单变量历史读取（`RIVERFLOW` 等） |
| `READHISTORYBYTEVARNC` | 1112 | 字节变量历史读取 |
| `READHISTORYVAR3DNC` | 1173 | 三维变量历史读取（`O18`） |
| `READHISTORYNCblock` | 1245 | 块状历史读取（restart 加速用） |

## 日 / 月输出

```fortran
subroutine WRITEOUTPUTFD(n2, n3, is, ie, js, je, nzg, smoi,
   waterdeficit, watext, qsrun, rech, pet, ppacum, filename, irec,
   istart, req, req2, req3)                       ! line 1461: free-drain GrADS
subroutine WRITEOUTPUTNC(n2, n3, is, ie, js, je, nzg, dz, smoi, ...)  ! line 1682: serial NC
subroutine WRITEOUTPUTNC_par(...)                 ! line 2115: parallel NC
subroutine WRITEOUTPUTNC_DAILY_par(...)           ! line 3095: 5-day diagnostic
subroutine WRITEOUTPUTNC_INFIL_par(...)           ! line 3311: infiltration + O18
```

并行版使用 `MPI_Type_CREATE_SUBARRAY`（由 `module_parallel` 提供 `arraysection`）实现每进程独立 NetCDF 文件写入。

## History 写入

```fortran
subroutine WRITEHISTORYNC_par(...)                ! line 2415
subroutine WRITEHISTORYNC(...)                    ! line 2661
subroutine WRITEEQSMOINC(...)                     ! line 3007
```

History 文件包含 `smoi`、`smoiwtd`、`intercepstore`、`wtd`、`inactivedays`、`riverflow`、`riverdepth`、`floodheight`、`o18` 等全部 restart 所需变量。

## 关键变量（输出维度）

- `nvar_out=15`（默认输出变量数，命令行参数控制）
- 输出维度 `(is:ie, js:je, nzg)` 三维，或 `(is:ie, js:je)` 二维
- 文件命名格式：`rootdaily_wt_<day><mon><year>_<pid>.nc`

## 关键依赖

- `module_parallel::arraysection` / `domblock` / `domblock2dint` 用于并行 I/O 的子数组类型
- `netcdf` 库（`nf90_open` / `nf90_get_var` / `nf90_put_var` / `nf90_def_var` 等）

## 典型调用顺序

**初始化**：`READINITIAL` → `READVEG` → `READHVEG` → `READLATLON` → `READTOPOERA5`（在 `module_forcings`）→ `READFLOWDIRECTION` → `READRIVERPARAMETERS` × 6 → `READWTDNC`（如 restart=0）。

**Restart**：`READHISTORYNC` → `READHISTORYVARNC` × 3（riverflow/depth/floodheight）→ `READHISTORYVAR3DNC`（o18）。

**输出**：每月 `WRITEHISTORYNC_par` + 每日 `WRITEOUTPUTNC_par` + `WRITEOUTPUTNC_INFIL_par`。
