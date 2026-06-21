# module_parallel.f90

**路径**：`/mnt/z/GitHub/jl-pkgs/ASAP-model/fortran/module_parallel.f90`
**行数**：658
**Module 声明**：`MODULE module_parallel`
**依赖**：`use mpi`（`include 'mpif.h'`）

本模块提供 MPI 并行框架：进程管理、域分解、halo 通信、自定义数据类型。`pid=0` 为 I/O 主进程，`pid=numtasks-1` 为收集/统计进程，中间 `pid=1..numtasks-2` 为计算进程。

## 全局参数与状态

```fortran
integer, parameter :: npmax = 1000    ! 最大进程数
integer, parameter :: n2big = 7320, n3big = 8520   ! 全球大网格
integer, parameter :: nw = 1, ne = n2big, ns = 1, nn = n3big
```

通过 `module_parallel` 共享的全局数组（在 `CONTAINS` 之前的模块体中定义）：

- `pid`、`numtasks`、`ierr`、`status(MPI_STATUS_SIZE)`
- `nini_x(0:npmax)`、`nend_x(0:npmax)`、`nini_y(0:npmax)`、`nend_y(0:npmax)`：每进程子域边界
- `n_right(0:npmax)`、`n_left(0:npmax)`、`n_up(0:npmax)`、`n_down(0:npmax)`：邻居进程号
- `domblock(0:npmax)`：子域 `MPI_Type_CONTIGUOUS(nmax_x*nmax_y, MPI_REAL, ...)`
- `domblock2dint` / `domblock3dint` / `domblocksmall` / `domblocksmallint` / `domblockbyte`：不同精度/字节版本
- `borderupdown(n)` / `borderleftright(n)`：单列/单行通信类型
- `updownhalo(n)` / `leftrighthalo(n)`：`MPI_TYPE_VECTOR` 提取的列/行
- `columntype`、`columntype2`：全局列类型
- `arraysection(n)` / `arraysectionint(n)`：`MPI_TYPE_CREATE_SUBARRAY` 描述每个子域在全球数组中的位置
- `rcountblocksmall(n)`、`rcountblock(n)`、`disp(n)`：通信计数与位移

## 关键 Subroutine 签名

### INITIALIZEDOMAIN

```fortran
SUBROUTINE INITIALIZEDOMAIN(n2, n3, nzg, filetopo)
```

域分解主入口：
1. `MPI_COMM_RANK` / `MPI_COMM_SIZE` 获得 `pid`、`numtasks`
2. 调用 `DIVIDEDOMAIN` 计算各进程子域边界与邻居关系
3. `pid=0` 通过 `MPI_send` 把分解结果广播到各 worker
4. 为每个子域创建自定义数据类型 `domblock`、`arraysection` 等（见上）
5. `MPI_TYPE_FREE` 释放临时 `tasktype`

### DIVIDEDOMAIN

```fortran
subroutine dividedomain(n2, n3, is, ie, js, je, n_right, n_left, n_up,
   n_down, nprocmax, filetopo)
```

按地形 `filetopo` 把全球划分为陆地块（跳过纯海洋），迭代缩小块大小直到陆地块数 `>= nprocmax-2`。最后给陆地块编号（`procid`）并填入 `nini_x/y`、`nend_x/y`。`pid=0` 与 `pid=numtasks-1` 保留给 I/O，不参与划分。

### DIVIDEDOMAINOLD（旧版）

```fortran
subroutine dividedomainold(n1, n2, ini, filetopo)
```

简单 1D 行向划分，已弃用但保留。

### SENDBORDERS4 / SENDBORDERS4blocking

```fortran
subroutine SENDBORDERS4(is, ie, js, je, wtd)
subroutine SENDBORDERS4blocking(is, ie, js, je, wtd)
```

**Halo 通信核心**：先用 `MPI_isend` + `MPI_irecv` 沿上下方向（tag 200/201）传递 `updownhalo` 类型列，再沿左右方向（tag 202/203）传递 `leftrighthalo` 类型行；阻塞版用 `MPI_send`/`MPI_recv`。先上下后左右以保证角点正确。

### SENDBORDERS / SENDBORDERSFLOOD（早期版本）

```fortran
subroutine SENDBORDERS(n2, js, je, wtd, reqsu, reqsd, reqru, reqrd)
subroutine SENDBORDERSFLOOD(n2, js, je, wtd, reqsu, reqsd, reqru, reqrd)
subroutine SENDBORDERSFLOOD4(is, ie, js, je, var, borderu, borderd,
   borderl, borderr)
```

早期一维列向 halo 通信，仅在南北方向（无角点处理）。`SENDBORDERSFLOOD4` 用于洪水模块，显式提供 `borderu/d/l/r` 缓冲。

## 典型调用顺序

```
program driver
  ├─ MPI_INIT
  ├─ INITIALIZEDOMAIN(n2,n3,nzg,filetopo)   ← 主入口
  │    └─ DIVIDEDOMAIN
  └─ 主循环中反复调用：
       SENDBORDERS4(is,ie,js,je,wtd)         ← 地下水位同步
       SENDBORDERSFLOOD4(is,ie,js,je,...)    ← 洪水同步
```

## 其他模块对其依赖

- `module_wtable::LATERAL` / `LATERALFLOW`：调 `SENDBORDERS4` 同步 `wtd`
- `module_rootdepth::ROOTDEPTH`：通过 halo 取邻居土壤水
- `module_forcings`：通过 `MPI_Type_CONTIGUOUS` + `MPI_ibcast` 分发强迫
- `module_io::WRITEOUTPUTNC_par` / `WRITEHISTORYNC_par`：用 `arraysection` 创建子域的 NetCDF 写入

## MPI 类型构造一览

| 类型 | 元素 | 用途 |
|------|------|------|
| `domblock(n)` | `nmax_x*nmax_y REAL` | 全子域实数通信 |
| `domblock2dint(n)` | `nmax_x*nmax_y INTEGER` | 全子域整数 |
| `domblock3dint(n)` | `2*nmax_x*nmax_y INTEGER` | 3D 整数 |
| `domblockbyte(n)` | `nmax_x*nmax_y BYTE` | 字节（`icefactor` 等） |
| `domblocksmall(n)` | `rcountblocksmall REAL` | 不含 halo 的核心区 |
| `updownhalo(n)` | `MPI_TYPE_VECTOR(nmax_x,1,1,...)` | 单列 halo |
| `leftrighthalo(n)` | `MPI_TYPE_VECTOR(nmax_y,1,nmax_x,...)` | 单行 halo |
| `arraysection(n)` | `MPI_TYPE_CREATE_SUBARRAY` | 全球数组中的子域视图 |
