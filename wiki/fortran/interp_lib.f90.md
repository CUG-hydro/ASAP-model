# interp_lib.f90

> 源文件：`fortran/interp_lib.f90`
> 行数：690
> Module 声明：`MODULE interp_lib`
> Julia 来源：**无对应**（Julia 端通过 `Interpolations.jl` / 闭包实现同等插值；ASAP Julia 主体代码未直接消费本模块）
> 测试：Fortran 侧无单元测试；功能与 RAMS v4.3.0.2 等价
> 状态：已摄取

## 1. 功能概述

本文件源自 **RAMS（Regional Atmospheric Modeling System）v4.3.0.2** 的大气模式插值工具库，迁移至 ASAP 后用于 ERA5 / ERA5-Land / MODIS LAI 等强迫数据从源网格到模式网格的水平与垂直重采样。提供 11 个 `CONTAINS` 子程序，覆盖：z ↔ z* 坐标系互转、线性插值、按对数气压插值、双线性/三次样条（BINOM）、重叠二次（GDTOST）、高程匹配（HTINT 家族）、六阶中心差分权重（AWTCMP）。所有子程序以**纯过程式接口**形式暴露，无模块级状态。

## 2. 子程序签名

### 2.1 坐标系互转

```fortran
SUBROUTINE TRNCL1(VCTRA, ZLEV, ZSURF, HTOP, VCTRB, VCTRC, VCTRD, NN)
!  Z* → Z  坐标重采样
SUBROUTINE TRNCL2(VCTRA, ZLEV, ZSURF, HTOP, VCTRB, VCTRC, VCTRD, NN)
!  Z → Z*  坐标重采样
```

### 2.2 通用插值

```fortran
SUBROUTINE INTRP(A, B, ZVAL, Z, NZP)
!  线性插值（自变量 z 单调）
SUBROUTINE INTRRAP(A, B, PLNVAL, PLN, NZP)
!  线性插值（自变量为 log p）
SUBROUTINE BINOM(X1, X2, X3, X4, Y1, Y2, Y3, Y4, XXX, YYY)
!  三次 Lagrange 插值（支持缺失数据掩码 1e30）
```

### 2.3 格点 → 站点（GDTOST 家族）

```fortran
SUBROUTINE GDTOST(A, IX, IY, STAX, STAY, STAVAL)
!  重叠二次（4×4 子网格 → 站点）
SUBROUTINE GDTOST2(A, IX, IY, STAX, STAY, STAVAL)
!  双线性（直接给 stax, stay 浮点位置）
SUBROUTINE GDTOST3(A, IX, IY, I, J, WT1, WT2, WT3, WT4, STAVAL)
!  双线性（直接给 i, j 与 4 个权重）
SUBROUTINE WEIGHTS(STAX, STAY, I, J, WT1, WT2, WT3, WT4)
!  计算双线性插值的 4 个权重
```

### 2.4 高程匹配（HTINT 家族）

```fortran
SUBROUTINE HTINTCP(NZZ1, VCTRA, ELEVA, NZZ2, VCTRB, ELEVB, &
                   VT1C, VT2C, VT1D, VT2D, VT1E, VT2E)
!  同时插 5 个变量 + 2 个辅助变量（VT1C/D/E → VT2C/D/E）
SUBROUTINE HTINT(NZZ1, VCTRA, ELEVA, NZZ2, VCTRB, ELEVB)
!  单变量高程插值（外推使用最近端线性外插）
SUBROUTINE HTINT2(NZZ1, VCTRA, ELEVA, NZZ2, VCTRB, ELEVB)
!  单变量高程插值（低于 ELEVA(1) 时保持 ELEVA(1) 值）
```

### 2.5 六阶中心差分权重

```fortran
SUBROUTINE AWTCMP(X, II, XX, IJ, ITYP, LLB, LRB, ICON, W, IORD)
!  计算六阶精度的有限差分权重（含 IORD=2 简化版）
```

## 3. 算法 / 公式

### 3.1 z* ↔ z 坐标系

定义 z* 地形追随坐标系（`ZSURF` 为地形高度，`HTOP` 为模式顶高度，`NN` 为层数）：

$$
r = 1 - \frac{Z_\text{SURF}}{H_\text{TOP}}, \quad \sigma_z(k) = r \cdot Z_\text{LEV}(k) + Z_\text{SURF}
$$

`TRNCL1`（z* → z）按 `Z_LEV(K)` 在 `σ_z` 数组中做线性反插；`TRNCL2`（z → z*）按 `σ_z(K)` 在 `Z_LEV` 数组中做线性插值。两者均使用 **一阶线性** 而非对数插值，适合气压、温度等连续场。

### 3.2 通用插值

```fortran
INTRP:   B = (A(K)*(ZVAL-Z(KABV-1)) + A(KABV-1)*(Z(KABV)-ZVAL)) / (Z(KABV)-Z(KABV-1))
INTRRAP: 同上，仅自变量由 Z 改为 PLN（log p）
```

`KABV` 是首个满足 `ZVAL < Z(K)` 的下标，若 `KABV==1` 则强制 `KABV=2`（保护下界）。

### 3.3 BINOM（带缺失值的三次 Lagrange）

`BINOM` 在 4 个节点上做**三次 Lagrange 插值**，但允许 `> 1e30` 表示"缺失"。完整分支：
- 若 `Y2` 或 `Y3` 缺失（`> 1e30`） → 直接 `YYY = 1e30`
- 若 4 点全部存在 → 标准三次 Lagrange（4 个 YZ1k 权重）
- 若端点 `Y1` 缺失 → 简化为 `YYY = WT1*Y2 + WT2*Y3`
- 若端点 `Y4` 缺失 → 简化为 `WT1*Y1 + WT2*Y2`

`WT1 = (XXX-X3)/(X2-X3)`、`WT2 = 1 - WT1` 是 X 方向的坐标插值权重；实际 YYY 综合 X 与 Y 两个方向的三次基函数：

$$
\Pi_k(X) = \prod_{j \neq k} \frac{X - X_j}{X_k - X_j},\quad k = 1,2,3,4
$$

### 3.4 GDTOST 家族

**GDTOST**（重叠二次）：

对 `STAX`/`STAY`（浮点网格坐标）取 4×4 子网格 `[IX1..IX2] × [IY1..IY2]`，逐列调用 `BINOM`（4 点三次）合成 4 个 `SCR`，再沿 X 方向调用 `BINOM` 得到 `STAVAL`。缺失值（`> 1e30`）由 `BINOM` 处理。

**GDTOST2 / GDTOST3**（双线性）：

$$
\text{STAVAL} = w_{tx,1}(w_{ty,1} a_{i,j} + w_{ty,2} a_{i,j+1}) + w_{tx,2}(w_{ty,1} a_{i+1,j} + w_{ty,2} a_{i+1,j+1})
$$

其中 `wtx2 = stax - int(stax)`，`wty2 = stay - int(stay)`；`wt_1 = 1 - wt_2`。`GDTOST3` 不再重新计算权重，直接接收 4 个权重 `WT1..WT4`。

### 3.5 HTINT 家族

按高层高程 `ELEVA(NZZ1)` 数组与目标高程 `ELEVB(NZZ2)` 数组做**单变量**线性插值。三种端点行为：

| 子程序 | 低于 ELEVA(1) | 高于 ELEVA(NZZ1) |
|---|---|---|
| `HTINT` | 视为 `ELEVB(K) < ELEVA(1)` → 用 `WT=(ELEVB(K)-ELEVA(L))/(ELEVA(L+1)-ELEVA(L))` | 端点外推：按 `WT=(ELEVB(K)-ELEVA(NZZ1))/(ELEVA(NZZ1-1)-ELEVA(NZZ1))` |
| `HTINT2` | **保持 ELEVA(1) 值不变**（用于"饱和层以上不再变化"的物理量） | 同 `HTINT` |
| `HTINTCP` | 同 `HTINT`，但同时插 5 个变量 | 同 `HTINT` |

### 3.6 AWTCMP（六阶中心差分权重）

`AWTCMP` 计算**六阶精度**有限差分权重矩阵 `W(6, 6)`，将 $\partial^{m} f / \partial x^{m}$（$m = 0 \mathrel{..} 5$）表示为 `X(II-3..II+3)` 处函数值的线性组合。

**分支条件**：
- `ITYP=0`（正常）或 `ITYP=1` 且 `II ≥ LLB+2`：
  - 若 `II-2 < LLB`、`II+3 > LRB`，或 `IORD=2`，使用 **简化 4 点线性权重**：
    $$
    W(1,3) = (X(II+1) - XX(IJ)) / (X(II+1) - X(II))
    $$
    $$
    W(2,3) = 0.5 / (X(II+1) - X(II))
    $$
    $W(1,4), W(2,4)$ 类似但符号反向
  - 否则使用 **六阶完整 6 点权重**（Lagrange 求导），含 $Q(L) = X(II-3+L), L=1..6$ 的 Lagrange 基函数展开

- `ITYP=1` 且 `II < LLB+2`（对称边界）：用对称镜像扩展：
  - `II == LLB`：`Q(2) = 2*X(II) - X(II+1)`、`Q(1) = 2*X(II) - X(II+2)`
  - `II == LLB+1`：`Q(2) = X(II-1)`、`Q(1) = 2*X(II) - X(II+1)`

- `ITYP=1` 且 `II < LLB+2` 的最终权重调整（L=4,5,2,1）确保在边界处保持对称性。

- `ICON=1`（退化分支）使用 **固定数值权重**（分母项 $W(L,LL)/(X(II+1)-X(II))^{LL-1}$）：

| (L, LL) | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| 1 | 1/60 | -8/60 | 37/60 | 37/60 | -8/60 | 1/60 |
| 2 | 2/360 | -25/360 | 245/360 | -245/360 | 25/360 | -2/360 |
| 3 | -1/48 | 7/48 | -6/48 | -6/48 | 7/48 | -1/48 |
| 4 | -1/144 | 11/144 | -28/144 | 28/144 | -11/144 | 1/144 |
| 5 | 1/240 | -3/240 | 2/240 | 2/240 | -3/240 | 1/240 |
| 6 | 1/720 | -5/720 | 10/720 | -10/720 | 5/720 | -1/720 |

最后除以 $(X(II+1)-X(II))^{L-1}$ 完成量纲归一。

## 4. 关键变量与单位

| 符号 | 含义 | 单位 |
|---|---|---|
| `VCTRA, VCTRB, VCTRC, VCTRD` | 输入 / 输出 / 倒数 / 坐标辅助向量 | 取决于上下文 |
| `ZLEV` | 高度坐标数组 | m |
| `ZSURF, HTOP` | 地表高度 / 模式顶高度 | m |
| `ZVAL, PLNVAL` | 目标高度 / 目标 log 气压值 | m 或 log hPa |
| `XXX, YYY` | BINOM 自变量 / 因变量 | - |
| `STAX, STAY` | 站点相对网格位置（浮点） | grid pt |
| `WT1..WT4` | 双线性插值的 4 个权重（满足 $w_{11}+w_{12}+w_{21}+w_{22}=1$） | - |
| `ELEVA, ELEVB` | 源 / 目标高程数组 | m |
| `W(L, LL)` | AWTCMP 输出：$m = LL-1$ 阶导数对应 $L$ 号节点的权重 | 单位由 $m$ 与 $X$ 单位共同决定 |
| `Q12` | `XX(IJ)`（目标坐标），用于六阶基函数展开 | 与 X 同单位 |
| `ITYP, ICON, IORD` | 边界类型（0=正常, 1=对称）/ 权重来源（0=六阶, 1=固定）/ 阶数（2=简化） | - |
| `LLB, LRB` | 网格左 / 右边界下标 | grid pt |
| 缺失值掩码 | `> 1.0e30`（`1E30`） | - |

## 5. 与 Julia 对应

| Fortran 子程序 | Julia 对应 | 差异 |
|---|---|---|
| `TRNCL1 / TRNCL2` | 无对应（Julia 不使用 z* 坐标系） | — |
| `INTRP` | `Interpolations.jl::linear_interpolation` 或内联闭包 | Julia 接口更通用；ASAP 未消费 |
| `INTRRAP` | 同上（自变量为 log p） | Julia 标准库无 `log pressure` 特化 |
| `BINOM` | `Interpolations.jl::Cubic`（Lagrange 等价） | Julia 缺失数据需用 `NaN` / `missing` 替代 `1e30` |
| `GDTOST` | `Interpolations.jl` + 4×4 局部采样 | Julia 标准库无"重叠二次"专门函数 |
| `GDTOST2 / GDTOST3` | 内联双线性公式 | Julia 主体代码未消费 |
| `WEIGHTS` | 内联权重计算 | 同上 |
| `HTINT / HTINT2 / HTINTCP` | Julia 主体代码未消费（ASAP 不用高程插值） | — |
| `AWTCMP` | `FiniteDiff.jl::derivative` 或 `DiffEqOperators` | ASAP Julia 端未直接使用 |
| 整体模块 `interp_lib` | `wiki/mapping/julia-fortran-对照.md` §3 显式标注 "Julia 端无对应" | Fortran 专用 |

## 6. 引用

- 行号：
  - `fortran/interp_lib.f90:L1-L13` 文件头注释 / Change Log / 版权
  - `fortran/interp_lib.f90:L17-L46` TRNCL1
  - `fortran/interp_lib.f90:L48-L77` TRNCL2
  - `fortran/interp_lib.f90:L83-L100` INTRP
  - `fortran/interp_lib.f90:L102-L121` INTRRAP
  - `fortran/interp_lib.f90:L123-L182` BINOM
  - `fortran/interp_lib.f90:L188-L223` GDTOST
  - `fortran/interp_lib.f90:L227-L255` GDTOST2
  - `fortran/interp_lib.f90:L259-L280` GDTOST3
  - `fortran/interp_lib.f90:L286-L312` WEIGHTS
  - `fortran/interp_lib.f90:L318-L347` HTINTCP
  - `fortran/interp_lib.f90:L353-L385` HTINT
  - `fortran/interp_lib.f90:L391-L420` HTINT2
  - `fortran/interp_lib.f90:L426-L669` AWTCMP（最复杂，含大量六阶 Lagrange 展开）
- 调用关系：
  - `module_forcings.f90::READFORCINGS` 内部使用 GDTOST2 / GDTOST3 / HTINT 系列做 ERA5 → 模式网格插值
  - `module_forcings.f90::READLAI` 使用 HTINTCP 同步插 5 个变量到目标高程
- 起源：RAMS v4.3.0.2（Mission Research Corporation / ASTeR Division, 1990–2000）
- 对照：`wiki/mapping/julia-fortran-对照.md §3（"网格插值（双线性 / 三次）"行）`

## 7. 已知问题与备注

1. **RAMS 版权头**：注释中标明"Copyright (C) 1990, 1995, 1999, 2000 - All Rights Reserved, Regional Atmospheric Modeling System - RAMS, Mission Research Corporation / ASTeR Division"；ASAP 仍保留此头但实际仅作为衍生工具库使用。
2. **`GDTOST` 与 `GDTOST2` 共存**：4.3.0.2 版本（010328 JHC）新增 `GDTOST2` 作为消除"ringing errors"的简化双线性方案；ASAP 主体使用 `GDTOST2/GDTOST3`，旧版 `GDTOST` 仅作历史保留。
3. **`HTINT` 在 `L == NZZ1` 时调用 `STOP`**：极端边界情形（`L == NZZ1`）会触发 `STOP 'htint'`，属于 Fortran 异常终止；调用前应保证 `ELEVB(K) ≤ ELEVA(NZZ1)`。
4. **`HTINT2` 在 `L == NZZ1` 时同样 `STOP`**：与 `HTINT` 一致，源 / 目标高程不匹配时会硬终止。
5. **`BINOM` 缺失值约定 `1E30`**：与 NetCDF 默认 `_FillValue` 不同（NetCDF 通常用 `9.96921e36`），需要 `READFORCINGS` 在读入时把 `_FillValue` 改写为 `1e30` 才能正确触发缺失分支。
6. **`AWTCMP` 的 `IORD=2` 路径仅算一阶导**：`W(2,3)` / `W(2,4)` 仅支持一阶导系数；当 `IORD=2` 时，`W(1,*)` 仍给出函数值（线性插值权重），高阶导数无意义。
7. **`AWTCMP` 缺少 `IMPLICIT NONE`**：早期 Fortran 子程序使用隐式类型，`I`/`J`/`K`/`L`/`II` 等默认为 `INTEGER`，调用者需自行确保一致。
8. **未在 Julia 端移植**：ASAP Julia 主体代码未直接消费本模块；`wiki/mapping/julia-fortran-对照.md §3` 显式标注"网格插值（双线性/三次）：Fortran 专用"。
9. **`TRNCL1` / `TRNCL2` 引用未声明的 `ID` / `JD`**：`WRITE(6,5) ID,JD` 与 `FORMAT(' STOP IN TRNCL1',2I5)` 中 `ID`/`JD` 是模块外全局变量（典型 RAMS 设计），ASAP 端编译时若无这两个全局变量会触发编译错误；本文件在 ASAP 编译环境中通常通过 `module_forcings` 的 include 与外部 `ID`/`JD` 提供，编译流程需谨慎。