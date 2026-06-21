# module_nrtype.f90

> 源文件：`fortran/module_nrtype.f90`
> 行数：36
> Module 声明：`MODULE nrtype`
> Julia 来源：**无对应**（Julia 依赖原生 `Float64`/`Int64` 与类型参数）
> 测试：Fortran 编译期即可校验
> 状态：已摄取

## 1. 功能概述

本文件提供 ASAP 模型各 Fortran 模块共享的**数值类型符号**与**常用数学常量**，并定义了 2 类稀疏矩阵派生类型。代码风格与 *Numerical Recipes* 系列一致（`I4B`/`SP`/`DP` 等命名），目的是把 `SELECTED_INT_KIND` / `SELECTED_REAL_KIND` 之类的 kind 值与可读的语义名称绑定，避免在主体代码中反复出现魔法数字。

## 2. 类型 / 参数签名

```fortran
MODULE nrtype
   INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)   ! 4 字节整数（kind=9）
   INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)   ! 2 字节整数（kind=4）
   INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)   ! 1 字节整数（kind=2）

   INTEGER, PARAMETER :: SP  = KIND(1.0)              ! 单精度实数
   INTEGER, PARAMETER :: DP  = KIND(1.0D0)            ! 双精度实数

   INTEGER, PARAMETER :: SPC = KIND((1.0, 1.0))       ! 单精度复数
   INTEGER, PARAMETER :: DPC = KIND((1.0D0, 1.0D0))   ! 双精度复数

   INTEGER, PARAMETER :: LGT = KIND(.true.)           ! 默认 logical

   REAL(SP), PARAMETER :: PI    = 3.141592653589793238462643383279502884197_sp
   REAL(SP), PARAMETER :: PIO2  = 1.57079632679489661923132169163975144209858_sp
   REAL(SP), PARAMETER :: TWOPI = 6.283185307179586476925286766559005768394_sp
   REAL(SP), PARAMETER :: SQRT2 = 1.41421356237309504880168872420969807856967_sp
   REAL(SP), PARAMETER :: EULER = 0.5772156649015328606065120900824024310422_sp

   REAL(DP), PARAMETER :: PI_D    = 3.141592653589793238462643383279502884197_dp
   REAL(DP), PARAMETER :: PIO2_D  = 1.57079632679489661923132169163975144209858_dp
   REAL(DP), PARAMETER :: TWOPI_D = 6.283185307179586476925286766559005768394_dp

   TYPE sprs2_sp
      INTEGER(I4B) :: n, len
      REAL(SP),    DIMENSION(:), POINTER :: val
      INTEGER(I4B),DIMENSION(:), POINTER :: irow
      INTEGER(I4B),DIMENSION(:), POINTER :: jcol
   END TYPE sprs2_sp

   TYPE sprs2_dp
      INTEGER(I4B) :: n, len
      REAL(DP),    DIMENSION(:), POINTER :: val
      INTEGER(I4B),DIMENSION(:), POINTER :: irow
      INTEGER(I4B),DIMENSION(:), POINTER :: jcol
   END TYPE sprs2_dp
END MODULE nrtype
```

## 3. 算法 / 公式

本模块**不含算法**，仅做 `KIND` 常量绑定与常量数值定义。

`PI` 等常量直接给出**长精度**（约 36 位有效数字），通过 `_sp` / `_dp` 后缀显式锚定精度，避免编译器将字面量截断到默认精度。

`EULER` 为欧拉-马歇罗尼常数 $\gamma \approx 0.5772$，其余常量均为 $\pi$ 及其派生。

`SELECTED_INT_KIND(R)` 返回能表示 $\pm 10^R$ 范围的最小整数 kind。`SELECTED_INT_KIND(9)` 通常对应 `INTEGER(4)`（gfortran/ifort 上等价于 `int32`），`SELECTED_INT_KIND(4)` 对应 `INTEGER(2)`，`SELECTED_INT_KIND(2)` 对应 `INTEGER(1)`。

## 4. 关键变量与单位

| 符号 | 含义 | 单位 |
|---|---|---|
| `I4B` | 4 字节整数 kind | - |
| `I2B` | 2 字节整数 kind | - |
| `I1B` | 1 字节整数 kind | - |
| `SP, DP` | 单/双精度实数 kind | - |
| `SPC, DPC` | 单/双精度复数 kind | - |
| `LGT` | 默认逻辑 kind | - |
| `PI, PIO2, TWOPI` | 圆周率 π、π/2、2π（SP） | - |
| `PI_D, PIO2_D, TWOPI_D` | 同上（DP） | - |
| `SQRT2` | $\sqrt{2}$（SP） | - |
| `EULER` | 欧拉常数 γ（SP） | - |
| `sprs2_sp / sprs2_dp` | 稀疏矩阵派生类型（COO 存储） | - |

## 5. 与 Julia 对应

| Fortran | Julia | 差异 |
|---|---|---|
| `I4B / I2B / I1B` | `Int64 / Int16 / Int8`（固定宽度） | Fortran 跨平台依赖 kind；Julia 标准宽度固定 |
| `SP / DP` | `Float32 / Float64` | Julia 标准宽度；ASAP 模型代码几乎全用 `Float64` |
| `SPC / DPC` | `ComplexF32 / ComplexF64` | Julia 不强制显式 kind |
| `LGT` | `Bool` | 同义 |
| `PI / TWOPI / SQRT2` | Julia 内置 `π / 2π / √2`（`Base.MathConstants`） | Julia 字面量由 IEEE 754 保证精度 |
| `EULER` | `Base.MathConstants.eulergamma` | 同义 |
| `sprs2_sp / sprs2_dp` | `SparseMatrixCSC{Float64,Int}` | Julia 用 CSC 压缩存储；Fortran 为稀疏三数组 (val, irow, jcol) COO |
| `module nrtype` 概念 | Julia 无对应（依赖原生类型系统） | Julia 通过 `const` / `@kwdef` 显式声明 |

## 6. 引用

- 行号：
  - `fortran/module_nrtype.f90:L1-L7` 整数 / 实数 / 复数 kind 定义
  - `fortran/module_nrtype.f90:L8-L9` logical kind
  - `fortran/module_nrtype.f90:L10-L20` 数学常量（SP + DP）
  - `fortran/module_nrtype.f90:L22-L36` 稀疏矩阵派生类型
- 调用关系：被 `module_rootdepth.f90`、`module_wtable.f90`、`module_io.f90`、`module_forcings.f90`、`module_initial.f90`、`module_parallel.f90` 等模块 `use nrtype`
- 命名参考：*Numerical Recipes* 系列（Press et al., 1992）的 `nrtype` 模块命名习惯
- 对照：`wiki/mapping/julia-fortran-对照.md §0`（明确"数据类型：module_nrtype.f90（无对应）"）

## 7. 已知问题与备注

1. **未在 ASAP 主体代码中广泛使用**：`sprs2_sp` / `sprs2_dp` 类型在 ASAP 实际代码中未被实例化，仅作为 *Numerical Recipes* 风格的占位继承；当前 ASAP 使用 `nf90_*` 与原生数组。
2. **`SELECTED_INT_KIND` 平台差异**：`SELECTED_INT_KIND(9)` 在 gfortran/ifort 上等价于 `INTEGER(4)`（int32），但 gfortran 上 `INTEGER(2)` 与 `INTEGER(1)` 在跨平台 ABI（特别是数组边界）有差异；ASAP 仅使用 `I4B`。
3. **`PI` 常量被 module_rootdepth 内部覆盖**：例如 `pi4=3.1415927*4.` 在 `module_wtable.f90` 中手动定义（约 7 位有效数字），而 `nrtype::PI` 长精度版未被引用。
4. **无独立测试**：本模块在 Fortran 编译期通过 `USE nrtype` 即可验证类型绑定，无需运行期断言。
5. **Julia 风格建议**：若要在 Julia 端对应地集中常量，可在 `src/ASAP.jl` 顶部声明 `const PI = π`、`const SQRT2 = √2` 等，无需新建模块。