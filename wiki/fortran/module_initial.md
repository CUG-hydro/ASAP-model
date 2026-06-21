# module_initial.f90

**路径**：`/mnt/z/GitHub/jl-pkgs/ASAP-model/fortran/module_initial.f90`
**行数**：1059
**Module 声明**：`MODULE module_initial`
**依赖**：`use module_parallel`、`use module_rootdepth`、`use module_wtable`、`use nrtype`

本模块负责模型初始化：调用 `LATERALFLOW` 计算稳态侧向流、构造初始土壤含水量、构造平衡地下水位以上的 `smoieq` 剖面，并提供 Brent 求根算法。

## 关键 Subroutine / Function 签名

### INITIALIZE

```fortran
SUBROUTINE INITIALIZE(n2, n3, is, ie, js, je, nzg, freedrain, slz, dz,
   soiltextures, wtd, smoi, smoieq, fdepth, topo, landmask, deltat, area)
```

主初始化入口。对每个 (i,j)：
1. 由 `wtd` 找到所在土壤层 `iwtd`；
2. `wtd` 以下填 `slmsts * exp((slz+1.5)/fdepth)`；
3. `wtd` 所在层按饱和与 `wgpmid` 线性插值；
4. `wtd` 以上各层用 `zbrent` 求平衡含水量（使稳态通量 = 0）。

若 `freedrain=0`，先调用 `LATERALFLOW` 计算稳态地下水。

### EQSOILMOISTUREtheor

```fortran
SUBROUTINE EQSOILMOISTUREtheor(is, ie, js, je, nzg, slz, dz,
   soiltextures, landmask, fdepth, smoieq)
```

**核心平衡含水量计算**：对地下水位以上的每一层，使用 `func(x, bexp, smoisat, hydcon, psisat, slz, vctr4, vctr6, dz, smoi1, flux)` 表示"重力通量 + 毛细通量 + 给定通量"的残差函数，通过 `ZBRENT` 求根（`x ∈ [0, 1]` 为归一化含水量），并夹紧到 `[soilcp, slmsts]`。

### EQSOILMOISTURE（数值法旧版）

```fortran
SUBROUTINE EQSOILMOISTURE(is, ie, js, je, nzg, slz, dz, dtll,
   soiltextures, landmask, smoieq)
```

通过对流扩散方程长时间积分至稳态（Newton 迭代），已弃用但保留。

### EQSOILMOISTUREiter

```fortran
SUBROUTINE EQSOILMOISTUREiter(is, ie, js, je, nzg, slz, dz, dtll,
   soiltxt, landmask, fdepth, smoieq)
```

迭代版（最大 500 步），适用于 `freedrain=0` 模式；调用 `SOILFLUXES_EQSMOI` 推进。

### SOILFLUXES_EQSMOI

```fortran
subroutine SOILFLUXES_EQSMOI(i, j, nzg, ztop, dtll, slz, dz, soiltxt,
   smoibotbc)
```

简化的 Richards 求解器（仅用于初始化），组装三对角矩阵并调用 `tridag`。

### ZBRENT（Brent 求根）

```fortran
FUNCTION zbrent(x1, x2, tol, bexp, smoisat, hydcon, psisat, slz, vctr4,
   vctr6, dz, smoi1, flux)
```

**Brent's method**：在 `[x1, x2]` 区间内对 `func` 求根，反二次插值 + 二分法保收敛。`ITMAX=100`，机器精度 `EPS=epsilon(x1)`，容差 `tol`。出错时调用 `nrerror` 终止。

### FUNC（被 ZBRENT 调用）

```fortran
FUNCTION func(x, bexp, smoisat, hydcon, psisat, slz, vctr4, vctr6, dz,
   smoi1, flux)
```

残差函数：返回 `d1*(x-smoi1)/dz + k1 + flux`，其中 `k1`、`d1` 为基于 Clapp-Hornberger 关系的水力传导度与扩散系数。

### ZBRENT1 / ZBRENT2 / ZBRENT3（变种）

针对不同情形（双层界面、三层界面）求根的 Brent 变种，分别用 `func1`、`func2`、`func3` 计算残差。

### NRERROR

```fortran
SUBROUTINE nrerror(string)
```

NR 系列统一错误处理：输出消息并 `STOP` 终止。

## 全局状态

- 无模块级变量，全部通过参数传递
- 通过 `use module_parallel` 共享 `pid`、`numtasks` 等 MPI 信息
- `use module_rootdepth` 获得 `slmsts`/`slcons`/`slpots`/`slbs`/`klatfactor` 等土壤参数
- `use module_wtable::LATERALFLOW` 用于稳态地下水计算
- `use nrtype` 提供 `SP`/`DP`/`I4B` 精度定义

## 数值细节

- `tol = 0.0001`（`INITIALIZE` 中设置）
- 含水量夹紧：`smoieq(k) = min(max(wmid, smoicp), smoisat)`
- 深度衰减：`smoisat = slmsts(nsoil) * max(min(exp((slz+1.5)/fdepth), 1.), 0.1)`

## 典型调用流程

```
READWTDNC → EQSOILMOISTUREtheor → INITIALIZE
                                      ├─ LATERALFLOW (freedrain=0)
                                      └─ ZBRENT (per cell per layer)
                                              └─ func (Clapp-Hornberger)
```
