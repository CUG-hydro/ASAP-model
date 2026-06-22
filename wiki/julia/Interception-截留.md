# 植被截留 (Interception)

> 源文件：src/Interception.jl:L1-L40
> Fortran 来源：fortran/soilfluxes.f90（拦截子程序 / 模块 ASAP7）
> 测试：test/test_interception.jl
> 状态：已摄取

## 1. 功能概述

本模块实现冠层降水截留的**水量平衡式更新**：依据当前截留库容 `intercepstore`、降水强度 `precip`、最小降水率阈值 `minpprate` 以及截留蒸发潜力 `pet_i`，将一次降水事件划分为穿透降水（throughfall）、截留蒸发与新的截留存储。截留容量由 `0.2 * LAI` 简单参数化，适合与 Shuttleworth-Wallace 输出的 `pet_i` 直接耦合，是 ASAP-model 陆面水循环中"冠层界面"环节。

## 2. 函数签名

```julia
function interception(precip::Float64, lai::Float64,
                     intercepstore::Float64, pet_i::Float64) ::
    Tuple{Float64, Float64, Float64}
```

返回 `(ppdrip, et_i, new_intercepstore)`，单位均为 mm：穿透降水量、截留蒸发量、更新后的冠层截留库容。

> 注意：原 Fortran 与 Julia 当前实现中存在未传入的参数 `minpprate`（最小降水率阈值，mm/timestep），文档字符串保留该参数说明，但函数体依赖的全局/外生约定需调用方确保。

## 3. 算法 / 公式

**截留容量**（LAI 线性参数化）：

$$S_{max} = 0.2 \cdot LAI \quad [\mathrm{mm}]$$

**库容不足量**：$D = S_{max} - S_t$

**(a) 降水量超过不足量（$P > D$）**

- 若 $P \geq P_{min}$：冠层饱和，无蒸发损失 $\Rightarrow E_i = 0$，穿透 $P_d = P - D$，更新 $S_{t+1} = S_{max}$
- 若 $P < P_{min}$：允许蒸发 $E_i = \min(S_{max}, PET_i)$，$S_{t+1} = S_{max} - E_i$

**(b) 降水量不超过不足量（$P \le D$）**

全部降水暂存于冠层：$P_d = 0$。当 $P < P_{min}$ 时按 $\min(S_t + P, PET_i)$ 蒸发；否则 $E_i = 0$。

$$S_{t+1} = S_t + P - E_i$$

## 4. 关键变量与单位

| 符号            | 含义                         | 单位        |
| --------------- | ---------------------------- | ----------- |
| `precip`        | 单步降水量                   | mm          |
| `lai`           | 叶面积指数                   | m²/m²       |
| `intercepstore` | 当前冠层截留量 $S_t$         | mm          |
| `pet_i`         | 截留蒸发潜力（来自 SW 模型） | mm          |
| `minpprate`     | 最小降水率阈值 $P_{min}$     | mm/timestep |
| `intercepmax`   | 截留容量 $S_{max}=0.2 LAI$   | mm          |
| `ppdrip`        | 穿透降水 $P_d$               | mm          |
| `et_i`          | 截留蒸发 $E_i$               | mm          |

## 5. 与 Fortran 对应

| Fortran 子程序                   | Julia 函数     | 差异                                                                                                                                     |
| -------------------------------- | -------------- | ---------------------------------------------------------------------------------------------------------------------------------------- |
| `INTERCEPTION`（soilfluxes.f90） | `interception` | Fortran 在子程序内通过 `MINPPRATE` 全局参数判定阈值；Julia 用 `const minpprate = 0.01` 模块常量（L3）+ 函数体 L33/L44 显式消费，接口对齐 |

## 6. 引用

- 函数定义：`src/Interception.jl:L13-L40`
- `intercepmax = 0.2 * lai` 公式：`src/Interception.jl:L15`
- 阈值分支：`src/Interception.jl:L19-L30`（PP > D）、`L31-L38`（PP ≤ D）
- 测试断言：`test/test_interception.jl:L20` `@test ppdrip ≈ precip - deficit`，`L42` `@test new_intercepstore ≈ intercepmax`
- 示例调用：`example/example_usage.jl:L62-L75`

## 7. 已知问题与备注

- **接口缺失**：`minpprate` 在文档字符串中作为参数列出（`src/Interception.jl:L5`），但函数签名与函数体均未引用该变量；调用方必须使用包内同名的全局常量或自行维护此阈值，否则 `if precip < minpprate` 会触发 `UndefVarError`。
- **裸地分支**：`lai = 0.0` 时 `intercepmax = 0.0`，所有降水穿透（`ppdrip = precip`），测试 `test_interception.jl:L75-L90` 验证此边界。
- **数值范围**：当前实现未对 `et_i` 做负值钳制，理论上极端参数组合可能返回极小负值，建议上游调用方加 `max(0, ...)` 防御。
- **未导出**：`export interception` 仅在行内注释中提及，实际模块无 `export` 语句，调用须显式 `ASAP.interception`。
- **耦合假设**：本模块假设 `pet_i` 由 Shuttleworth-Wallace 提供（截留蒸发项），单独使用时需传入经验值。
