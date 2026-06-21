# SoilParameters - 土壤参数

> 源文件：src/SoilParameters.jl:L1-L97
> Fortran 来源：fortran/module_initial.f90（参数数组定义）
> 测试：test/test_soil_parameters.jl
> 状态：已摄取

## 1. 功能概述
该模块定义了 ASAP 模型的土壤类型参数系统，包括 13 种土壤类型的物理水力属性（饱和含水量、残差含水量、基质势、饱和导水率、侧向流因子、Clapp–Hornberger b 参数），并提供 `get_soil_params` 查询函数、`init_soil_param` 凋萎点计算函数以及 Campbell 经验公式的导水率 `cal_K` 函数。模块以 `SoilType` struct 为核心，配合若干模块级常量数组以全局 const 形式存储查找表。

## 2. 函数签名
```julia
struct SoilType
    θ_sat::Float64
    θ_cp::Float64
    θ_wilt::Float64
    ψsat::Float64
    Ksat::Float64
    K_latfactor::Float64
    b::Float64
end

function get_soil_params(soil_type::Int)::SoilType
function init_soil_param(nzg::Int)::Tuple{Array{Float64,2},Array{Float64,1}}
function cal_K(θ::Float64, θ_sat::Float64, Ksat::Float64, b::Float64)::Float64
```

## 3. 算法 / 公式

### Campbell 导水率公式（cal_K）
$$K(\theta) = K_{sat} \left(\frac{\theta}{\theta_{sat}}\right)^{2b+3}$$

### 凋萎点含水量（init_soil_param）
$$\theta_{wilt} = \theta_{sat} \left(\frac{\psi_{sat}}{\psi_{wilt}}\right)^{1/b}$$
其中 $\psi_{wilt} = -153.0\ \text{m}$（POTWILT_LOCAL）。

## 4. 关键变量与单位
| 符号 | 含义 | 单位 |
|---|---|---|
| θ_sat | 饱和体积含水量 | m³/m³ |
| θ_cp | 残余/凋萎含水量下限 | m³/m³ |
| θ_wilt | 凋萎点含水量 | m³/m³ |
| ψsat | 饱和基质势 | m |
| Ksat | 饱和水力导水率 | m/s |
| K_latfactor | 侧向流比例因子（Klat/Ksat） | — |
| b | Clapp–Hornberger 经验指数 | — |
| NVTYP | 植被类型数量（=30） | — |
| NSTYP | 土壤类型数量（=13） | — |

## 5. 与 Fortran 对应
| Fortran 子程序 | Julia 函数 | 差异 |
|---|---|---|
| module_initial 中的 `θSAT/SOILCP/SLBS/KSAT/ΨSAT/KLATFACTOR` 全局数组 | `const θSAT/SOILCP/SLBS/KSAT/ΨSAT/KLATFACTOR` | Julia 使用 Unicode 标识符，Fortran 数组索引从 1，值一致 |
| `calK` / Campbell 导水率子程序 | `cal_K` | 参数顺序、公式完全一致；Julia 显式声明 Float64 类型 |
| 凋萎点初始化块 | `init_soil_param` | 公式相同；Julia 返回元组 `(fieldcp, θ_wilt)` |

## 6. 引用
- 常量定义：src/SoilParameters.jl:L8-L10（NVTYP、NSTYP）
- 查找表数组：src/SoilParameters.jl:L25-L31
- struct SoilType：src/SoilParameters.jl:L13-L21
- 函数 get_soil_params：src/SoilParameters.jl:L43-L57
- 函数 init_soil_param：src/SoilParameters.jl:L68-L79
- 函数 cal_K：src/SoilParameters.jl:L94-L96
- 测试断言：test/test_soil_parameters.jl:L6-L22（L26 边界异常、L51 饱和自洽、L65 公式验证）

## 7. 已知问题与备注
- `get_soil_params` 中 `θ_wilt` 字段写死为 0.0（L51），实际值需通过 `init_soil_param` 单独计算后回填，存在两阶段初始化耦合。
- `export SoilType, get_soil_params, init_soil_param, cal_K`（L6）在模块内 export，但 `ASAP.jl` 未将 SoilParameters 子模块 `@reexport`，外部需 `using ASAP` 直接可见；调用者更常通过 `get_soil_params(:)` 取值。
- `init_soil_param` 返回的 `fieldcp` 数组当前始终为 0（L71），未被填充，疑似占位实现。
- 文档：docs/土壤水运动_ASAP.typ 包含 Campbell 公式背景说明。
- 备注：DataFrames 在 SoilInitialization.jl 中被 `using`（L1）但当前实现并未使用（虽然 SoilInitialization 在 ASAP.jl L20 中被 include 后，DataFrames 仅作为传递依赖），属于悬空 import。
