# ASAP 模型 Julia 版本

<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://CUG-hydro.github.io/ASAP-model/stable) -->
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://CUG-hydro.github.io/ASAP-model/dev)
[![CI](https://github.com/CUG-hydro/ASAP-model/actions/workflows/CI.yml/badge.svg)](https://github.com/CUG-hydro/ASAP-model/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/CUG-hydro/ASAP-model/branch/main/graph/badge.svg)](https://app.codecov.io/gh/CUG-hydro/ASAP-model/tree/main)

这是 ASAP (Agricultural Systems Analysis and Prediction) 模型的 Julia 实现版本，从原始的 Fortran `module_rootdepth.f90` 翻译而来。

## 模块结构

Julia 版本将原始的单一 Fortran 模块分解为多个独立的、可测试的模块：

### 核心模块

1. **SoilParameters.jl** - 土壤参数模块
   - 定义13种土壤类型的参数
   - 计算土壤导水率
   - 初始化土壤参数

2. **Evapotranspiration.jl** - 蒸散发计算模块
   - Priestley-Taylor 方法
   - Penman-Monteith 方法
   - Shuttleworth-Wallace 双源法

3. **Interception.jl** - 截留模块
   - 植被截留计算
   - 截留蒸发计算

4. **WaterExtraction.jl** - 水分提取模块
   - 植物根系水分提取
   - 土壤水分胁迫计算

5. **SoilFluxes.jl** - 土壤水流模块
   - 土壤水分运动方程
   - 三对角矩阵求解器
   - 氧18同位素追踪

6. **WaterTableDynamics.jl** - 地下水位动态模块
   - 浅层地下水位更新
   - 侧向流影响

7. **SoilInitialization.jl** - 土壤初始化模块
   - CLM 方法土壤层设置
   - 固定层厚方法

8. **RootDepth.jl** - 主模块
   - 整合所有子模块
   - 提供主要计算接口

## 使用方法

```julia
# 获取土壤参数
soil_params = get_soil_params(5)  # 第5种土壤类型

# 计算蒸散发
pet = potevap_penman_monteith(1, 1, 298.15, 200.0, 300.0, 101325.0, 0.01, 2.0, 3.0, 5.0, 15.0)

# 计算截留
ppdrip, et_i, new_store = interception(0.01, 5.0, 3.0, 0.2, 0.5)
```
<!-- 
### 主要特性：

- **完整的物理过程**：保留了原始模型的所有物理计算

- **多种蒸散发方法**：支持三种不同的蒸散发计算方法

- **同位素追踪**：包含氧18同位素的追踪计算

- **灵活的土壤配置**：支持多种土壤层配置方案

- **数值稳定性**：改进的数值算法确保计算稳定 -->
