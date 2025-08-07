# ASAP 模型 Julia 版本

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

### 测试模块

每个核心模块都有对应的测试文件：

- `test_soil_parameters.jl`
- `test_evapotranspiration.jl`
- `test_interception.jl`
- `test_soil_initialization.jl`

## 使用方法

### 快速开始

```julia
# 运行示例
include("src/example_usage.jl")

# 运行所有测试
include("src/run_tests.jl")
```

### 基本使用

```julia
using .SoilParameters
using .Evapotranspiration
using .Interception

# 获取土壤参数
soil_params = get_soil_params(5)  # 第5种土壤类型

# 计算蒸散发
pet = potevap_penman_monteith(1, 1, 298.15, 200.0, 300.0, 101325.0, 0.01, 2.0, 3.0, 5.0, 15.0)

# 计算截留
ppdrip, et_i, new_store = interception(0.01, 5.0, 3.0, 0.2, 0.5)
```

## 优势

### 相比原始 Fortran 版本的改进：

1. **模块化设计**：每个功能模块独立，便于理解和维护
2. **类型安全**：Julia 的类型系统提供更好的错误检查
3. **可测试性**：每个模块都有完整的单元测试
4. **文档化**：详细的中文注释和文档字符串
5. **易于扩展**：面向对象的设计便于添加新功能
6. **性能优化**：Julia 的即时编译器提供接近 C/Fortran 的性能

### 主要特性：

- **完整的物理过程**：保留了原始模型的所有物理计算
- **多种蒸散发方法**：支持三种不同的蒸散发计算方法
- **同位素追踪**：包含氧18同位素的追踪计算
- **灵活的土壤配置**：支持多种土壤层配置方案
- **数值稳定性**：改进的数值算法确保计算稳定

## 文件结构

```
src/
├── SoilParameters.jl          # 土壤参数模块
├── Evapotranspiration.jl      # 蒸散发计算模块
├── Interception.jl            # 截留模块
├── WaterExtraction.jl         # 水分提取模块
├── SoilFluxes.jl             # 土壤水流模块
├── WaterTableDynamics.jl      # 地下水位动态模块
├── SoilInitialization.jl     # 土壤初始化模块
├── RootDepth.jl              # 主模块
├── test_*.jl                 # 测试文件
├── run_tests.jl              # 测试运行脚本
└── example_usage.jl          # 使用示例
```

## 开发指南

### 添加新功能

1. 在相应模块中添加新函数
2. 添加详细的文档字符串
3. 编写相应的测试
4. 更新 export 列表

### 运行测试

```julia
# 运行所有测试
include("src/run_tests.jl")

# 运行特定模块测试
include("src/test_soil_parameters.jl")
```

### 性能优化

- 使用类型注解提高性能
- 避免全局变量
- 使用 `@inbounds` 和 `@simd` 宏优化循环
- 预分配数组避免内存分配

## 依赖关系

- Julia 1.6+
- Test.jl (标准库)

## 许可证

本项目遵循与原始 ASAP 模型相同的许可证。

## 贡献

欢迎提交问题报告和功能请求。请确保新代码包含适当的测试和文档。

## 联系信息

如有疑问或建议，请联系项目维护者。
