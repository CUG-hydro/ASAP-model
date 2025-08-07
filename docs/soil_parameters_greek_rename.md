"""
土壤参数希腊字母重命名总结
==============================

## 结构体字段重命名

在 `SoilParameters.jl` 模块中，`SoilType` 结构体的字段已更新为希腊字母表示法：

### 字段名称对应表

| 原字段名     | 新字段名（希腊字母） | 物理含义     | 单位  |
| ------------ | -------------------- | ------------ | ----- |
| `slmsts`     | `θ_sat`              | 饱和含水量   | -     |
| `soilcp`     | `ρb`                 | 土壤容重     | g/cm³ |
| `slbs`       | `b`                  | 土壤b参数    | -     |
| `slcons`     | `Ksat`               | 饱和导水率   | m/s   |
| `slpots`     | `ψ`                  | 饱和基质势   | m     |
| `slwilt`     | `θ_wilt`             | 凋萎点含水量 | -     |
| `klatfactor` | `K_latfactor`        | 侧向流因子   | -     |

## 更新的模块

### 1. SoilParameters.jl ✅
- [x] 结构体字段重命名
- [x] `get_soil_params()` 函数更新
- [x] `khyd()` 函数签名更新

### 2. WaterExtraction.jl ✅
- [x] 所有 `soil_params.slmsts` → `soil_params.θ_sat`
- [x] 所有 `soil_params.slpots` → `soil_params.ψ`
- [x] 所有 `soil_params.slbs` → `soil_params.b`

### 3. SoilFluxes.jl ✅
- [x] 所有 `soil_params.soilcp` → `soil_params.ρb`
- [x] 所有 `soil_params.slcons` → `soil_params.Ksat`
- [x] 所有 `soil_params.slmsts` → `soil_params.θ_sat`
- [x] 所有 `soil_params.slpots` → `soil_params.ψ`
- [x] 所有 `soil_params.slbs` → `soil_params.b`

### 4. WaterTableDynamics.jl ✅
- [x] 所有 `soil_params.slmsts` → `soil_params.θ_sat`

## 函数签名更新

### khyd() 函数
```julia
# 修改前
function khyd(smoi::Float64, nsoil::Int)

# 修改后  
function khyd(smoi::Float64, θ_sat::Float64, Ksat::Float64, b::Float64)
```

## 希腊字母输入方法

在 Julia REPL 或支持的编辑器中：
- `\theta` + Tab → `θ`
- `\rho` + Tab → `ρ`
- `\psi` + Tab → `ψ`

## 优势

1. **符合科学标准**：使用国际通用的土壤物理学符号
2. **代码简洁性**：希腊字母比英文缩写更短更清晰
3. **避免歧义**：希腊字母具有明确的物理含义
4. **国际化友好**：不依赖特定语言的词汇

## 测试验证

更新后需要运行以下测试确保兼容性：
- [ ] `test_soil_parameters.jl`
- [ ] `test_water_extraction.jl`
- [ ] `test_soil_fluxes.jl`
- [ ] `test_water_table_dynamics.jl`
- [ ] 集成测试

## 注意事项

1. **向后兼容性**：这是一个破坏性更改，可能需要更新调用代码
2. **文档更新**：确保所有文档和注释使用新的字段名
3. **团队培训**：确保团队成员了解新的命名约定

这些更改使 ASAP 模型的土壤参数表示更加符合土壤物理学的标准符号体系。
"""
