"""
希腊字母变量名使用总结
====================

Julia 支持直接使用 Unicode 希腊字母作为变量名，这使得数学代码更加清晰和符合标准数学符号。

## 已修正的希腊字母变量名

### 1. 时间相关
- `deltat` → `Δt` (delta t, 时间步长)

### 2. 蒸散发相关希腊字母
- `alpha` → `α` (alpha, Priestley-Taylor 系数)
- `delta` → `δ` (delta, 饱和蒸汽压斜率)
- `gamma` → `γ` (gamma, 干湿表常数)
- `lambda` → `λ` (lambda, 潜热)

### 3. 在函数中的使用示例

#### Evapotranspiration.jl
```julia
# 修正前
function potevap_priestly_taylor(i::Int, j::Int, tempk::Float64, rad::Float64, presshp::Float64)
    alpha = 1.26
    delta = 0.2 * (0.00738 * tempc + 0.8072)^7 - 0.000116
    lambda = 2.501 - 0.002361 * tempc
    gamma = (CP_LOCAL * presskp) / (0.622 * lambda)
    pet = alpha * rad_mj * delta / (delta + gamma)
    pet = pet / lambda
end

# 修正后
function potevap_priestly_taylor(i::Int, j::Int, tempk::Float64, rad::Float64, presshp::Float64)
    α = 1.26
    δ = 0.2 * (0.00738 * tempc + 0.8072)^7 - 0.000116
    λ = 2.501 - 0.002361 * tempc
    γ = (CP_LOCAL * presskp) / (0.622 * λ)
    pet = α * rad_mj * δ / (δ + γ)
    pet = pet / λ
end
```

#### WaterExtraction.jl
```julia
# 函数参数中使用希腊字母
function extraction(i::Int, j::Int, nzg::Int, slz::Vector{Float64}, dz::Vector{Float64}, 
                   Δt::Float64, soiltxt::Int, wtd::Float64, smoi::Vector{Float64}, 
                   smoiwtd::Float64, δ::Float64, γ::Float64, λ::Float64, ...)
    # 函数体中直接使用希腊字母变量
    pet_c = max(Δt * pet_c / λ, 0.0)
    pet_s = max(Δt * pet_s / λ, 0.0)
end
```

#### ASAP.jl (主模块)
```julia
# 变量赋值使用希腊字母
δ, γ, λ, ra_a, ra_c, rs_c, R_a, R_s, petfactor_s, petfactor_c, petstep_w, petstep_i =
    potevap_shutteworth_wallace(Δt, temp[i, j], netrad[i, j], rshort[i, j], ...)
```

## 修正的优势

### 1. 数学清晰性
- 符合标准数学文献中的符号约定
- 代码与数学公式直接对应，减少理解负担

### 2. 代码简洁性
- 希腊字母通常比对应的英文单词更短
- 减少变量名冲突的可能性

### 3. 国际化友好
- 希腊字母是国际通用的数学符号
- 不依赖特定语言的词汇

### 4. Julia 语言特性
- Julia 原生支持 Unicode 字符
- 可以直接在编辑器中输入希腊字母（如 \alpha + Tab → α）

## 常用希腊字母输入方法

在支持 Julia 的编辑器中：
- `\alpha` + Tab → `α`
- `\beta` + Tab → `β`  
- `\gamma` + Tab → `γ`
- `\delta` + Tab → `δ`
- `\lambda` + Tab → `λ`
- `\Delta` + Tab → `Δ`

## 注意事项

1. **一致性**: 在整个项目中保持希腊字母使用的一致性
2. **文档**: 在文档中明确说明希腊字母变量的含义
3. **团队约定**: 确保团队成员都理解希腊字母变量的含义
4. **编辑器支持**: 确保使用的编辑器支持 Unicode 字符输入

## 测试验证

所有修改后的函数都已通过单元测试验证：
- ✅ Priestley-Taylor 蒸散发计算
- ✅ Penman-Monteith 蒸散发计算  
- ✅ Shuttleworth-Wallace 双源蒸散发计算
- ✅ 水分提取计算
- ✅ 主模块集成计算

这些修改使 ASAP 模型的 Julia 版本更加符合数学建模的标准规范，同时保持了代码的可读性和可维护性。
"""
