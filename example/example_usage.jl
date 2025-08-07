"""
ASAP 模型 Julia 版本使用示例
展示如何使用各个模块的基本功能
"""

# 添加当前目录到路径
push!(LOAD_PATH, ".")

# 导入各个模块

println("ASAP 模型 Julia 版本使用示例")
println("=" ^ 40)

# 1. 土壤参数示例
println("\n1. 土壤参数模块示例:")
println("-" ^ 20)
soil_type = 5  # 第5种土壤类型
soil_params = get_soil_params(soil_type)
println("土壤类型 $soil_type 的参数:")
println("  饱和含水量: $(soil_params.slmsts)")
println("  饱和导水率: $(soil_params.slcons)")
println("  土壤b参数: $(soil_params.slbs)")

# 计算导水率
smoi = 0.3  # 当前含水量
k = khyd(smoi, soil_type)
println("  当前导水率: $k")

# 2. 土壤层初始化示例
println("\n2. 土壤层初始化示例:")
println("-" ^ 20)
nzg = 10  # 10层土壤
slz, dz = initializesoildepth(nzg)
println("土壤层深度 (m):")
for i in 1:nzg
    println("  第$(i)层: $(slz[i]) 到 $(slz[i+1]), 厚度: $(dz[i])")
end

# 3. 蒸散发计算示例
println("\n3. 蒸散发计算示例:")
println("-" ^ 20)
# 设置气象条件
tempk = 298.15      # 温度 25°C
rad = 200.0         # 净辐射 200 W/m²
rshort = 300.0      # 短波辐射 300 W/m²
press = 101325.0    # 大气压力 101325 Pa
qair = 0.01         # 比湿 0.01
wind = 2.0          # 风速 2 m/s
lai = 3.0           # 叶面积指数 3
veg = 5.0           # 植被类型（落叶阔叶林）
hveg = 15.0         # 植被高度 15 m

# Priestley-Taylor 方法
pet_pt = potevap_priestly_taylor(tempk, rad * 0.1, press * 0.01)
println("Priestley-Taylor 蒸散发: $(round(pet_pt, digits=3)) mm")

# Penman-Monteith 方法
pet_pm = potevap_penman_monteith(tempk, rad, rshort, press, qair, wind, lai, veg, hveg)
println("Penman-Monteith 蒸散发: $(round(pet_pm, digits=3)) mm/3h")

# Shuttleworth-Wallace 方法
result = potevap_shutteworth_wallace(3600.0, tempk, rad, rshort, press, qair, wind, lai, veg, hveg, 0)
delta, gamma, lambda, ra_a, ra_c, rs_c, R_a, R_s, pet_s, pet_c, pet_w, pet_i = result
println("Shuttleworth-Wallace 截留蒸发: $(round(pet_i, digits=3)) mm/h")

# 4. 植被截留示例
println("\n4. 植被截留示例:")
println("-" ^ 20)
minpprate = 0.01    # 最小降水率阈值
precip = 5.0        # 降水量 5 mm
lai = 3.0           # 叶面积指数 3
intercepstore = 0.2 # 当前截留量 0.2 mm
pet_i = 0.5         # 截留蒸发潜力 0.5 mm

ppdrip, et_i, new_intercepstore = interception(minpprate, precip, lai, intercepstore, pet_i)
println("输入:")
println("  降水量: $precip mm")
println("  截留存储: $intercepstore mm")
println("  蒸发潜力: $pet_i mm")
println("输出:")
println("  穿透降水: $(round(ppdrip, digits=3)) mm")
println("  截留蒸发: $(round(et_i, digits=3)) mm")
println("  新截留量: $(round(new_intercepstore, digits=3)) mm")

# 5. 参数初始化示例
println("\n5. 土壤参数初始化示例:")
println("-" ^ 20)
fieldcp, slwilt = init_soil_param(nzg)
println("田间持水量矩阵大小: $(size(fieldcp))")
println("凋萎点含水量向量长度: $(length(slwilt))")
println("前3种土壤的凋萎点含水量:")
for i in 1:3
    println("  土壤类型 $i: $(round(slwilt[i], digits=4))")
end

println("\n" * "=" ^ 40)
println("示例运行完成！")
println("所有模块都可以正常使用。")
println("请参考各模块的文档了解更多功能。")
