"""
ASAP 模型完整使用示例
展示如何设置和运行完整的根系深度计算
"""

println("ASAP 模型 Julia 完整版本使用示例")
println("=" * 50)

using ASAP

# 导入所需的子模块用于初始化

# 1. 设置模型网格
println("\n1. 设置模型参数")
println("-" * 20)

# 网格设置
is, ie = 1, 3  # 较小的测试网格
js, je = 1, 3
nzg = 15       # 15层土壤
freedrain = 0  # 非自由排水
maxinactivedays = 240  # 最大非活跃天数

# 时间设置
Δt = 3600.0    # 1小时时间步长
steps = 1.0    # 每次调用的步数
hour = 12      # 当前小时

println("网格大小: $(ie)×$(je), 土壤层数: $nzg")
println("时间步长: $(Δt/3600) 小时")

# 2. 初始化土壤层
println("\n2. 初始化土壤层")
println("-" * 20)

slz, dz = initializesoildepth(nzg)
println("土壤层深度范围: $(slz[1]) 到 $(slz[end]) m")
println("表层厚度: $(dz[end]) m, 底层厚度: $(dz[1]) m")

# 3. 初始化土壤参数
fieldcp, slwilt = init_soil_param(nzg)
println("田间持水量矩阵大小: $(size(fieldcp))")

# 4. 创建输入数据
println("\n3. 创建输入数据")
println("-" * 20)

# 陆地掩码 (1=陆地, 0=海洋)
landmask = ones(Int, ie, je)

# 植被参数
veg = ones(Float64, ie, je) * 5.0      # 落叶阔叶林
hveg = ones(Float64, ie, je) * 20.0    # 植被高度 20m
lai = ones(Float64, ie, je) * 4.0      # 叶面积指数 4

# 土壤类型 (第一层和第二层)
soiltxt = ones(Int, 2, ie, je) * 5     # 第5种土壤类型

# 气象强迫数据
wind = ones(Float64, ie, je) * 3.0     # 风速 3 m/s
temp = ones(Float64, ie, je) * 298.15  # 温度 25°C
qair = ones(Float64, ie, je) * 0.012   # 比湿
press = ones(Float64, ie, je) * 101325.0  # 大气压力
netrad = ones(Float64, ie, je) * 250.0    # 净辐射 250 W/m²
rshort = ones(Float64, ie, je) * 400.0    # 短波辐射 400 W/m²
precip = ones(Float64, ie, je) * 2.0      # 降水 2 mm

# 5. 初始化状态变量
println("\n4. 初始化状态变量")
println("-" * 20)

# 土壤含水量相关
smoi = ones(Float64, nzg, ie, je) * 0.25      # 初始含水量
smoieq = ones(Float64, nzg, ie, je) * 0.20    # 平衡含水量
smoiwtd = ones(Float64, ie, je) * 0.30        # 地下水含水量
wtd = ones(Float64, ie, je) * -1.5            # 地下水位深度 -1.5m

# 蒸散发相关
et_s = zeros(Float64, ie, je)          # 土壤蒸发
et_i = zeros(Float64, ie, je)          # 截留蒸发
et_c = zeros(Float64, ie, je)          # 冠层蒸腾
et_s_daily = zeros(Float64, ie, je)    # 日土壤蒸发
et_c_daily = zeros(Float64, ie, je)    # 日冠层蒸腾

# 水分平衡相关
qsrun = zeros(Float64, ie, je)         # 地表径流
waterdeficit = zeros(Float64, ie, je)  # 水分不足量
watext = zeros(Float64, nzg, ie, je)   # 各层水分提取
watextdeep = zeros(Float64, ie, je)    # 深层水分提取
rech = zeros(Float64, ie, je)          # 补给量
deeprech = zeros(Float64, ie, je)      # 深层补给
transptop = zeros(Float64, ie, je)     # 表层蒸腾

# 截留相关
intercepstore = zeros(Float64, ie, je) # 截留存储
ppacum = zeros(Float64, ie, je)        # 累计降水

# 侧向流相关
qlat = zeros(Float64, ie, je)          # 侧向流
qlatsum = zeros(Float64, ie, je)       # 累计侧向流
qsprings = zeros(Float64, ie, je)      # 泉水流

# 渗透相关
pppendepth = zeros(Float64, ie, je)           # 渗透深度
pppendepthold = zeros(Int8, ie, je)           # 渗透深度记录
infilk = zeros(Int8, ie, je)                  # 入渗层指标
infilflux = zeros(Float64, nzg, ie, je)       # 入渗通量
infilfluxday = zeros(Float64, nzg, ie, je)    # 日入渗通量
infilcounter = zeros(Int16, nzg, ie, je)      # 入渗计数器
upflux = zeros(Float64, nzg, ie, je)          # 上升通量

# 根系活性
inactivedays = zeros(Int, nzg+2, ie, je)      # 非活跃天数

# 其他变量
fdepth = ones(Float64, ie, je) * 2.0          # 根系深度因子
floodheight = zeros(Float64, ie, je)          # 洪水高度
qrf = zeros(Float64, ie, je)                  # 河流反馈流
delsfcwat = zeros(Float64, ie, je)            # 地表水变化
icefactor = zeros(Int8, ie, je, 15)           # 冰冻因子 (简化为15层)
wtdflux = zeros(Float64, ie, je)              # 地下水通量

# 同位素追踪
o18 = ones(Float64, nzg, ie, je) * 0.05       # 氧18含量
o18ratiopp = ones(Float64, ie, je) * 0.001    # 降水氧18比值
tempsfc = copy(temp)                          # 地表温度
qlato18 = zeros(Float64, ie, je)              # 侧向流氧18
transpo18 = zeros(Float64, ie, je)            # 蒸腾氧18

println("状态变量初始化完成")
println("土壤含水量范围: $(minimum(smoi)) - $(maximum(smoi))")
println("地下水位深度: $(minimum(wtd)) - $(maximum(wtd)) m")

# 6. 运行模型
println("\n5. 运行 ASAP 模型")
println("-" * 20)

try
    # 调用主函数
    rootdepth_main(freedrain, is, ie, js, je, nzg, slz, dz, Δt,
                  landmask, veg, hveg, soiltxt, wind, temp, qair, press, netrad,
                  rshort, lai, precip, qsrun, smoi, smoieq, smoiwtd, wtd,
                  waterdeficit, watext, watextdeep, rech, deeprech, et_s, et_i,
                  et_c, intercepstore, ppacum, pppendepth, pppendepthold, qlat,
                  qlatsum, qsprings, inactivedays, maxinactivedays, fieldcp,
                  fdepth, steps, floodheight, qrf, delsfcwat, icefactor,
                  wtdflux, et_s_daily, et_c_daily, transptop, infilk, infilflux,
                  infilfluxday, infilcounter, hour, o18, o18ratiopp, tempsfc,
                  qlato18, transpo18, upflux)
    
    println("✓ 模型运行成功！")
    
    # 7. 输出结果
    println("\n6. 输出结果")
    println("-" * 20)
    
    println("蒸散发结果 (mm):")
    println("  土壤蒸发: $(round(mean(et_s), digits=3))")
    println("  截留蒸发: $(round(mean(et_i), digits=3))")
    println("  冠层蒸腾: $(round(mean(et_c), digits=3))")
    println("  总蒸散发: $(round(mean(et_s .+ et_i .+ et_c), digits=3))")
    
    println("\n水分平衡 (mm):")
    println("  地表径流: $(round(mean(qsrun) * 1000, digits=3))")
    println("  水分不足: $(round(mean(waterdeficit), digits=3))")
    println("  补给量: $(round(mean(rech), digits=3))")
    
    println("\n土壤状态:")
    println("  表层含水量: $(round(mean(smoi[nzg,:,:]), digits=3))")
    println("  底层含水量: $(round(mean(smoi[1,:,:]), digits=3))")
    println("  地下水位: $(round(mean(wtd), digits=3)) m")
    
    println("\n截留状态 (mm):")
    println("  截留存储: $(round(mean(intercepstore), digits=3))")
    println("  累计降水: $(round(mean(ppacum), digits=3))")
    
    if any(o18 .> 0)
        println("\n同位素追踪:")
        println("  平均氧18含量: $(round(mean(o18), digits=6))")
        println("  蒸腾氧18: $(round(mean(transpo18), digits=3))")
    end
    
catch e
    println("✗ 模型运行出现错误:")
    println("  错误类型: $(typeof(e))")
    println("  错误信息: $e")
    
    if isa(e, MethodError)
        println("  这可能是由于子模块函数签名不匹配造成的")
        println("  请检查各子模块的函数参数是否与主函数调用一致")
    end
end

println("\n" * "=" * 50)
println("ASAP 模型示例运行完成")
println("此示例展示了完整的模型设置和调用过程")
println("在实际应用中，您需要:")
println("1. 从文件读取真实的气象和地理数据")
println("2. 设置适当的初始条件")
println("3. 在时间循环中调用模型")
println("4. 保存和分析输出结果")
