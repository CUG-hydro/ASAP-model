"""
ASAP 模型完整使用示例
展示如何设置和运行完整的根系深度计算

本示例的调用签名与 `test/test_rootdepth.jl::make_rootdepth_inputs` 完全对齐：
- nzg = 40（标准层数，与 icefactor[:, :, 26:40] 索引匹配）
- fieldcp 不作为形参传入 rootdepth_main（由 SoilParameters 模块内部维护）
"""

println("ASAP 模型 Julia 完整版本使用示例")
println(repeat("=", 50))

using ASAP, Statistics

# 1. 设置模型网格
println("\n1. 设置模型参数")
println(repeat("-", 20))

# 网格设置
is, ie = 1, 3  # 较小的测试网格
js, je = 1, 3
nzg    = 40    # 标准层数（与 icefactor[:, :, 26:40] 索引匹配）
freedrain = 1  # 自由排水（测试覆盖全路径）
maxinactivedays = 240

# 时间设置
Δt    = 3600.0   # 1 小时时间步长
steps = 1.0
hour  = 0

println("网格大小: $(ie)×$(je), 土壤层数: $nzg")
println("时间步长: $(Δt/3600) 小时")

# 2. 初始化土壤层
println("\n2. 初始化土壤层")
println(repeat("-", 20))

slz, dz = initializesoildepth(nzg)
println("土壤层深度范围: $(slz[1]) 到 $(slz[end]) m")
println("表层厚度: $(dz[end]) m, 底层厚度: $(dz[1]) m")

# 3. 初始化土壤参数（验证但不传入 rootdepth_main）
fieldcp, slwilt = init_soil_param(nzg)
println("田间持水量矩阵大小: $(size(fieldcp))")

# 4. 创建输入数据（与 make_rootdepth_inputs 一致）
println("\n3. 创建输入数据")
println(repeat("-", 20))

# 陆地掩码 (1=陆地, 0=海洋)
landmask = ones(Int, ie, je)

# 植被参数（中间格点 = 7 = 短草，覆盖完整路径）
veg  = fill(7.0,  ie, je)
hveg = fill(1.0,  ie, je)
lai  = fill(2.0,  ie, je)

# 土壤类型
soiltxt = ones(Int, 2, ie, je)

# 气象强迫
wind   = fill(3.0,      ie, je)
temp   = fill(293.0,    ie, je)
qair   = fill(0.01,     ie, je)
press  = fill(101325.0, ie, je)
netrad = fill(100.0,    ie, je)
rshort = fill(150.0,    ie, je)
precip = fill(5.0,      ie, je)

# 5. 初始化状态变量
println("\n4. 初始化状态变量")
println(repeat("-", 20))

qsrun        = zeros(Float64, ie, je)
smoi         = fill(0.25, nzg, ie, je)
smoieq       = fill(0.20, nzg, ie, je)
smoiwtd      = fill(0.35, ie, je)
wtd          = fill(-1.0, ie, je)
waterdeficit = zeros(Float64, ie, je)
watext       = zeros(Float64, nzg, ie, je)
watextdeep   = zeros(Float64, ie, je)
rech         = zeros(Float64, ie, je)
deeprech     = zeros(Float64, ie, je)
et_s         = zeros(Float64, ie, je)
et_i         = zeros(Float64, ie, je)
et_c         = zeros(Float64, ie, je)
et_s_daily   = zeros(Float64, ie, je)
et_c_daily   = zeros(Float64, ie, je)
transptop    = zeros(Float64, ie, je)

intercepstore   = zeros(Float64, ie, je)
ppacum          = zeros(Float64, ie, je)
pppendepth      = zeros(Float64, ie, je)
pppendepthold   = zeros(Int8,    ie, je)

qlat      = zeros(Float64, ie, je)
qlatsum   = zeros(Float64, ie, je)
qsprings  = zeros(Float64, ie, je)

inactivedays    = zeros(Int, nzg+2, ie, je)

fdepth      = fill(2.0, ie, je)
floodheight = zeros(Float64, ie, je)
qrf         = zeros(Float64, ie, je)
delsfcwat   = zeros(Float64, ie, je)
icefactor   = zeros(Int8,  ie, je, nzg)  # 注：函数内部使用 icefactor[i, j, 26:40]
wtdflux     = zeros(Float64, ie, je)

infilk       = zeros(Int8,    ie, je)
infilflux    = zeros(Float64, nzg, ie, je)
infilfluxday = zeros(Float64, nzg, ie, je)
infilcounter = zeros(Int16,   nzg, ie, je)
upflux       = zeros(Float64, nzg, ie, je)

# 同位素追踪（占位，函数签名要求但 o18 段当前未启用）
o18        = zeros(Float64, nzg, ie, je)
o18ratiopp = zeros(Float64, ie, je)
tempsfc    = fill(293.0, ie, je)
qlato18    = zeros(Float64, ie, je)
transpo18  = zeros(Float64, ie, je)

println("状态变量初始化完成")
println("土壤含水量范围: $(minimum(smoi)) - $(maximum(smoi))")
println("地下水位深度: $(minimum(wtd)) - $(maximum(wtd)) m")

# 6. 运行模型
println("\n5. 运行 ASAP 模型")
println(repeat("-", 20))

try
    rootdepth_main(freedrain, is, ie, js, je, nzg, slz, dz, Δt,
                   landmask, veg, hveg,
                   soiltxt, wind, temp,
                   qair, press, netrad,
                   rshort, lai, precip,
                   qsrun, smoi, smoieq,
                   smoiwtd, wtd, waterdeficit,
                   watext,
                   watextdeep, rech,
                   deeprech,
                   et_s, et_i, et_c,
                   intercepstore, ppacum,
                   pppendepth, pppendepthold, qlat,
                   qlatsum, qsprings, inactivedays,
                   maxinactivedays, fdepth,
                   steps, floodheight, qrf,
                   delsfcwat, icefactor, wtdflux,
                   et_s_daily, et_c_daily, transptop,
                   infilk, infilflux, infilfluxday,
                   infilcounter, hour,
                   o18, o18ratiopp, tempsfc,
                   qlato18, transpo18, upflux)

    println("✓ 模型运行成功！")

    # 7. 输出结果
    println("\n6. 输出结果")
    println(repeat("-", 20))

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

catch e
    println("✗ 模型运行出现错误:")
    println("  错误类型: $(typeof(e))")
    println("  错误信息: ", sprint(showerror, e))

    if isa(e, MethodError)
        println("  这可能是由于子模块函数签名不匹配造成的")
        println("  请检查各子模块的函数参数是否与主函数调用一致")
    end
end

println("\n" * repeat("=", 50))
println("ASAP 模型示例运行完成")
println("此示例展示了完整的模型设置和调用过程")
println("在实际应用中，您需要:")
println("1. 从文件读取真实的气象和地理数据")
println("2. 设置适当的初始条件")
println("3. 在时间循环中调用模型")
println("4. 保存和分析输出结果")