"""
ASAP 模型泛型类型参数示例
展示如何使用 V (Vector), M (Matrix), A3 (Array3D) 类型参数
"""

println("🧮 ASAP 模型泛型类型参数示例")
println("=" * 50)

# 导入主模块
push!(LOAD_PATH, ".")
include("ASAP.jl")
using .ASAP

# 导入初始化模块
include("SoilInitialization.jl")
include("SoilParameters.jl")
using .SoilInitialization
using .SoilParameters

println("\n📊 演示不同数值精度类型的使用")
println("-" * 40)

# 模型参数设置
is, ie = 1, 2
js, je = 1, 2
nzg = 5  # 使用较少的层数以简化示例
freedrain = 0
maxinactivedays = 240
Δt = 3600.0
steps = 1.0
hour = 12

println("网格: $(ie)×$(je), 土壤层数: $nzg")

# 初始化土壤层
slz_f64, dz_f64 = initializesoildepth(nzg)
fieldcp_f64, _ = init_soil_param(nzg)

println("\n1️⃣  使用 Float64 类型 (高精度)")
println("-" * 25)

# Float64 类型数据
slz_64 = Vector{Float64}(slz_f64)
dz_64 = Vector{Float64}(dz_f64)
fieldcp_64 = Matrix{Float64}(fieldcp_f64)

# 创建 Float64 类型的输入数据
function create_data_f64()
    # 标量输入
    landmask = ones(Int, ie, je)
    soiltxt = ones(Int, 2, ie, je) * 5
    icefactor = zeros(Int8, ie, je, 15)  # 简化版本
    pppendepthold = zeros(Int8, ie, je)
    infilk = zeros(Int8, ie, je)
    inactivedays = zeros(Int, nzg+2, ie, je)
    infilcounter = zeros(Int16, nzg, ie, je)
    
    # Float64 矩阵和数组
    veg = ones(Float64, ie, je) * 5.0
    hveg = ones(Float64, ie, je) * 20.0
    wind = ones(Float64, ie, je) * 3.0
    temp = ones(Float64, ie, je) * 298.15
    qair = ones(Float64, ie, je) * 0.012
    press = ones(Float64, ie, je) * 101325.0
    netrad = ones(Float64, ie, je) * 250.0
    rshort = ones(Float64, ie, je) * 400.0
    lai = ones(Float64, ie, je) * 4.0
    precip = ones(Float64, ie, je) * 2.0
    
    # 初始化状态变量
    qsrun = zeros(Float64, ie, je)
    smoi = ones(Float64, nzg, ie, je) * 0.25
    smoieq = ones(Float64, nzg, ie, je) * 0.20
    smoiwtd = ones(Float64, ie, je) * 0.30
    wtd = ones(Float64, ie, je) * -1.5
    waterdeficit = zeros(Float64, ie, je)
    watext = zeros(Float64, nzg, ie, je)
    watextdeep = zeros(Float64, ie, je)
    rech = zeros(Float64, ie, je)
    deeprech = zeros(Float64, ie, je)
    et_s = zeros(Float64, ie, je)
    et_i = zeros(Float64, ie, je)
    et_c = zeros(Float64, ie, je)
    intercepstore = zeros(Float64, ie, je)
    ppacum = zeros(Float64, ie, je)
    pppendepth = zeros(Float64, ie, je)
    qlat = zeros(Float64, ie, je)
    qlatsum = zeros(Float64, ie, je)
    qsprings = zeros(Float64, ie, je)
    fdepth = ones(Float64, ie, je) * 2.0
    floodheight = zeros(Float64, ie, je)
    qrf = zeros(Float64, ie, je)
    delsfcwat = zeros(Float64, ie, je)
    wtdflux = zeros(Float64, ie, je)
    et_s_daily = zeros(Float64, ie, je)
    et_c_daily = zeros(Float64, ie, je)
    transptop = zeros(Float64, ie, je)
    infilflux = zeros(Float64, nzg, ie, je)
    infilfluxday = zeros(Float64, nzg, ie, je)
    o18 = ones(Float64, nzg, ie, je) * 0.05
    o18ratiopp = ones(Float64, ie, je) * 0.001
    tempsfc = copy(temp)
    qlato18 = zeros(Float64, ie, je)
    transpo18 = zeros(Float64, ie, je)
    upflux = zeros(Float64, nzg, ie, je)
    
    return (landmask, veg, hveg, soiltxt, wind, temp, qair, press, netrad, rshort, lai, precip,
            qsrun, smoi, smoieq, smoiwtd, wtd, waterdeficit, watext, watextdeep, rech, deeprech,
            et_s, et_i, et_c, intercepstore, ppacum, pppendepth, pppendepthold, qlat, qlatsum,
            qsprings, inactivedays, fdepth, floodheight, qrf, delsfcwat, icefactor, wtdflux,
            et_s_daily, et_c_daily, transptop, infilk, infilflux, infilfluxday, infilcounter,
            o18, o18ratiopp, tempsfc, qlato18, transpo18, upflux)
end

# 测试 Float64 版本
println("创建 Float64 类型数据...")
data_f64 = create_data_f64()

println("类型验证:")
println("  slz_64: $(typeof(slz_64))")
println("  fieldcp_64: $(typeof(fieldcp_64))")
println("  smoi: $(typeof(data_f64[14]))")

println("调用 rootdepth_main (Float64)...")
try
    @time rootdepth_main(freedrain, is, ie, js, je, nzg, slz_64, dz_64, Δt,
                        data_f64...)
    println("✅ Float64 版本运行成功!")
catch e
    println("❌ Float64 版本运行失败: $e")
end

println("\n2️⃣  使用 Float32 类型 (高性能)")
println("-" * 25)

# Float32 类型数据 
slz_32 = Vector{Float32}(slz_f64)
dz_32 = Vector{Float32}(dz_f64)
fieldcp_32 = Matrix{Float32}(fieldcp_f64)

# 创建 Float32 类型的输入数据
function create_data_f32()
    # 重用之前的函数，但转换为 Float32
    data = create_data_f64()
    
    # 转换相关的数组为 Float32
    converted = []
    for item in data
        if isa(item, Matrix{Float64})
            push!(converted, Matrix{Float32}(item))
        elseif isa(item, Array{Float64,3})
            push!(converted, Array{Float32,3}(item))
        else
            push!(converted, item)  # 保持整数类型不变
        end
    end
    
    return tuple(converted...)
end

println("创建 Float32 类型数据...")
data_f32 = create_data_f32()

println("类型验证:")
println("  slz_32: $(typeof(slz_32))")
println("  fieldcp_32: $(typeof(fieldcp_32))")
println("  smoi: $(typeof(data_f32[14]))")

println("调用 rootdepth_main (Float32)...")
try
    @time rootdepth_main(freedrain, is, ie, js, je, nzg, slz_32, dz_32, Δt,
                        data_f32...)
    println("✅ Float32 版本运行成功!")
    
    # 比较内存使用
    mem_f64 = sizeof(data_f64[14]) + sizeof(slz_64) + sizeof(fieldcp_64)
    mem_f32 = sizeof(data_f32[14]) + sizeof(slz_32) + sizeof(fieldcp_32)
    
    println("\n📈 性能对比:")
    println("  Float64 内存使用: $(mem_f64) 字节")
    println("  Float32 内存使用: $(mem_f32) 字节")
    println("  内存节省: $(round((1 - mem_f32/mem_f64) * 100, digits=1))%")
    
catch e
    println("❌ Float32 版本运行失败: $e")
end

println("\n🎯 泛型类型参数的优势:")
println("• 类型安全: 编译时检查类型一致性")
println("• 性能优化: Float32 可减少内存使用和提高计算速度")  
println("• 代码复用: 同一函数支持多种数值精度")
println("• 灵活性: 可根据精度需求选择合适的类型")

println("\n" * "=" * 50)
println("🎉 泛型类型参数示例完成")
