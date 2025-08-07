using Test

# 包含所有wtable测试
println("开始运行水位模块完整测试套件...")

@testset "Water Table Module Complete Tests" begin
    # 基础功能测试
    include("test_watertable.jl")
    
    # 河流路由测试
    include("test_river_routing.jl")
    
    # 同位素追踪测试
    include("test_isotope_tracing.jl")
end

println("✅ 水位模块所有测试通过!")
println()
println("测试总结：")
println("- 基础水位计算: ✓")
println("- 侧向流计算: ✓") 
println("- 河流-地下水交换: ✓")
println("- 运动波洪水路由: ✓")
println("- 扩散波洪水路由: ✓")
println("- 洪泛区漫流: ✓")
println("- 同位素传输: ✓")
println("- 质量守恒验证: ✓")
println()
println("模块状态: 已完成翻译并通过测试 🎉")
