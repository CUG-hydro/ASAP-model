"""
土壤初始化模块测试
"""
module TestSoilInitialization

using Test
include("../src/SoilInitialization.jl")
using .SoilInitialization

function test_soil_initialization()
    @testset "土壤初始化测试" begin
        # 测试CLM方法
        @testset "CLM土壤层初始化" begin
            nzg = 10
            slz, dz = initializesoildepthclm(nzg)
            
            @test length(slz) == nzg + 1
            @test length(dz) == nzg
            @test slz[nzg + 1] ≈ 0.0  # 地表深度为0
            @test all(slz[1:nzg] .< 0.0)  # 所有土壤层深度为负值
            @test all(dz .> 0.0)  # 所有层厚度为正值
            
            # 检查深度单调性
            for k in 1:nzg
                @test slz[k] < slz[k + 1]
            end
            
            # 检查层厚度一致性
            for k in 1:nzg
                @test abs(slz[k + 1] - slz[k] - dz[k]) < 1e-10
            end
        end
        
        # 测试固定层厚方法
        @testset "固定层厚土壤初始化" begin
            nzg = 20
            slz, dz = initializesoildepth(nzg)
            
            @test length(slz) == nzg + 1
            @test length(dz) == nzg
            @test slz[nzg + 1] ≈ 0.0  # 地表深度为0
            @test all(slz[1:nzg] .< 0.0)  # 所有土壤层深度为负值
            @test all(dz .> 0.0)  # 所有层厚度为正值
            
            # 检查深度单调性
            for k in 1:nzg
                @test slz[k] < slz[k + 1]
            end
            
            # 检查层厚度一致性
            for k in 1:nzg
                @test abs(slz[k + 1] - slz[k] - dz[k]) < 1e-10
            end
            
            # 检查预定义的层厚度
            expected_dz = [0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2, 
                          0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4, 0.5, 0.5, 0.6]
            @test dz ≈ expected_dz[end-nzg+1:end]
        end
        
        # 测试边界条件
        @testset "边界条件测试" begin
            # 测试最小层数
            nzg = 1
            slz, dz = initializesoildepthclm(nzg)
            @test length(slz) == 2
            @test length(dz) == 1
            @test slz[2] ≈ 0.0
            @test slz[1] < 0.0
            @test dz[1] > 0.0
            
            # 测试最大允许层数
            nzg = 40
            slz, dz = initializesoildepth(nzg)
            @test length(slz) == 41
            @test length(dz) == 40
            
            # 测试超过最大层数
            @test_throws ErrorException initializesoildepth(41)
        end
        
        # 测试两种方法的一致性检查
        @testset "方法比较" begin
            nzg = 15
            slz1, dz1 = initializesoildepthclm(nzg)
            slz2, dz2 = initializesoildepth(nzg)
            
            # 两种方法应该给出不同的结果（因为算法不同）
            @test !(slz1 ≈ slz2)
            @test !(dz1 ≈ dz2)
            
            # 但都应该满足基本要求
            @test all(slz1[1:nzg] .< 0.0)
            @test all(slz2[1:nzg] .< 0.0)
            @test all(dz1 .> 0.0)
            @test all(dz2 .> 0.0)
        end
        
        # 测试深度分布特性
        @testset "深度分布特性" begin
            nzg = 20
            slz_clm, dz_clm = initializesoildepthclm(nzg)
            slz_fix, dz_fix = initializesoildepth(nzg)
            
            # CLM方法：表层应该更细
            @test dz_clm[end] < dz_clm[1]  # 表层（最后一个元素）比底层更薄
            
            # 固定层厚方法：表层应该是预定义的0.1m
            @test dz_fix[end] ≈ 0.1
            
            # 检查总深度
            total_depth_clm = abs(slz_clm[1])
            total_depth_fix = abs(slz_fix[1])
            @test total_depth_clm > 0.0
            @test total_depth_fix > 0.0
        end
    end
end

# 运行测试
if abspath(PROGRAM_FILE) == @__FILE__
    test_soil_initialization()
    println("土壤初始化模块测试通过！")
end

end # module TestSoilInitialization
