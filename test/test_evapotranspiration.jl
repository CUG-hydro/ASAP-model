"""
蒸散发模块测试
"""
module TestEvapotranspiration

using Test
include("../src/Evapotranspiration.jl")
using .Evapotranspiration

function test_evapotranspiration()
    @testset "蒸散发计算测试" begin
        # 测试 Priestley-Taylor 方法
        @testset "Priestley-Taylor 方法" begin
            i, j = 1, 1
            tempk = 298.15  # 25°C
            rad = 200.0     # W/m²
            presshp = 1013.25  # hPa
            
            pet = potevap_priestly_taylor(i, j, tempk, rad, presshp)
            @test pet > 0.0
            @test pet < 50.0  # 合理范围检查
        end
        
        # 测试 Penman-Monteith 方法
        @testset "Penman-Monteith 方法" begin
            i, j = 1, 1
            tempk = 298.15  # 25°C
            rad = 200.0     # W/m²
            rshort = 300.0  # W/m²
            press = 101325.0  # Pa
            qair = 0.01     # kg/kg
            wind = 2.0      # m/s
            lai = 3.0
            veg = 5.0       # 落叶阔叶林
            hveg = 15.0     # m
            
            pet = potevap_penman_monteith(i, j, tempk, rad, rshort, press, qair, wind, lai, veg, hveg)
            @test pet > 0.0
            @test pet < 100.0  # 合理范围检查
        end
        
        # 测试 Shuttleworth-Wallace 方法
        @testset "Shuttleworth-Wallace 方法" begin
            i, j = 1, 1
            deltat = 3600.0  # 1小时
            tempk = 298.15   # 25°C
            rad = 200.0      # W/m²
            rshort = 300.0   # W/m²
            press = 101325.0 # Pa
            qair = 0.01      # kg/kg
            wind = 2.0       # m/s
            lai = 3.0
            veg = 5.0        # 落叶阔叶林
            hhveg = 15.0     # m
            floodflag = 0
            
            result = potevap_shutteworth_wallace(i, j, deltat, tempk, rad, rshort, press, 
                                               qair, wind, lai, veg, hhveg, floodflag)
            
            delta, gamma, lambda, ra_a, ra_c, rs_c, R_a, R_s, pet_s, pet_c, pet_w, pet_i = result
            
            @test delta > 0.0
            @test gamma > 0.0
            @test lambda > 0.0
            @test ra_a > 0.0
            @test ra_c >= 0.0
            @test rs_c > 0.0
            @test pet_i >= 0.0
            @test pet_w == 0.0  # 植被情况下应为0
        end
        
        # 测试水体情况
        @testset "水体蒸发" begin
            i, j = 1, 1
            deltat = 3600.0
            tempk = 298.15
            rad = 200.0
            rshort = 300.0
            press = 101325.0
            qair = 0.01
            wind = 2.0
            lai = 0.0
            veg = 1.0        # 水体
            hhveg = 0.1
            floodflag = 0
            
            result = potevap_shutteworth_wallace(i, j, deltat, tempk, rad, rshort, press, 
                                               qair, wind, lai, veg, hhveg, floodflag)
            
            delta, gamma, lambda, ra_a, ra_c, rs_c, R_a, R_s, pet_s, pet_c, pet_w, pet_i = result
            
            @test pet_w > 0.0
            @test pet_s == 0.0
            @test pet_c == 0.0
            @test pet_i == 0.0
        end
    end
end

# 运行测试
if abspath(PROGRAM_FILE) == @__FILE__
    test_evapotranspiration()
    println("蒸散发模块测试通过！")
end

end # module TestEvapotranspiration
