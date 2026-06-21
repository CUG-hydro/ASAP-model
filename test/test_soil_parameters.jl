using ASAP, Test

# 测试获取土壤参数
@testset "获取土壤参数" begin
    soil = get_soil_params(1)
    @test soil.θ_sat ≈ 0.395
    @test soil.θ_cp ≈ 0.050
    @test soil.b ≈ 4.05
    @test soil.Ksat ≈ 0.000176
    @test soil.ψsat ≈ -0.121
    @test soil.K_latfactor ≈ 2.0

    # 测试所有13种土壤类型
    for i in 1:13
        soil_i = get_soil_params(i)
        @test soil_i.θ_sat > 0.0
        @test soil_i.θ_cp > 0.0
        @test soil_i.θ_sat > soil_i.θ_cp
        @test soil_i.Ksat > 0.0
        @test soil_i.b > 0.0
        @test soil_i.K_latfactor > 0.0
        @test soil_i.ψsat < 0.0  # 基质势为负值
    end

    # 测试边界条件
    @test_throws ErrorException get_soil_params(0)
    @test_throws ErrorException get_soil_params(14)
end


# 测试初始化土壤参数
@testset "初始化土壤参数" begin
    nzg = 10
    fieldcp, θ_wilt = init_soil_param(nzg)
    @test size(fieldcp) == (nzg, ASAP.NSTYP)
    @test length(θ_wilt) == ASAP.NSTYP
    @test all(θ_wilt .> 0)
    # 凋萎点应小于饱和含水量
    for i in 1:ASAP.NSTYP
        soil = get_soil_params(i)
        @test θ_wilt[i] < soil.θ_sat
    end

    # fieldcp 应使用 Campbell 公式 θ_fc = θ_sat · (ψ_sat / -3.366)^(1/b)
    @test all(fieldcp .> 0.0)
    @test all(fieldcp .< [soil.θ_sat for soil in get_soil_params.(1:ASAP.NSTYP)]')
    for nsoil in 1:ASAP.NSTYP
        soil = get_soil_params(nsoil)
        expected_fc = soil.θ_sat * (soil.ψsat / -3.366)^(1.0 / soil.b)
        @test all(fieldcp[:, nsoil] .≈ expected_fc)
        # 田间持水量应大于凋萎点
        @test expected_fc > θ_wilt[nsoil]
    end
end

# 测试导水率计算
@testset "导水率计算" begin
    soil = get_soil_params(1)

    # 完全饱和时，K = Ksat
    k_sat = cal_K(soil.θ_sat, soil.θ_sat, soil.Ksat, soil.b)
    @test k_sat ≈ soil.Ksat

    # 含水量越高，导水率越大
    k1 = cal_K(0.2, soil.θ_sat, soil.Ksat, soil.b)
    k2 = cal_K(0.3, soil.θ_sat, soil.Ksat, soil.b)
    @test k2 > k1

    # 导水率应为正值
    @test k1 > 0.0
    @test k2 > 0.0

    # 验证公式: K = Ksat * (θ/θ_sat)^(2b+3)
    θ = 0.25
    k_expected = soil.Ksat * (θ / soil.θ_sat)^(2.0 * soil.b + 3.0)
    @test cal_K(θ, soil.θ_sat, soil.Ksat, soil.b) ≈ k_expected
end
