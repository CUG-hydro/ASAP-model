# 测试获取土壤参数
@testset "获取土壤参数" begin
    soil_params = get_soil_params(1)
    @test soil_params.slmsts ≈ 0.395
    @test soil_params.soilcp ≈ 0.050
    @test soil_params.slbs ≈ 4.05
    @test soil_params.slcons ≈ 0.000176
    @test soil_params.slpots ≈ -0.121
    @test soil_params.klatfactor ≈ 2.0

    # 测试边界条件
    @test_throws BoundsError get_soil_params(0)
    @test_throws BoundsError get_soil_params(14)
end


# 测试初始化土壤参数
@testset "初始化土壤参数" begin
    nzg = 10
    fieldcp, slwilt = init_soil_param(nzg)
    @test size(fieldcp) == (nzg, NSTYP)
    @test length(slwilt) == NSTYP
    @test all(slwilt .> 0)
end

# 测试导水率计算
@testset "导水率计算" begin
    smoi = 0.3
    nsoil = 1
    k = khyd(smoi, nsoil)
    @test k > 0.0
    @test k ≈ SLCONS[nsoil] * (smoi / SLMSTS[nsoil])^(2.0 * SLBS[nsoil] + 3.0)

    # 测试边界条件
    @test_throws BoundsError khyd(0.3, 0)
    @test_throws BoundsError khyd(0.3, 14)
end
