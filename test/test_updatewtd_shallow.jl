using ASAP, Test

@testset "find_jwt基本功能" begin
    # 使用 nzg=10 的标准层深度
    z₋ₕ, _ = initializesoildepth(10)
    # z₋ₕ 从最深到表面 (负值到0)

    # 水位深于所有层（默认返回1）
    jwt_deep = find_jwt(z₋ₕ[1] - 1.0, z₋ₕ)
    @test jwt_deep == 1

    # 水位在顶层（接近地表）
    jwt_top = find_jwt(z₋ₕ[10] + 0.01, z₋ₕ)
    @test jwt_top >= 1

    # 水位正好在第一层中心
    z_center1 = 0.5 * (z₋ₕ[1] + z₋ₕ[2])  # 第1层中心
    jwt1 = find_jwt(z_center1, z₋ₕ)
    @test jwt1 == 2  # 水位在第1层，上边界是 z₋ₕ[2]

    # 水位在某个中间层
    z_center5 = 0.5 * (z₋ₕ[5] + z₋ₕ[6])  # 第5层中心
    jwt5 = find_jwt(z_center5, z₋ₕ)
    # 应该返回6（z₋ₕ[6] 是第5层的上边界）
    @test jwt5 == 6
end

@testset "find_jwt边界条件" begin
    z₋ₕ = [-3.0, -2.0, -1.0, -0.5, -0.1, 0.0]  # nzg=5 手工构造

    # wtd=-1.5: z₋ₕ[3]=-1.0 > -1.5? YES → jwt=3
    @test find_jwt(-1.5, z₋ₕ) == 3

    # wtd=-0.8: z₋ₕ[4]=-0.5 > -0.8? YES → jwt=4  
    @test find_jwt(-0.8, z₋ₕ) == 4

    # wtd=-3.5 (低于所有层边界): 无 z₋ₕ[k] > -3.5 → jwt=1 (default)
    @test find_jwt(-3.5, z₋ₕ) == 1

    # wtd=-0.05 (在最浅层内部): z₋ₕ[5]=-0.1 > -0.05? No → jwt=1 (default)
    # 因为 nzg=5 的 z₋ₕ 中没有比 -0.05 更大的元素（除了z₋ₕ[6]=0，但不检查）
    @test find_jwt(-0.05, z₋ₕ) == 1
end

@testset "updatewtd_shallow水位在有效层内" begin
    nzg = 20
    z₋ₕ, dz = initializesoildepth(nzg)
    soiltxt = 1
    soil = get_soil_params(soiltxt)
    fdepth = 10.0

    # 水位在某一中间层内
    # 选择一个合适的 wtd
    k_target = 10  # 第10层
    wtd_init = 0.5 * (z₋ₕ[k_target] + z₋ₕ[k_target+1])  # 层中心

    # 构造平衡含水量 (略低于饱和含水量)
    θ_sat_approx = soil.θ_sat * 0.9
    θ_eq = fill(θ_sat_approx * 0.6, nzg)  # 平衡含水量 < 饱和量

    # 土壤接近饱和
    θ = fill(θ_sat_approx * 0.95, nzg)
    θ_wtd = soil.θ_sat

    wtd_new, rech = updatewtd_shallow(nzg, 0, z₋ₕ, dz, soiltxt, θ_eq, θ_wtd, θ, wtd_init, fdepth)

    # 基本检查：函数应当返回有限值
    @test isfinite(wtd_new)
    @test isfinite(rech)
end

@testset "updatewtd_shallow自由排水模式" begin
    nzg = 20
    z₋ₕ, dz = initializesoildepth(nzg)
    soiltxt = 1
    soil = get_soil_params(soiltxt)
    fdepth = 10.0

    θ_eq = fill(soil.θ_sat * 0.6, nzg)
    θ = fill(soil.θ_sat * 0.3, nzg)
    θ_wtd = soil.θ_sat

    # 水位在最深层中心
    wtd_init = 0.5 * (z₋ₕ[1] + z₋ₕ[2])

    # 自由排水 (freedrain=1)
    wtd_new, rech = updatewtd_shallow(nzg, 1, z₋ₕ, dz, soiltxt, θ_eq, θ_wtd, θ, wtd_init, fdepth)

    @test isfinite(wtd_new)
    @test isfinite(rech)
end

@testset "cal_factor指数衰减函数" begin
    fdepth = 2.0
    # z = -1.5 (即模型中的 z + 1.5 = 0): factor = exp(0) = 1.0
    @test ASAP.cal_factor(-1.5, fdepth) ≈ 1.0

    # z = -0.5 (z + 1.5 = 1.0): factor = exp(0.5) > 1, 但 clamp 到 1.0
    f1 = ASAP.cal_factor(-0.5, fdepth)
    @test f1 ≈ 1.0  # clamped to 1.0

    # z = -3.5 (z + 1.5 = -1.0): factor = exp(-0.5) ≈ 0.607
    f2 = ASAP.cal_factor(-3.5, fdepth)
    @test f2 ≈ exp(-1.0) atol=1e-10

    # 非常深的层 → 接近下限 0.1
    f3 = ASAP.cal_factor(-100.0, fdepth)
    @test f3 ≈ 0.1

    # fdepth 增大 → 衰减更慢
    f_shallow = ASAP.cal_factor(-5.0, 1.0)
    f_deep = ASAP.cal_factor(-5.0, 10.0)
    @test f_deep > f_shallow
end
