using ASAP, Test

# ===========================================================================
# zbrent 基本求根
# ===========================================================================
@testset "zbrent basic root finding" begin
  # f(x) = x^2 - 2, 根为 √2 ≈ 1.4142135623730951
  root = zbrent(x -> x^2 - 2.0, 0.0, 2.0; tol=1e-12)
  @test root ≈ sqrt(2.0) atol=1e-10
  @test (root^2 - 2.0) ≈ 0.0 atol=1e-10

  # 在 [1, 2] 区间
  root2 = zbrent(x -> x^2 - 2.0, 1.0, 2.0; tol=1e-12)
  @test root2 ≈ sqrt(2.0) atol=1e-10

  # 不同函数: f(x) = sin(x) - 0.5, 根为 π/6
  root3 = zbrent(x -> sin(x) - 0.5, 0.0, 1.0; tol=1e-12)
  @test root3 ≈ π/6 atol=1e-10

  # 线性函数: f(x) = 3x - 6, 根为 2
  root4 = zbrent(x -> 3.0*x - 6.0, 0.0, 5.0; tol=1e-12)
  @test root4 ≈ 2.0 atol=1e-10
end

# ===========================================================================
# zbrent 异常处理
# ===========================================================================
@testset "zbrent error handling" begin
  # 同号端点 → ArgumentError
  @test_throws ArgumentError zbrent(x -> x^2 + 1.0, 0.0, 2.0)
  @test_throws ArgumentError zbrent(x -> x^2 + 1.0, -2.0, 0.0)
  @test_throws ArgumentError zbrent(x -> exp(x), 0.0, 1.0)
end

# ===========================================================================
# eqsoilmoisturetheor 烟测
# ===========================================================================
@testset "eqsoilmoisturetheor smoke test" begin
  nzg = 40
  nsoil = 5
  fdepth = 2.0
  wtd = -1.0

  z₋ₕ, dz = initializesoildepth(nzg)
  smoieq = eqsoilmoisturetheor(nzg, nsoil, z₋ₕ, dz, fdepth, wtd)

  # 形状
  @test length(smoieq) == nzg
  @test eltype(smoieq) == Float64

  # 物理范围:含深度衰减的 [smoicp_local, smoisat_local]
  # 深层由于 exp((vctr4+1.5)/fdepth) 衰减,smoicp 可能下探到
  # `soilcp(nsoil) * 0.1`,smoisat 同样下探。
  θ_cp_floor = soilcp(nsoil) * 0.05   # 留点裕度
  θ_sat_ceil = slmsts(nsoil)
  @test all(θ_cp_floor .≤ smoieq .≤ θ_sat_ceil + 1e-9)

  # 没有 NaN / Inf
  @test all(isfinite.(smoieq))

  # 不是常量(应当随深度变化)
  @test length(unique(round.(smoieq, digits=6))) > 1
end

# ===========================================================================
# eqsoilmoisturetheor 守恒:5 种土壤类型
# ===========================================================================
@testset "eqsoilmoisturetheor conservation across soil types" begin
  nzg = 40
  fdepth = 2.0
  wtd = -1.0
  z₋ₕ, dz = initializesoildepth(nzg)

  for nsoil in (1, 3, 5, 8, 11)
    smoieq = eqsoilmoisturetheor(nzg, nsoil, z₋ₕ, dz, fdepth, wtd)

    @test length(smoieq) == nzg
    @test all(isfinite.(smoieq))
    @test all(0.0 .≤ smoieq .≤ slmsts(nsoil) + 1e-9)

    # 平均含水量应在合理范围
    @test 0.0 < sum(smoieq)/length(smoieq) < slmsts(nsoil)
  end
end
