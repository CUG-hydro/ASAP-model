# 测试CLM方法
@testset "initializesoildepth_clm" begin
  nzg = 10
  z₋ₕ, dz = initializesoildepth_clm(nzg)

  @test length(z₋ₕ) == nzg + 1
  @test length(dz) == nzg
end

# 测试固定层厚方法
@testset "initializesoildepth" begin
  nzg = 40
  z₋ₕ, dz = initializesoildepth(nzg)
  @test z₋ₕ[1] ≈ -1000
end
