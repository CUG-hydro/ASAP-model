# 测试CLM方法
@testset "initializesoildepth_clm" begin
    nzg = 10
    slz, dz = initializesoildepth_clm(nzg)

    @test length(slz) == nzg + 1
    @test length(dz) == nzg
end

# 测试固定层厚方法
@testset "initializesoildepth" begin
    nzg = 40
    slz, dz = initializesoildepth(nzg)
    @test slz[1] ≈ -1000
end
