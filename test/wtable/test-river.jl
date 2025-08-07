@testset "River-Groundwater Interaction" begin
    # 测试河流-地下水交换
    imax, jmax = 3, 3
    nzg = 2

    slz = [-0.1, -1.0, -2.0]
    soiltxt = ones(Int, 2, imax, jmax)
    landmask = ones(Int, imax, jmax)
    wtd = fill(-0.5, imax, jmax)  # 浅水位
    maxdepth = ones(imax, jmax) * 2.0
    riverdepth = ones(imax, jmax) * 1.0
    width = ones(imax, jmax) * 10.0
    length = ones(imax, jmax) * 100.0
    area = ones(imax, jmax) * 1000.0
    fdepth = ones(imax, jmax) * 3.0
    qrf = zeros(imax, jmax)

    Δt = 3600.0

    @test_nowarn gw2river!(
        imax, jmax, 1, imax, 1, jmax, nzg,
        slz, Δt, soiltxt, landmask, wtd, maxdepth,
        riverdepth, width, length, area, fdepth, qrf
    )

    # 检查结果
    @test all(isfinite.(qrf))
    # 由于水位高于河床，应该有向河流的排水
    @test any(qrf .> 0.0)
end
