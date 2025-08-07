using Test, ASAP


@testset "Diffusive Wave Flood Routing" begin
    # 类似的测试但使用扩散波方程
    imax, jmax = 3, 4

    fd = [1 1 0;
        1 1 0;
        1 1 0;
        0 0 0]

    bfd = zeros(Int, size(fd))

    qnew = [5.0 4.0 0.0;
        6.0 5.0 0.0;
        4.0 3.0 0.0;
        0.0 0.0 0.0]

    qs = zeros(size(fd))
    qrf = zeros(size(fd))
    delsfcwat = zeros(size(fd))
    slope = fill(0.001, size(fd))
    depth = ones(size(fd)) * 1.5
    width = ones(size(fd)) * 8.0
    length = ones(size(fd)) * 800.0
    maxdepth = ones(size(fd)) * 2.5
    area = ones(size(fd)) * 6.4e5
    riverarea = width .* length
    floodarea = area .- riverarea
    riverchannel = maxdepth .* riverarea
    qmean = zeros(size(fd))
    floodheight = zeros(size(fd))
    topo = zeros(size(fd))

    Δt = 3600.0
    dtlr = Δt / 4.0

    @test_nowarn rivers_dw_flood!(
        imax, 1, jmax, Δt, dtlr, fd, bfd, qnew,
        qs, qrf, delsfcwat, slope, depth, width, length,
        maxdepth, area, riverarea, floodarea, riverchannel,
        qmean, floodheight, topo
    )

    @test all(isfinite.(qnew))
    @test all(qnew .>= 0.0)
end
