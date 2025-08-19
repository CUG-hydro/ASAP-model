using ASAP, Test


@testset "Kinematic Wave Flood Routing" begin
    # 创建简单的河流网络测试
    imax, jmax = 5, 5

    # 流向矩阵(简单的从左到右流动)
    fd = [1 1 1 1 0;
        1 1 1 1 0;
        1 1 1 1 0;
        1 1 1 1 0;
        0 0 0 0 0]

    bfd = zeros(Int, size(fd))  # 无反向流

    # 河流参数
    qnew = [10.0 8.0 6.0 4.0 0.0;
        12.0 10.0 8.0 6.0 0.0;
        8.0 6.0 4.0 2.0 0.0;
        5.0 4.0 3.0 2.0 0.0;
        0.0 0.0 0.0 0.0 0.0]

    qs = zeros(size(fd))
    qrf = zeros(size(fd))
    delsfcwat = zeros(size(fd))
    slope = fill(0.001, size(fd))  # 1‰坡度
    depth = ones(size(fd)) * 2.0
    width = ones(size(fd)) * 10.0
    length = ones(size(fd)) * 1000.0
    maxdepth = ones(size(fd)) * 3.0
    area = ones(size(fd)) * 1e6  # 1 km²
    riverarea = width .* length
    floodarea = area .- riverarea
    riverchannel = maxdepth .* riverarea
    qmean = zeros(size(fd))
    floodheight = zeros(size(fd))
    topo = zeros(size(fd))

    Δt = 3600.0  # 1小时
    dtlr = Δt / 4.0  # 子时间步

    # 测试运动波洪水路由
    qnew_old = copy(qnew)
    @test_nowarn rivers_kw_flood!(
        imax, jmax, 1, imax, 1, jmax,
        Δt, dtlr, fd, bfd, qnew, qs, qrf, delsfcwat,
        slope, depth, width, length, maxdepth, area,
        riverarea, floodarea, riverchannel, qmean, floodheight, topo
    )

    # 检查结果合理性
    @test all(isfinite.(qnew))
    @test all(qnew .>= 0.0)  # 流量非负
    @test all(isfinite.(depth))
    @test all(depth .>= 0.0)  # 深度非负
end


@testset "Flooding Calculation" begin
    # 测试洪泛区漫流
    imax, jmax = 4, 4

    fd = [1 1 1 0;
        1 1 1 0;
        1 1 1 0;
        0 0 0 0]

    bfd = zeros(Int, size(fd))

    # 地形(中心较低)
    topo = [100.0 99.5 99.0 98.5;
        99.5 98.0 97.5 98.0;
        99.0 97.5 98.0 98.5;
        98.5 98.0 98.5 99.0]

    area = ones(size(fd)) * 1e6
    riverwidth = ones(size(fd)) * 10.0
    riverlength = ones(size(fd)) * 1000.0
    riverdepth = ones(size(fd)) * 2.0

    # 初始洪水高度(中心较高)
    floodheight = [0.1 0.2 0.1 0.0;
        0.2 0.5 0.3 0.1;
        0.1 0.3 0.2 0.1;
        0.0 0.1 0.1 0.0]

    delsfcwat = zeros(size(fd))

    Δt = 3600.0

    floodheight_old = copy(floodheight)
    @test_nowarn flooding!(
        imax, jmax, 1, imax, 1, jmax, Δt,
        fd, bfd, topo, area, riverwidth, riverlength,
        riverdepth, floodheight, delsfcwat
    )

    @test all(isfinite.(floodheight))
    @test all(floodheight .>= 0.0)
    @test all(isfinite.(delsfcwat))

    # 验证水量保守性(总的地表水变化应该接近0)
    @test abs(sum(delsfcwat)) < 1e-10
end

@testset "MoveQRF Function" begin
    # 测试小河流通量重分配
    imax = 4
    js, je = 1, 3

    fd = [1 1 0; 1 1 0; 1 1 0; 0 0 0]

    qrf = [0.1 0.05 0.0;
        0.08 0.03 0.0;
        0.06 0.02 0.0;
        0.0 0.0 0.0]

    area = ones(size(fd)) * 1e6

    # 小河流宽度(< 1m)
    width = [0.5 0.8 1.5;
        0.6 0.7 1.2;
        0.9 1.1 1.0;
        1.0 1.0 1.0]

    qrf_old = copy(qrf)
    @test_nowarn moveqrf!(imax, js, je, fd, qrf, area, width)

    @test all(isfinite.(qrf))

    # 检查小河流的通量已被重新分配
    for j in (js+1):(je-1)
        for i in 2:(imax-1)
            if fd[i, j] > 0 && width[i, j] < 1.0
                @test qrf[i, j] == 0.0  # 小河流通量应为0
            end
        end
    end
end
