using ASAP, Test


@testset "Flow Direction Function" begin
    # 测试流向计算
    fdir = [1 2 4; 8 16 32; 64 128 0]

    # 测试各个方向
    @test flowdir(fdir, 1, 1) == (2, 1)  # 向东
    @test flowdir(fdir, 1, 2) == (2, 1)  # 向东南
    @test flowdir(fdir, 1, 3) == (1, 2)  # 向南
    @test flowdir(fdir, 2, 2) == (1, 2)  # 向西
    @test flowdir(fdir, 3, 1) == (3, 2)  # 向北
    @test flowdir(fdir, 3, 2) == (4, 3)  # 向东北
    @test flowdir(fdir, 2, 1) == (1, 0)  # 8, 向西南
    @test flowdir(fdir, 2, 3) == (1, 4)  # 32, 向西北
end


@testset "Lateral Flow Calculation" begin
    # 简单的3x3网格测试侧向流
    imax, jmax = 3, 3

    wtd = [-1.0 -1.2 -1.4;
        -1.1 -1.3 -1.5;
        -1.2 -1.4 -1.6]

    qlat = zeros(size(wtd))
    fdepth = ones(size(wtd)) * 2.0
    topo = zeros(size(wtd))
    landmask = ones(Int, size(wtd))
    area = ones(size(wtd)) * 1000.0
    κlat = ones(size(wtd)) * 1e-5

    Δt = 3600.0

    @test_nowarn lateral_flow!(
        imax, jmax, 1, imax, 1, jmax,
        wtd, qlat, fdepth, topo, landmask, Δt, area, κlat
    )

    # 检查侧向流的基本性质
    @test all(isfinite.(qlat))
    # 中心点应该有流出(因为水位梯度)
    @test qlat[2, 2] != 0.0
end



@testset "Basic Water Table Calculations" begin
    # 创建简单的测试网格
    imax, jmax = 5, 5
    nzg = 3

    # 土壤层设置
    slz = [-0.1, -0.5, -1.0, -2.0]  # nzg+1 = 4 layers
    dz = [0.4, 0.5, 1.0]  # nzg = 3 thickness

    # 创建测试数组
    area = ones(imax, jmax) * 1000.0  # 1000 m²
    soiltxt = ones(Int, 2, imax, jmax)  # 土壤类型1
    wtd = fill(-1.5, imax, jmax)  # 水位深度
    bottomflux = zeros(imax, jmax)
    rech = zeros(imax, jmax)
    qslat = zeros(imax, jmax)
    fdepth = ones(imax, jmax) * 2.0
    topo = zeros(imax, jmax)
    landmask = ones(Int, imax, jmax)
    smoi = fill(0.3, nzg, imax, jmax)
    smoieq = fill(0.2, nzg, imax, jmax)
    smoiwtd = fill(0.25, imax, jmax)
    qsprings = zeros(imax, jmax)

    Δt = 3600.0  # 1小时

    # # 测试wtable!函数运行不出错
    # @test_nowarn wtable!(
    #     imax, jmax, 1, imax, 1, jmax, nzg,
    #     slz, dz, area, soiltxt, wtd, bottomflux,
    #     rech, qslat, fdepth, topo, landmask, Δt,
    #     smoi, smoieq, smoiwtd, qsprings
    # )

    # # 检查结果的基本属性
    # @test all(isfinite.(wtd))
    # @test all(isfinite.(qsprings))
    # @test all(isfinite.(rech))
end
