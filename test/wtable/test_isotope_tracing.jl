using ASAP, Test

@testset "Lateral Isotope Transport" begin
    # 创建简单的网格用于同位素传输测试
    imax, jmax = 4, 4
    nzg = 3

    slz = [-0.1, -0.5, -1.0, -2.0]  # nzg+1层边界

    # 土壤类型
    soiltxt = ones(Int, 2, imax, jmax)

    # 水位分布(创建梯度)
    wtd = [-0.8 -0.9 -1.0 -1.1;
        -0.7 -0.8 -0.9 -1.0;
        -0.6 -0.7 -0.8 -0.9;
        -0.5 -0.6 -0.7 -0.8]

    qlat = zeros(size(wtd))
    fdepth = ones(size(wtd)) * 2.0
    topo = zeros(size(wtd))
    landmask = ones(Int, size(wtd))
    area = ones(size(wtd)) * 1e6
    lats = fill(30.0, size(wtd))  # 纬度30度
    dxy = 0.01  # 网格间距

    # 土壤湿度(假设均匀分布)
    smoi = fill(0.3, nzg, imax, jmax)

    # 同位素浓度(创建梯度)
    o18 = zeros(nzg, imax, jmax)
    for k in 1:nzg
        for j in 1:jmax
            for i in 1:imax
                o18[k, i, j] = -8.0 + i * 0.5 + j * 0.3  # δ18O值
            end
        end
    end

    # 输出数组
    qlato18 = zeros(size(wtd))
    qlatin = zeros(size(wtd))
    qlatout = zeros(size(wtd))
    qlatino18 = zeros(size(wtd))
    qlatouto18 = zeros(size(wtd))
    qlatinsum = zeros(size(wtd))
    qlatoutsum = zeros(size(wtd))
    qlatino18sum = zeros(size(wtd))
    qlatouto18sum = zeros(size(wtd))

    Δt = 3600.0

    @test_nowarn lateral_isotope!(
        imax, jmax, 1, imax, 1, jmax, nzg,
        soiltxt, wtd, qlat, fdepth, topo, landmask,
        Δt, area, lats, dxy, slz, o18, smoi,
        qlato18, qlatin, qlatout, qlatino18, qlatouto18,
        qlatinsum, qlatoutsum, qlatino18sum, qlatouto18sum
    )

    # 检查结果
    @test all(isfinite.(qlat))
    @test all(isfinite.(qlato18))
    @test all(isfinite.(qlatin))
    @test all(isfinite.(qlatout))
    @test all(isfinite.(qlatino18))
    @test all(isfinite.(qlatouto18))

    # 检查累积量是否合理
    @test all(isfinite.(qlatinsum))
    @test all(isfinite.(qlatoutsum))
    @test all(isfinite.(qlatino18sum))
    @test all(isfinite.(qlatouto18sum))

    # 由于有水位梯度，应该有侧向流
    @test any(abs.(qlat) .> 1e-10)

    # 由于有同位素浓度，应该有同位素传输
    @test any(abs.(qlato18) .> 1e-15)
end

@testset "Deep Water Table Update with Isotope" begin
    # 测试深层水位更新功能
    imax, jmax = 3, 3
    nzg = 2

    slz = [-0.1, -1.0, -3.0]
    dz = [0.9, 2.0]

    soiltxt = ones(Int, 2, imax, jmax)
    wtd = fill(-2.0, imax, jmax)  # 深水位
    bottomflux = zeros(imax, jmax)
    rech = zeros(imax, jmax)
    qslat = fill(0.001, imax, jmax)  # 小的侧向流入
    qlat = fill(0.0015, imax, jmax)  # 略大的侧向流出
    landmask = ones(Int, imax, jmax)
    smoi = fill(0.25, nzg, imax, jmax)
    smoieq = fill(0.15, nzg, imax, jmax)
    smoiwtd = fill(0.2, imax, jmax)
    qsprings = zeros(imax, jmax)

    Δt = 3600.0

    # 保存初始值
    wtd_initial = copy(wtd)
    rech_initial = copy(rech)

    @test_nowarn updatedeepwtable!(
        imax, jmax, 1, jmax, nzg, slz, dz, soiltxt,
        wtd, bottomflux, rech, qslat, qlat, landmask,
        Δt, smoi, smoieq, smoiwtd, qsprings
    )

    # 检查结果
    @test all(isfinite.(wtd))
    @test all(isfinite.(rech))
    @test all(isfinite.(qsprings))
    @test all(isfinite.(smoiwtd))

    # 由于有净流出(qlat > qslat)，水位应该下降
    # @test all(wtd .<= wtd_initial)

    # 补给量应该发生变化
    @test any(rech .!= rech_initial)
end

@testset "Isotope Conservation" begin
    # 测试同位素质量守恒
    imax, jmax = 3, 3
    nzg = 2

    # 简单的均匀场用于测试守恒性
    slz = [-0.1, -0.6, -1.5]
    soiltxt = ones(Int, 2, imax, jmax)
    wtd = fill(-0.8, imax, jmax)
    fdepth = ones(size(wtd)) * 2.0
    topo = zeros(size(wtd))
    landmask = ones(Int, size(wtd))
    area = ones(size(wtd)) * 1e6
    lats = fill(45.0, size(wtd))
    dxy = 0.01

    smoi = fill(0.3, nzg, imax, jmax)

    # 均匀的同位素浓度
    o18 = fill(-10.0, nzg, imax, jmax)

    # 由于浓度均匀且地形平坦，同位素净传输应该很小
    qlato18 = zeros(size(wtd))
    qlatin = zeros(size(wtd))
    qlatout = zeros(size(wtd))
    qlatino18 = zeros(size(wtd))
    qlatouto18 = zeros(size(wtd))
    qlatinsum = zeros(size(wtd))
    qlatoutsum = zeros(size(wtd))
    qlatino18sum = zeros(size(wtd))
    qlatouto18sum = zeros(size(wtd))

    qlat = zeros(size(wtd))

    Δt = 3600.0

    lateral_isotope!(
        imax, jmax, 1, imax, 1, jmax, nzg,
        soiltxt, wtd, qlat, fdepth, topo, landmask,
        Δt, area, lats, dxy, slz, o18, smoi,
        qlato18, qlatin, qlatout, qlatino18, qlatouto18,
        qlatinsum, qlatoutsum, qlatino18sum, qlatouto18sum
    )

    # 由于场均匀，侧向流应该很小
    @test maximum(abs.(qlat)) < 1e-8

    # 同位素传输也应该很小
    @test maximum(abs.(qlato18)) < 1e-12

    # 质量守恒检查：入流 + 出流 ≈ 净变化
    for i in 2:(imax-1)
        for j in 2:(jmax-1)
            if landmask[i, j] == 1
                net_flow = qlatin[i, j] + qlatout[i, j]
                @test abs(net_flow - qlat[i, j]) < 1e-10

                net_isotope = qlatino18[i, j] + qlatouto18[i, j]
                @test abs(net_isotope - qlato18[i, j]) < 1e-12
            end
        end
    end
end




@testset "Isotope Gradient Transport" begin
    # 测试同位素浓度梯度传输
    imax, jmax = 3, 3
    nzg = 1

    slz = [-0.1, -2.0]
    soiltxt = ones(Int, 2, imax, jmax)

    # 平坦地形但有水位梯度
    wtd = [-0.5 -0.8 -1.1;
        -0.5 -0.8 -1.1;
        -0.5 -0.8 -1.1]

    fdepth = ones(size(wtd)) * 3.0
    topo = zeros(size(wtd))
    landmask = ones(Int, size(wtd))
    area = ones(size(wtd)) * 1e6
    lats = fill(30.0, size(wtd))
    dxy = 0.01

    smoi = fill(0.3, nzg, imax, jmax)

    # 创建同位素浓度梯度 (左到右递增)
    o18 = zeros(nzg, imax, jmax)
    for j in 1:jmax
        for i in 1:imax
            o18[1, i, j] = -12.0 + (i - 1) * 2.0  # -12, -10, -8
        end
    end

    qlat = zeros(size(wtd))
    qlato18 = zeros(size(wtd))
    qlatin = zeros(size(wtd))
    qlatout = zeros(size(wtd))
    qlatino18 = zeros(size(wtd))
    qlatouto18 = zeros(size(wtd))
    qlatinsum = zeros(size(wtd))
    qlatoutsum = zeros(size(wtd))
    qlatino18sum = zeros(size(wtd))
    qlatouto18sum = zeros(size(wtd))

    Δt = 3600.0

    lateral_isotope!(
        imax, jmax, 1, imax, 1, jmax, nzg,
        soiltxt, wtd, qlat, fdepth, topo, landmask,
        Δt, area, lats, dxy, slz, o18, smoi,
        qlato18, qlatin, qlatout, qlatino18, qlatouto18,
        qlatinsum, qlatoutsum, qlatino18sum, qlatouto18sum
    )

    # 由于水位梯度，应该有从左到右的流动
    @test any(qlat .> 1e-10)

    # 由于同位素浓度梯度，应该有显著的同位素传输
    @test any(abs.(qlato18) .> 1e-12)

    # 中间列应该有净同位素流入（因为从高浓度区域流入）
    @test qlato18[2, 2] != 0.0
end
