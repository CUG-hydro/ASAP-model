using Test
using ASAP.WaterTable

@testset "Diffusion Wave Flood Routing - flood scenario" begin
    # Grid and flow setup
    imax, jmax = 5, 5
    js, je = 1, jmax

    # Flow direction: eastward in interior, 0 at last column/boundaries
    fd = [1 1 1 1 0;
          1 1 1 1 0;
          1 1 1 1 0;
          1 1 1 1 0;
          0 0 0 0 0]
    bfd = zeros(Int, size(fd))

    # Fields
    qnew = [10.0 8.0 6.0 4.0 0.0;
            12.0 10.0 8.0 6.0 0.0;
            8.0 6.0 4.0 2.0 0.0;
            5.0 4.0 3.0 2.0 0.0;
            0.0 0.0 0.0 0.0 0.0]

    qs = zeros(size(fd))
    qrf = zeros(size(fd))
    delsfcwat = zeros(size(fd))
    slope = fill(0.001, size(fd))
    depth = ones(size(fd)) * 2.0
    width = ones(size(fd)) * 10.0
    length = ones(size(fd)) * 1000.0
    maxdepth = ones(size(fd)) * 3.0
    area = ones(size(fd)) * 1e6
    riverarea = width .* length
    floodarea = area .- riverarea
    riverchannel = maxdepth .* riverarea
    qmean = zeros(size(fd))

    # Create flood conditions to trigger diffusion wave branch
    floodheight = zeros(size(fd))
    floodheight[2, 2] = 0.1
    floodheight[3, 3] = 0.2

    # Also test width==0 path (uses sqrt(area) as flow width)
    width[2, 2] = 0.0

    topo = zeros(size(fd))

    Δt = 3600.0
    dtlr = Δt / 4.0

    qnew_before = copy(qnew)
    qmean_before = copy(qmean)

    @test_nowarn rivers_dw_flood!(
        imax, js, je, Δt, dtlr,
        fd, bfd, qnew, qs, qrf, delsfcwat,
        slope, depth, width, length, maxdepth, area,
        riverarea, floodarea, riverchannel, qmean, floodheight, topo
    )

    # Basic sanity
    @test all(isfinite.(qnew))
    # @test all(qnew .>= 0.0)
    @test all(isfinite.(depth)) && all(depth .>= 0.0)
    @test qnew != qnew_before  # some change expected

    # qmean accumulated correctly
    @test qmean ≈ qmean_before .+ qnew .* dtlr atol=1e-8

    # Cells with fd==0 should have zero discharge
    for j in js:je, i in 1:imax
        if fd[i, j] == 0
            @test qnew[i, j] == 0.0
        end
    end
end

@testset "Diffusion Wave Flood Routing - non-flood scenario" begin
    imax, jmax = 5, 5
    js, je = 1, jmax

    fd = [1 1 1 1 0;
          1 1 1 1 0;
          1 1 1 1 0;
          1 1 1 1 0;
          0 0 0 0 0]
    bfd = zeros(Int, size(fd))

    qnew = [2.0 2.0 2.0 2.0 0.0;
            2.0 2.0 2.0 2.0 0.0;
            2.0 2.0 2.0 2.0 0.0;
            2.0 2.0 2.0 2.0 0.0;
            0.0 0.0 0.0 0.0 0.0]

    qs = zeros(size(fd))
    qrf = zeros(size(fd))
    delsfcwat = zeros(size(fd))
    slope = fill(0.001, size(fd))
    depth = ones(size(fd)) * 1.0
    width = ones(size(fd)) * 5.0
    length = ones(size(fd)) * 500.0
    maxdepth = ones(size(fd)) * 3.0
    area = ones(size(fd)) * 1e6
    riverarea = width .* length
    floodarea = area .- riverarea
    riverchannel = maxdepth .* riverarea
    qmean = zeros(size(fd))
    floodheight = zeros(size(fd))  # no flood
    topo = zeros(size(fd))

    Δt = 3600.0
    dtlr = Δt / 6.0

    @test_nowarn rivers_dw_flood!(
        imax, js, je, Δt, dtlr,
        fd, bfd, qnew, qs, qrf, delsfcwat,
        slope, depth, width, length, maxdepth, area,
        riverarea, floodarea, riverchannel, qmean, floodheight, topo
    )

    @test all(isfinite.(qnew))
    @test all(qnew .>= 0.0)
    @test all(isfinite.(depth)) && all(depth .>= 0.0)
end
