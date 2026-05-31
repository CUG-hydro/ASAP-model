using ASAP, Test

# rootdepth_main使用nzg=40（标准层数，与icefactor[26:40]索引匹配）
const NZG_RD = 40

# 构建rootdepth_main所需的最小参数集（网格大小3×3，内层循环只访问(2,2)）
function make_rootdepth_inputs(; freedrain=1, veg_val=1.0, landmask_val=1, nzg=NZG_RD)
    is, ie, js, je = 1, 3, 1, 3
    nx, ny = ie, je   # 3×3 grid
    z₋ₕ, dz = initializesoildepth(nzg)
    Δt = 3600.0

    landmask        = zeros(Int,     nx, ny);  landmask[2, 2] = landmask_val
    veg             = fill(veg_val,  nx, ny)
    hveg            = fill(1.0,      nx, ny)
    soiltxt         = ones(Int,   2, nx, ny)
    wind            = fill(3.0,      nx, ny)
    temp            = fill(293.0,    nx, ny)
    qair            = fill(0.01,     nx, ny)
    press           = fill(101325.0, nx, ny)
    netrad          = fill(100.0,    nx, ny)
    rshort          = fill(150.0,    nx, ny)
    lai             = fill(2.0,      nx, ny)
    precip          = fill(5.0,      nx, ny)
    qsrun           = zeros(Float64, nx, ny)
    θ               = fill(0.25,  nzg, nx, ny)
    θ_eq            = fill(0.20,  nzg, nx, ny)
    θ_wtd           = fill(0.35,     nx, ny)
    wtd             = fill(-1.0,     nx, ny)
    waterdeficit    = zeros(Float64, nx, ny)
    watext          = zeros(Float64, nzg, nx, ny)
    watextdeep      = zeros(Float64, nx, ny)
    rech            = zeros(Float64, nx, ny)
    deeprech        = zeros(Float64, nx, ny)
    et_s            = zeros(Float64, nx, ny)
    et_i            = zeros(Float64, nx, ny)
    et_c            = zeros(Float64, nx, ny)
    intercepstore   = zeros(Float64, nx, ny)
    ppacum          = zeros(Float64, nx, ny)
    pInfiltDepth    = zeros(Float64, nx, ny)
    pInfiltDepthK_old = zeros(Int8,  nx, ny)
    qlat            = zeros(Float64, nx, ny)
    qlatsum         = zeros(Float64, nx, ny)
    qsprings        = zeros(Float64, nx, ny)
    inactivedays    = zeros(Int,   nzg+2, nx, ny)
    maxinactivedays = 10
    fdepth          = fill(2.0,      nx, ny)
    steps           = 1.0
    floodheight     = zeros(Float64, nx, ny)
    qrf             = zeros(Float64, nx, ny)
    delsfcwat       = zeros(Float64, nx, ny)
    icefactor       = zeros(Int8,    nx, ny, nzg)
    wtdflux         = zeros(Float64, nx, ny)
    et_s_daily      = zeros(Float64, nx, ny)
    et_c_daily      = zeros(Float64, nx, ny)
    transptop       = zeros(Float64, nx, ny)
    infilk          = zeros(Int8,    nx, ny)
    infilflux       = zeros(Float64, nzg, nx, ny)
    infilfluxday    = zeros(Float64, nzg, nx, ny)
    infilcounter    = zeros(Int16,   nzg, nx, ny)
    hour            = 0
    o18             = zeros(Float64, nzg, nx, ny)
    o18ratiopp      = zeros(Float64, nx, ny)
    tempsfc         = fill(293.0,    nx, ny)
    qlato18         = zeros(Float64, nx, ny)
    transpo18       = zeros(Float64, nx, ny)
    upflux          = zeros(Float64, nzg, nx, ny)

    return (freedrain, is, ie, js, je, nzg,
            z₋ₕ, dz, Δt,
            landmask, veg, hveg,
            soiltxt, wind, temp,
            qair, press, netrad,
            rshort, lai, precip,
            qsrun, θ, θ_eq,
            θ_wtd, wtd, waterdeficit,
            watext,
            watextdeep, rech,
            deeprech,
            et_s, et_i, et_c,
            intercepstore, ppacum,
            pInfiltDepth, pInfiltDepthK_old, qlat,
            qlatsum, qsprings, inactivedays,
            maxinactivedays, fdepth,
            steps, floodheight, qrf,
            delsfcwat, icefactor, wtdflux,
            et_s_daily, et_c_daily, transptop,
            infilk, infilflux, infilfluxday,
            infilcounter, hour,
            o18, o18ratiopp, tempsfc,
            qlato18, transpo18, upflux)
end

@testset "rootdepth_main: 全零陆地掩码（无活跃格点）" begin
    # landmask全零：内层循环立即跳过，函数应返回nothing且不修改任何数组
    args = make_rootdepth_inputs(landmask_val=0)
    et_s_before = copy(args[33])   # et_s is the 33rd arg
    result = rootdepth_main(args...)
    @test result === nothing
    # 所有输出应保持不变
    @test args[33] == et_s_before
end

@testset "rootdepth_main: 水体格点（veg=1）" begin
    # veg=1 (水体)：执行潜在蒸散发计算，但跳过植被相关代码
    args = make_rootdepth_inputs(veg_val=1.0, landmask_val=1)
    et_s = args[33]
    et_s_init = et_s[2, 2]

    result = rootdepth_main(args...)
    @test result === nothing
    # 水体格点应更新土壤蒸发
    @test et_s[2, 2] >= et_s_init
end

@testset "rootdepth_main: 自由排水+植被格点（veg=7，短草）" begin
    # freedrain=1（自由排水），veg=7（短草）：执行全部子模块（extraction、soilfluxes、updatewtd_shallow）
    args = make_rootdepth_inputs(freedrain=1, veg_val=7.0, landmask_val=1)
    et_s    = args[33]
    et_c    = args[35]
    rech    = args[30]

    result = rootdepth_main(args...)
    @test result === nothing
    # 植被格点应更新蒸腾和补给
    @test isfinite(et_s[2, 2])
    @test isfinite(et_c[2, 2])
    @test isfinite(rech[2, 2])
end

@testset "rootdepth_main: 非自由排水+植被格点（veg=7）" begin
    # freedrain=0（非自由排水）：测试soilfluxes的!freedrain分支在全流程中的集成
    args = make_rootdepth_inputs(freedrain=0, veg_val=7.0, landmask_val=1)
    et_s = args[33]
    et_c = args[35]

    result = rootdepth_main(args...)
    @test result === nothing
    @test isfinite(et_s[2, 2])
    @test isfinite(et_c[2, 2])
end
