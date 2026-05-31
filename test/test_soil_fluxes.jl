using ASAP, Test

# 测试基本功能
@testset "soilfluxes基本功能" begin
    nzg = 10
    dt = 3600.0    # 1小时时间步长
    soiltxt = 1    # 沙土

    z₋ₕ, Δz = initializesoildepth(nzg)

    θ_wtd = 0.3
    transp = zeros(Float64, nzg)
    transpdeep = 0.0
    θ = fill(0.25, nzg)
    wtd = -2.0
    precip = 5.0
    pet_s = 2.0
    fdepth = 2.0
    qlat = 0.0
    qrf = 0.0
    flood = 0.0
    icefactor = zeros(Int8, nzg)

    et_s, runoff, rech, flux, qrfcorrect, θ_new =
        soilfluxes(nzg, dt, z₋ₕ, Δz, soiltxt,
            θ_wtd, transp, transpdeep, copy(θ), wtd,
            precip, pet_s, fdepth, qlat, qrf, flood, icefactor)

    @test et_s >= 0.0
    @test runoff >= 0.0
    @test length(θ_new) == nzg
    @test length(flux) == nzg + 1
    soil = get_soil_params(soiltxt)
    for k in 1:nzg
        @test θ_new[k] >= 0.0
        @test θ_new[k] <= soil.θ_sat * 1.01
    end
end

# 测试不同降水条件
@testset "不同降水条件测试" begin
    nzg = 5
    z₋ₕ, Δz = initializesoildepth(nzg)
    dt = 3600.0
    soiltxt = 1
    θ_wtd = 0.3
    transp = zeros(Float64, nzg)
    transpdeep = 0.0
    θ = fill(0.25, nzg)
    wtd = -2.0
    pet_s = 2.0
    fdepth = 2.0
    icefactor = zeros(Int8, nzg)

    # 无降水
    _, runoff1, _, _, _, _ =
        soilfluxes(nzg, dt, z₋ₕ, Δz, soiltxt,
            θ_wtd, transp, transpdeep, copy(θ), wtd,
            0.0, pet_s, fdepth, 0.0, 0.0, 0.0, icefactor)

    # 大降水
    _, runoff2, _, _, _, _ =
        soilfluxes(nzg, dt, z₋ₕ, Δz, soiltxt,
            θ_wtd, transp, transpdeep, copy(θ), wtd,
            50.0, pet_s, fdepth, 0.0, 0.0, 0.0, icefactor)

    # 大降水应该产生更多径流
    @test runoff2 >= runoff1
end

# 测试自由排水条件
@testset "自由排水条件" begin
    nzg = 5
    z₋ₕ, Δz = initializesoildepth(nzg)
    dt = 3600.0
    soiltxt = 1
    θ_wtd = 0.3
    transp = zeros(Float64, nzg)
    θ = fill(0.25, nzg)
    wtd = -2.0
    icefactor = zeros(Int8, nzg)

    et_s, runoff, rech, flux, qrfcorrect, θ_new =
        soilfluxes(nzg, dt, z₋ₕ, Δz, soiltxt,
            θ_wtd, transp, 0.0, copy(θ), wtd,
            10.0, 2.0, 2.0, 0.0, 0.0, 0.0, icefactor;
            freedrain=true)

    @test et_s >= 0.0
    @test runoff >= 0.0
    @test length(θ_new) == nzg
end

# 测试冰冻因子影响
@testset "冰冻因子测试" begin
    nzg = 5
    z₋ₕ, Δz = initializesoildepth(nzg)
    dt = 3600.0
    soiltxt = 1
    θ_wtd = 0.3
    transp = zeros(Float64, nzg)
    θ = fill(0.25, nzg)
    wtd = -2.0
    precip = 10.0

    # 无冰冻
    icefactor_none = zeros(Int8, nzg)
    _, runoff1, _, _, _, _ =
        soilfluxes(nzg, dt, z₋ₕ, Δz, soiltxt,
            θ_wtd, transp, 0.0, copy(θ), wtd,
            precip, 2.0, 2.0, 0.0, 0.0, 0.0, icefactor_none)

    # 全部冰冻
    icefactor_all = ones(Int8, nzg)
    _, runoff2, _, _, _, _ =
        soilfluxes(nzg, dt, z₋ₕ, Δz, soiltxt,
            θ_wtd, transp, 0.0, copy(θ), wtd,
            precip, 2.0, 2.0, 0.0, 0.0, 0.0, icefactor_all)

    # 冰冻会阻碍水分运动，导致更多径流
    @test runoff2 >= runoff1
end

# 测试饱和土壤产生径流
@testset "饱和土壤径流" begin
    nzg = 5
    z₋ₕ, Δz = initializesoildepth(nzg)
    # 使用粘土 (type 11, Ksat=0.000001283 m/s ≈ 4.6mm/hr) — 50mm降水大幅超过入渗能力
    soiltxt = 11
    soil = get_soil_params(soiltxt)
    icefactor = zeros(Int8, nzg)

    θ_sat = fill(soil.θ_sat, nzg)

    _, runoff, _, _, _, _ =
        soilfluxes(nzg, 3600.0, z₋ₕ, Δz, soiltxt,
            soil.θ_sat, zeros(nzg), 0.0, θ_sat, -2.0,
            50.0, 2.0, 2.0, 0.0, 0.0, 0.0, icefactor)

    @test runoff > 0.0
end

# 测试非自由排水条件：水位在土柱以下（jwt <= 2）
@testset "非自由排水-水位在土柱以下" begin
    nzg = 10
    z₋ₕ, Δz = initializesoildepth(nzg)
    dt = 3600.0
    soiltxt = 4    # 壤土
    θ_wtd = 0.4
    transp = zeros(Float64, nzg)
    θ = fill(0.25, nzg)
    # wtd低于最底层边界，jwt=1 → jwt<=2分支
    wtd = z₋ₕ[1] - 0.5
    icefactor = zeros(Int8, nzg)

    et_s, runoff, rech, flux, qrfcorrect, θ_new =
        soilfluxes(nzg, dt, z₋ₕ, Δz, soiltxt,
            θ_wtd, transp, 0.0, copy(θ), wtd,
            5.0, 1.0, 2.0, 0.0, 0.0, 0.0, icefactor;
            freedrain=false)

    @test et_s >= 0.0
    @test runoff >= 0.0
    @test length(θ_new) == nzg
    @test length(flux) == nzg + 1
    @test all(isfinite.(θ_new))
    # 非自由排水且水位极低，底层无渗漏
    @test flux[1] == 0.0
end

# 测试非自由排水条件：水位在中间层（jwt > 3）
@testset "非自由排水-水位在中间层" begin
    nzg = 10
    z₋ₕ, Δz = initializesoildepth(nzg)
    dt = 3600.0
    soiltxt = 4    # 壤土
    soil = get_soil_params(soiltxt)
    θ_wtd = soil.θ_sat
    transp = zeros(Float64, nzg)
    # 土柱接近饱和，水位在层5和层6之间（wtd=-0.55 → jwt=6）
    θ = fill(soil.θ_sat * 0.9, nzg)
    wtd = 0.5 * (z₋ₕ[5] + z₋ₕ[6])  # between layers 5 and 6
    icefactor = zeros(Int8, nzg)

    et_s, runoff, rech, flux, qrfcorrect, θ_new =
        soilfluxes(nzg, dt, z₋ₕ, Δz, soiltxt,
            θ_wtd, transp, 0.0, copy(θ), wtd,
            3.0, 1.0, 2.0, 0.0, 0.0, 0.0, icefactor;
            freedrain=false)

    @test et_s >= 0.0
    @test runoff >= 0.0
    @test length(θ_new) == nzg
    @test all(isfinite.(θ_new))
    # 非自由排水，底层通量为零
    @test flux[1] == 0.0
end

# 对比自由排水与非自由排水的底层渗漏差异
@testset "自由排水与非自由排水底层渗漏对比" begin
    nzg = 10
    z₋ₕ, Δz = initializesoildepth(nzg)
    dt = 3600.0
    soiltxt = 1    # 沙土（Ksat大，容易产生底层渗漏）
    θ_wtd = 0.35
    transp = zeros(Float64, nzg)
    θ = fill(0.3, nzg)
    wtd = -2.0
    icefactor = zeros(Int8, nzg)

    # 自由排水
    _, _, rech_free, flux_free, _, _ =
        soilfluxes(nzg, dt, z₋ₕ, Δz, soiltxt,
            θ_wtd, transp, 0.0, copy(θ), wtd,
            5.0, 1.0, 2.0, 0.0, 0.0, 0.0, icefactor;
            freedrain=true)

    # 非自由排水
    _, _, rech_nonfree, flux_nonfree, _, _ =
        soilfluxes(nzg, dt, z₋ₕ, Δz, soiltxt,
            θ_wtd, transp, 0.0, copy(θ), wtd,
            5.0, 1.0, 2.0, 0.0, 0.0, 0.0, icefactor;
            freedrain=false)

    # 非自由排水时底层通量为零（无重力排水）
    @test flux_nonfree[1] == 0.0
    @test rech_nonfree == 0.0
    # 自由排水时底层有向下渗漏（Q[1]<0 表示向下排水）
    @test rech_free <= 0.0
end

# 测试三对角矩阵求解器
@testset "tridag! 函数测试" begin
    # 测试简单3x3系统: [2 1 0; 1 2 1; 0 1 2] * x = [5; 8; 5]  → x = [1; 2; 1]+adjust
    # 使用一个有解析解的简单系统
    n = 3
    a = [0.0, -1.0, -1.0]  # 下对角线 (a[1]未使用)
    b = [2.0, 2.0, 2.0]    # 主对角线
    c = [-1.0, -1.0, 0.0]  # 上对角线 (c[n]未使用)
    r = [1.0, 2.0, 3.0]    # 右端向量
    u = zeros(Float64, n)

    # 方程组: 2u1 - u2 = 1; -u1 + 2u2 - u3 = 2; -u2 + 2u3 = 3
    # 解: u1 = 2.5, u2 = 4.0, u3 = 3.5
    ASAP.tridag!(a, b, c, r, u)

    @test u[1] ≈ 2.5 atol=1e-10
    @test u[2] ≈ 4.0 atol=1e-10
    @test u[3] ≈ 3.5 atol=1e-10

    # 测试更大的系统
    n2 = 5
    a2 = [0.0, 1.0, 1.0, 1.0, 1.0]
    b2 = [3.0, 3.0, 3.0, 3.0, 3.0]
    c2 = [1.0, 1.0, 1.0, 1.0, 0.0]
    r2 = [1.0, 2.0, 3.0, 4.0, 5.0]
    u2 = zeros(Float64, n2)

    ASAP.tridag!(a2, b2, c2, r2, u2)

    @test length(u2) == n2
    @test all(isfinite.(u2))
    # 验证 A*u ≈ r
    Av = [b2[1]*u2[1] + c2[1]*u2[2],
          a2[2]*u2[1] + b2[2]*u2[2] + c2[2]*u2[3],
          a2[3]*u2[2] + b2[3]*u2[3] + c2[3]*u2[4],
          a2[4]*u2[3] + b2[4]*u2[4] + c2[4]*u2[5],
          a2[5]*u2[4] + b2[5]*u2[5]]
    @test Av ≈ r2 atol=1e-10
end

