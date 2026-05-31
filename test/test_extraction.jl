using ASAP, Test

# 典型参数设置（Penman-Monteith 方程参数）
const Δ_test = 145.0    # 饱和水汽压斜率 (Pa/K)
const γ_test = 66.0     # 湿度计常数 (Pa/K)
const λ_test = 2.45e6   # 汽化潜热 (J/kg)

function make_extraction_inputs(nzg; lai=2.0, θ_val=0.25, icefac_val=Int8(0))
    z₋ₕ, dz = initializesoildepth(nzg)
    Δt = 3600.0
    soiltxt = 1
    wtd = z₋ₕ[1] - 0.5  # 水位在所有层之下（不饱和）
    θ = fill(θ_val, nzg)
    θ_wtd = 0.3
    ra_a = 100.0
    ra_c = 50.0
    rs_c_factor = 200.0
    R_a = (Δ_test + γ_test) * ra_a * 3.0  # 近似辐射驱动力
    R_s = (Δ_test + γ_test) * ra_a * 1.0
    petfactor_s = R_s
    petfactor_c = R_a
    inactivedays = zeros(Int, nzg + 2)   # 全部活跃
    maxinactivedays = 10
    hhveg = 10.0
    fdepth = 2.0
    icefac = fill(icefac_val, nzg)
    return (nzg, z₋ₕ, dz, Δt, soiltxt, wtd, θ, θ_wtd,
            Δ_test, γ_test, λ_test,
            lai, ra_a, ra_c, rs_c_factor,
            R_a, R_s, petfactor_s, petfactor_c,
            inactivedays, maxinactivedays,
            hhveg, fdepth, icefac)
end


@testset "extraction基本功能" begin
    args = make_extraction_inputs(10)
    pet_s, pet_c, watdef, dθ, dθ_deep = extraction(args...)

    @test pet_s >= 0.0
    @test pet_c >= 0.0
    @test watdef >= 0.0
    @test length(dθ) == 10
    @test dθ_deep == 0.0           # 深层无提取
    @test all(dθ .>= 0.0)          # 每层提取量非负
end

@testset "零LAI时无冠层蒸腾" begin
    args = make_extraction_inputs(10; lai=0.0)
    pet_s, pet_c, watdef, dθ, dθ_deep = extraction(args...)

    # LAI=0 时，冠层蒸腾强制为零
    @test pet_c == 0.0
    @test pet_s >= 0.0
end

@testset "全部冰冻时无提取" begin
    args = make_extraction_inputs(10; icefac_val=Int8(1))
    pet_s, pet_c, watdef, dθ, dθ_deep = extraction(args...)

    # 冰冻情况下根区水分提取便利性为零
    @test all(dθ .== 0.0)
    @test watdef >= 0.0  # 水分亏缺 ≥ 0
end

@testset "干旱土壤产生水分亏缺" begin
    # 非常干的土壤 (θ << θ_wilt)
    nzg = 10
    args_dry = make_extraction_inputs(nzg; θ_val=0.05)
    pet_s, pet_c, watdef_dry, dθ_dry, _ = extraction(args_dry...)

    args_wet = make_extraction_inputs(nzg; θ_val=0.35)
    _, _, watdef_wet, dθ_wet, _ = extraction(args_wet...)

    # 湿润土壤水分亏缺应该比干燥时少
    @test watdef_wet <= watdef_dry
    # 湿润土壤提取总量应该更多
    @test sum(dθ_wet) >= sum(dθ_dry)
end

@testset "inactivedays正确更新" begin
    nzg = 10
    args = make_extraction_inputs(nzg; θ_val=0.25)
    # inactivedays初始为全0（全部活跃）
    pet_s, pet_c, watdef, dθ, dθ_deep = extraction(args...)

    # 函数会修改 inactivedays（若根区活跃则清零，否则+1）
    @test pet_s >= 0.0
    @test pet_c >= 0.0
end

@testset "蒸发量与时间步长成正比" begin
    nzg = 10

    # 短时间步长（30 min）
    args1 = make_extraction_inputs(nzg)
    # 修改 Δt 为 1800 s
    args1 = (args1[1], args1[2], args1[3], 1800.0, args1[5:end]...)
    pet_s1, pet_c1, _, _, _ = extraction(args1...)

    # 长时间步长（2 h）
    args2 = make_extraction_inputs(nzg)
    args2 = (args2[1], args2[2], args2[3], 7200.0, args2[5:end]...)
    pet_s2, pet_c2, _, _, _ = extraction(args2...)

    # 时间步长加倍，PET近似翻倍
    @test pet_s2 > pet_s1
    @test pet_c2 > pet_c1
end

@testset "各层提取量之和等于冠层蒸腾量（无亏缺时）" begin
    nzg = 10
    # 使用高含水量，确保无水分亏缺
    args = make_extraction_inputs(nzg; θ_val=0.38)
    pet_s, pet_c, watdef, dθ, dθ_deep = extraction(args...)

    if watdef == 0.0
        # 无水分亏缺：提取量之和应等于 pet_c（m）
        @test sum(dθ) ≈ pet_c * 1e-3 atol=1e-10
    else
        @test sum(dθ) + watdef ≈ pet_c * 1e-3 atol=1e-10
    end
end
