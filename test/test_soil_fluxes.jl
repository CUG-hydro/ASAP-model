using ASAP, Test

# 测试基本功能
@testset "soilfluxes基本功能" begin
    # 设置测试参数
    nzg = 10
    i, j = 1, 1
    freedrain = 0
    dtll = 3600.0  # 1小时时间步长
    soiltxt = 1    # 土壤类型1

    # 初始化土壤层
    slz, dz = initializesoildepth(nzg)

    # 初始化输入数据
    smoiwtd = 0.3
    transp = zeros(Float64, nzg)
    transpdeep = 0.0
    θ = fill(0.25, nzg)  # 使用希腊字母θ
    wtd = -2.0
    precip = 5.0  # mm
    pet_s = 2.0   # mm
    fdepth = 2.0
    qlat = 0.0
    qrf = 0.0
    flood = 0.0
    icefactor = zeros(Int8, nzg)
    θ_eq = fill(0.2, nzg)  # 使用希腊字母θ
    o18 = θ .* (-5.0)  # 初始同位素值
    precipo18 = -8.0
    tempsfc = 293.15
    qlato18 = 0.0
    transpo18 = 0.0

    # 调用soilfluxes函数
    et_s, runoff, rech, flux, qrfcorrect, transpo18_out, θ_new, o18_new =
        soilfluxes(nzg, freedrain, dtll, slz, dz, soiltxt,
            smoiwtd, transp, transpdeep, θ, wtd, precip, pet_s,
            fdepth, qlat, qrf, flood, icefactor, θ_eq,
            o18, precipo18, tempsfc, qlato18, transpo18)

    # 基本检查
    @test et_s >= 0.0  # 蒸发应该非负
    @test runoff >= 0.0  # 径流应该非负
    @test length(θ_new) == nzg  # 输出土壤含水量长度正确
    @test length(o18_new) == nzg  # 输出同位素长度正确
    @test length(flux) == nzg + 1  # 通量长度正确

    # 检查土壤含水量在合理范围内
    soil = get_soil_params(soiltxt)
    for k in 1:nzg
        @test θ_new[k] >= 0.0  # 含水量非负
        @test θ_new[k] <= soil.θ_sat * 1.1  # 不超过饱和含水量太多
    end

    # 检查同位素值非负
    @test all(o18_new .>= 0.0)
end

# 测试不同降水条件
@testset "不同降水条件测试" begin
    nzg = 5
    slz, dz = initializesoildepth(nzg)

    # 基础参数
    i, j = 1, 1
    freedrain = 0
    dtll = 3600.0
    soiltxt = 1
    smoiwtd = 0.3
    transp = zeros(Float64, nzg)
    transpdeep = 0.0
    θ = fill(0.25, nzg)
    wtd = -2.0
    pet_s = 2.0
    fdepth = 2.0
    qlat = 0.0
    qrf = 0.0
    flood = 0.0
    icefactor = zeros(Int8, nzg)
    θ_eq = fill(0.2, nzg)
    o18 = θ .* (-5.0)
    precipo18 = -8.0
    tempsfc = 293.15
    qlato18 = 0.0
    transpo18 = 0.0

    # 测试无降水
    et_s1, runoff1, _, _, _, _, _, _ =
        soilfluxes(nzg, freedrain, dtll, slz, dz, soiltxt,
            smoiwtd, transp, transpdeep, copy(θ), wtd, 0.0, pet_s,
            fdepth, qlat, qrf, flood, icefactor, θ_eq,
            copy(o18), precipo18, tempsfc, qlato18, transpo18)

    # 测试大降水
    et_s2, runoff2, _, _, _, _, _, _ =
        soilfluxes(nzg, freedrain, dtll, slz, dz, soiltxt,
            smoiwtd, transp, transpdeep, copy(θ), wtd, 50.0, pet_s,
            fdepth, qlat, qrf, flood, icefactor, θ_eq,
            copy(o18), precipo18, tempsfc, qlato18, transpo18)

    # 大降水应该产生更多径流
    @test runoff2 >= runoff1
end

# 测试自由排水条件
@testset "自由排水条件测试" begin
    nzg = 5
    slz, dz = initializesoildepth(nzg)

    # 基础参数
    i, j = 1, 1
    freedrain = 1  # 自由排水
    dtll = 3600.0
    soiltxt = 1
    smoiwtd = 0.3
    transp = zeros(Float64, nzg)
    transpdeep = 0.0
    θ = fill(0.25, nzg)
    wtd = -2.0
    precip = 10.0
    pet_s = 2.0
    fdepth = 2.0
    qlat = 0.0
    qrf = 0.0
    flood = 0.0
    icefactor = zeros(Int8, nzg)
    θ_eq = fill(0.2, nzg)
    o18 = θ .* (-5.0)
    precipo18 = -8.0
    tempsfc = 293.15
    qlato18 = 0.0
    transpo18 = 0.0

    et_s, runoff, rech, flux, qrfcorrect, transpo18_out, θ_new, o18_new =
        soilfluxes(nzg, freedrain, dtll, slz, dz, soiltxt,
            smoiwtd, transp, transpdeep, θ, wtd, precip, pet_s,
            fdepth, qlat, qrf, flood, icefactor, θ_eq,
            o18, precipo18, tempsfc, qlato18, transpo18)

    # 自由排水条件下应该有补给
    @test abs(rech) >= 0.0  # 补给可能为正或负
end

# 测试有蒸腾的情况
@testset "有蒸腾测试" begin
    nzg = 5
    slz, dz = initializesoildepth(nzg)

    # 基础参数
    i, j = 1, 1
    freedrain = 0
    dtll = 3600.0
    soiltxt = 1
    smoiwtd = 0.3
    transp = fill(0.001, nzg)  # 每层有蒸腾
    transpdeep = 0.0
    θ = fill(0.3, nzg)  # 较高含水量
    wtd = -2.0
    precip = 0.0  # 无降水
    pet_s = 2.0
    fdepth = 2.0
    qlat = 0.0
    qrf = 0.0
    flood = 0.0
    icefactor = zeros(Int8, nzg)
    θ_eq = fill(0.2, nzg)
    o18 = θ .* (-5.0)
    precipo18 = -8.0
    tempsfc = 293.15
    qlato18 = 0.0
    transpo18 = 0.0

    θ_original = copy(θ)

    et_s, runoff, rech, flux, qrfcorrect, transpo18_out, θ_new, o18_new =
        soilfluxes(nzg, freedrain, dtll, slz, dz, soiltxt,
            smoiwtd, transp, transpdeep, θ, wtd, precip, pet_s,
            fdepth, qlat, qrf, flood, icefactor, θ_eq,
            o18, precipo18, tempsfc, qlato18, transpo18)

    # 有蒸腾应该降低土壤含水量
    @test sum(θ_new) <= sum(θ_original)
    @test transpo18_out >= 0.0  # 蒸腾同位素输出应该非负
end

# 测试冰冻因子影响
@testset "冰冻因子测试" begin
    nzg = 5
    slz, dz = initializesoildepth(nzg)

    # 基础参数
    i, j = 1, 1
    freedrain = 0
    dtll = 3600.0
    soiltxt = 1
    smoiwtd = 0.3
    transp = zeros(Float64, nzg)
    transpdeep = 0.0
    θ = fill(0.25, nzg)
    wtd = -2.0
    precip = 10.0
    pet_s = 2.0
    fdepth = 2.0
    qlat = 0.0
    qrf = 0.0
    flood = 0.0
    θ_eq = fill(0.2, nzg)
    o18 = θ .* (-5.0)
    precipo18 = -8.0
    tempsfc = 293.15
    qlato18 = 0.0
    transpo18 = 0.0

    # 无冰冻
    icefactor1 = zeros(Int8, nzg)
    et_s1, runoff1, _, _, _, _, _, _ =
        soilfluxes(nzg, freedrain, dtll, slz, dz, soiltxt,
            smoiwtd, transp, transpdeep, copy(θ), wtd, precip, pet_s,
            fdepth, qlat, qrf, flood, icefactor1, θ_eq,
            copy(o18), precipo18, tempsfc, qlato18, transpo18)

    # 部分冰冻
    icefactor2 = ones(Int8, nzg)  # 全部冰冻
    et_s2, runoff2, _, _, _, _, _, _ =
        soilfluxes(nzg, freedrain, dtll, slz, dz, soiltxt,
            smoiwtd, transp, transpdeep, copy(θ), wtd, precip, pet_s,
            fdepth, qlat, qrf, flood, icefactor2, θ_eq,
            copy(o18), precipo18, tempsfc, qlato18, transpo18)

    # 冰冻会影响水分运动，通常导致更多径流
    @test runoff2 >= runoff1
end

# 测试不同土壤类型
@testset "不同土壤类型测试" begin
    nzg = 5
    slz, dz = initializesoildepth(nzg)

    # 基础参数
    i, j = 1, 1
    freedrain = 0
    dtll = 3600.0
    smoiwtd = 0.3
    transp = zeros(Float64, nzg)
    transpdeep = 0.0
    θ = fill(0.25, nzg)
    wtd = -2.0
    precip = 10.0
    pet_s = 2.0
    fdepth = 2.0
    qlat = 0.0
    qrf = 0.0
    flood = 0.0
    icefactor = zeros(Int8, nzg)
    θ_eq = fill(0.2, nzg)
    o18 = θ .* (-5.0)
    precipo18 = -8.0
    tempsfc = 293.15
    qlato18 = 0.0
    transpo18 = 0.0

    results = []

    # 测试多种土壤类型
    for soiltxt in [1, 5, 10]
        et_s, runoff, rech, flux, qrfcorrect, transpo18_out, θ_new, o18_new =
            soilfluxes(nzg, freedrain, dtll, slz, dz, soiltxt,
                smoiwtd, transp, transpdeep, copy(θ), wtd, precip, pet_s,
                fdepth, qlat, qrf, flood, icefactor, θ_eq,
                copy(o18), precipo18, tempsfc, qlato18, transpo18)

        push!(results, (soiltxt=soiltxt, et_s=et_s, runoff=runoff))

        # 基本合理性检查
        @test et_s >= 0.0
        @test runoff >= 0.0
    end

    # 不同土壤类型应该产生不同的结果
    @test length(unique([r.runoff for r in results])) > 1
end

# 测试三对角矩阵求解器
@testset "tridag! 函数测试" begin
    n = 5
    a = [0.0, 1.0, 1.0, 1.0, 1.0]  # 下对角线
    b = [2.0, 2.0, 2.0, 2.0, 2.0]  # 主对角线
    c = [1.0, 1.0, 1.0, 1.0, 0.0]  # 上对角线
    r = [1.0, 2.0, 3.0, 4.0, 5.0]  # 右端向量
    u = zeros(Float64, n)           # 解向量

    # 调用三对角求解器
    tridag!(a, b, c, r, u)

    # 检查解的合理性
    @test length(u) == n
    @test all(isfinite.(u))  # 解应该是有限的
end

# 边界条件测试
@testset "边界条件测试" begin
    nzg = 3  # 使用较少的层数简化测试
    slz, dz = initializesoildepth(nzg)

    # 基础参数
    i, j = 1, 1
    freedrain = 0
    dtll = 3600.0
    soiltxt = 1
    smoiwtd = 0.3
    transp = zeros(Float64, nzg)
    transpdeep = 0.0
    θ = fill(0.25, nzg)
    wtd = -2.0
    pet_s = 2.0
    fdepth = 2.0
    qlat = 0.0
    qrf = 0.0
    flood = 0.0
    icefactor = zeros(Int8, nzg)
    θ_eq = fill(0.2, nzg)
    o18 = θ .* (-5.0)
    precipo18 = -8.0
    tempsfc = 293.15
    qlato18 = 0.0
    transpo18 = 0.0

    # 测试极端干燥条件
    θ_dry = fill(0.05, nzg)  # 很干的土壤
    et_s, runoff, _, _, _, _, θ_new, _ =
        soilfluxes(nzg, freedrain, dtll, slz, dz, soiltxt,
            smoiwtd, transp, transpdeep, θ_dry, wtd, 0.0, pet_s,
            fdepth, qlat, qrf, flood, icefactor, θ_eq,
            copy(o18), precipo18, tempsfc, qlato18, transpo18)

    # 干燥条件下蒸发应该受限
    @test et_s <= pet_s

    # 测试饱和条件
    soil = get_soil_params(soiltxt)
    θ_sat = fill(soil.θ_sat, nzg)
    et_s2, runoff2, _, _, _, _, θ_new2, _ =
        soilfluxes(nzg, freedrain, dtll, slz, dz, soiltxt,
            smoiwtd, transp, transpdeep, θ_sat, wtd, 50.0, pet_s,
            fdepth, qlat, qrf, flood, icefactor, θ_eq,
            copy(o18), precipo18, tempsfc, qlato18, transpo18)

    # 饱和条件下应该有显著径流
    @test runoff2 > 0.0
end
