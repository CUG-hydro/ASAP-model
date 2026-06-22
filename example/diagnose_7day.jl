#!/usr/bin/env julia
# 诊断 7 天模拟: 收集水均衡与状态变量时序, 验证模型行为
#
# 用法: julia --project example/diagnose_7day.jl [--days N]
#
# 与 regional_example.jl::main 共享数据生成与初始化, 但每 24h
# 输出水均衡与状态变量空间分布, 用于诊断模型是否物理合理。
using ASAP, NCDatasets, Statistics, Dates, Printf

# ---- CLI ----
days = 7
for (i, a) in enumerate(ARGS)
    if a == "--days" && i < length(ARGS)
        days = parse(Int, ARGS[i+1])
    end
end
println("ASAP 7-day diagnose run: days=$days")

# ---- 复用区域示例的辅助 ----
# 这些函数在 regional_example.jl 中定义, 通过 include 复用
include(joinpath(@__DIR__, "regional_example.jl"))

# ---- 生成 mock 数据 ----
nx, ny, nzg = 12, 10, 40
paths = generate_mock_dataset("20200701", nx, ny, nzg)

# ---- 初始化 (与 main 一致) ----
is, ie = 1, nx; js, je = 1, ny
initial = read_initial(paths.static)
soiltxt_raw, topo, fdepth, landmask = initial.soiltxt, initial.topo,
                                          initial.fdepth, initial.landmask
wtd_in = read_wtdnc(paths.wtd)
z₋ₕ, dz = initializesoildepth(nzg)

smoieq = zeros(Float64, nzg, nx, ny)
soiltxt = ones(Int, 2, nx, ny)
for j in 1:ny, i in 1:nx
    nsoil = max(1, min(13, Int(soiltxt_raw[i, j])))
    soiltxt[1, i, j] = nsoil
    soiltxt[2, i, j] = nsoil
    smoieq[:, i, j] .= eqsoilmoisturetheor(nzg, nsoil, z₋ₕ, dz, fdepth[i, j], wtd_in[i, j])
end

freedrain, maxinactivedays, Δt, steps = 1, 240, 3600.0, 1.0
hveg = fill(1.0, nx, ny)
veg  = fill(7.0, nx, ny)

# 启动用平衡态 smoi → 无 dS 初始激扰
smoi = copy(smoieq)
sm0_top = copy(smoi[nzg, :, :])   # 表层初始

smoiwtd      = fill(0.35, nx, ny)
wtd          = copy(wtd_in)
waterdeficit = zeros(Float64, nx, ny)
watext       = zeros(Float64, nzg, nx, ny)
watextdeep   = zeros(Float64, nx, ny)
rech         = zeros(Float64, nx, ny)
deeprech     = zeros(Float64, nx, ny)
et_s = zeros(Float64, nx, ny)
et_i = zeros(Float64, nx, ny)
et_c = zeros(Float64, nx, ny)
intercepstore = zeros(Float64, nx, ny)
ppacum         = zeros(Float64, nx, ny)
pppendepth     = zeros(Float64, nx, ny)
pppendepthold  = zeros(Int8, nx, ny)
qlat         = zeros(Float64, nx, ny)
qlatsum      = zeros(Float64, nx, ny)
qsprings     = zeros(Float64, nx, ny)
inactivedays = zeros(Int, nzg+2, nx, ny)
floodheight  = zeros(Float64, nx, ny)
qrf          = zeros(Float64, nx, ny)
delsfcwat    = zeros(Float64, nx, ny)
icefactor    = zeros(Int8, nx, ny, nzg)
wtdflux      = zeros(Float64, nx, ny)
et_s_daily   = zeros(Float64, nx, ny)
et_c_daily   = zeros(Float64, nx, ny)
transptop    = zeros(Float64, nx, ny)
infilk       = zeros(Int8, nx, ny)
infilflux    = zeros(Float64, nzg, nx, ny)
infilfluxday = zeros(Float64, nzg, nx, ny)
infilcounter = zeros(Int16, nzg, nx, ny)
upflux       = zeros(Float64, nzg, nx, ny)
o18          = zeros(Float64, nzg, nx, ny)
o18ratiopp   = zeros(Float64, nx, ny)
tempsfc      = fill(293.0, nx, ny)
qlato18      = zeros(Float64, nx, ny)
transpo18    = zeros(Float64, nx, ny)

# 水均衡累积
ΣP = 0.0; ΣET_s = 0.0; ΣET_c = 0.0; ΣET_i = 0.0
ΣRch = 0.0

println("\n=== 初始状态 ===")
println("  soiltxt 范围: ", extrema(soiltxt_raw))
println("  topo 范围 (m): ", round.(extrema(topo); digits=1))
println("  wtd 范围 (m): ", round.(extrema(wtd); digits=2))
println("  smoi[顶层] 范围: ", round.(extrema(smoi[nzg, :, :]); digits=3))
println("  smoi[底层] 范围: ", round.(extrema(smoi[1, :, :]); digits=3))

function run_simulation(days::Int)
    local ΣP = 0.0
    local ΣET_s = 0.0
    local ΣET_c = 0.0
    local ΣET_i = 0.0
    local ΣRch = 0.0

    for hour in 1:(days * 24)
        hour_in_day = mod1(hour, 24)
        day_offset  = div(hour - 1, 24)
        cur_date    = "202007" * lpad(string(1 + day_offset), 2, '0')

        paths_h = (day_offset > 0) ? paths : paths  # mock 模式: 复用单日路径
        f = read_hourly_forcings(hour_in_day, paths_h)

        icefactor .= 0

        et_s0 = copy(et_s); et_c0 = copy(et_c); et_i0 = copy(et_i)
        rech0 = copy(rech);   wtdflux0 = copy(wtdflux)

        rootdepth_main(freedrain, is, ie, js, je, nzg, z₋ₕ, dz, Δt,
            landmask, veg, hveg,
            soiltxt, f.wind, f.temp, f.qair, f.press, f.netrad,
            f.rshort, f.lai, f.precip,
            fill(0.0, nx, ny), smoi, smoieq,
            smoiwtd, wtd, waterdeficit,
            watext, watextdeep, rech, deeprech,
            et_s, et_i, et_c,
            intercepstore, ppacum,
            pppendepth, pppendepthold, qlat,
            qlatsum, qsprings, inactivedays,
            maxinactivedays, fdepth,
            steps, floodheight, qrf,
            delsfcwat, icefactor, wtdflux,
            et_s_daily, et_c_daily, transptop,
            infilk, infilflux, infilfluxday,
            infilcounter, hour,
            o18, o18ratiopp, tempsfc,
            qlato18, transpo18, upflux)

        # 增量(空间平均, mm/h)
        dP  = mean(f.precip)
        dET = (mean(et_s) - mean(et_s0)) + (mean(et_c) - mean(et_c0)) + (mean(et_i) - mean(et_i0))
        dR  = mean(rech) - mean(rech0)
        ΣP += dP; ΣET_s += mean(et_s) - mean(et_s0); ΣET_c += mean(et_c) - mean(et_c0)
        ΣET_i += mean(et_i) - mean(et_i0); ΣRch = mean(rech)

        # 每 24h 报一次
        if hour_in_day == 24
            dS_top = sum((smoi[nzg, :, :] .- sm0_top)) / (nx * ny) * dz[nzg] * 1000
            println("\n--- Day $(day_offset + 1) end ---")
            @printf("  日均 P       : %6.3f mm/h  → %6.1f mm/day\n", dP, dP * 24)
            @printf("  日均 ET_s    : %6.3f mm/h  → %6.1f mm/day\n",
                    mean(et_s) - mean(et_s0), (mean(et_s) - mean(et_s0)) * 24)
            @printf("  日均 ET_c    : %6.3f mm/h  → %6.1f mm/day\n",
                    mean(et_c) - mean(et_c0), (mean(et_c) - mean(et_c0)) * 24)
            @printf("  日均 ET_i    : %6.3f mm/h  → %6.1f mm/day\n",
                    mean(et_i) - mean(et_i0), (mean(et_i) - mean(et_i0)) * 24)
            @printf("  水均衡检查 P - ET = ΔS + R:  %6.2f - %6.2f = %6.2f ; ΔS_top + R = %6.2f + %6.2f = %6.2f\n",
                    dP * 24, dET * 24, dP*24 - dET*24, dS_top, dR, dS_top + dR)
            println("  smoi[顶层] 范围: ", round.(extrema(smoi[nzg, :, :]); digits=3))
            println("  smoi[中层 nzg÷2] 范围: ", round.(extrema(smoi[nzg ÷ 2, :, :]); digits=3))
            println("  smoi[底层] 范围: ", round.(extrema(smoi[1, :, :]); digits=3))
            println("  wtd 范围 (m): ", round.(extrema(wtd); digits=2))
            println("  waterdeficit 范围: ", round.(extrema(waterdeficit); digits=1), " mm")
            println("  ----- 累计: P=", round(ΣP*24; digits=1), " ET=", round((ΣET_s+ΣET_c+ΣET_i)*24; digits=1))
        end
    end
    return ΣP, ΣET_s + ΣET_c + ΣET_i, ΣRch
end

ΣP_total, ΣET_total, ΣR_total = run_simulation(days)

println("\n=== 模拟完成 (共 $days 天) ===")
println("最终 smoi[顶层] 范围: ", round.(extrema(smoi[nzg, :, :]); digits=3))
println("最终 smoi[底层] 范围: ", round.(extrema(smoi[1, :, :]); digits=3))
println("最终 wtd 范围 (m): ", round.(extrema(wtd); digits=2))
println("最终 waterdeficit 范围: ", round.(extrema(waterdeficit); digits=1), " mm")
println("最终累计水均衡: P=$(round(ΣP*24; digits=1)) mm, ET=$(round((ΣET_s+ΣET_c+ΣET_i)*24; digits=1)) mm, R=$(round(ΣRch; digits=1)) mm")
