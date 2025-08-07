export soilfluxes, tridag

"""
计算土壤中的水分通量和运动

# 参数
- `i::Int`, `j::Int`: 网格坐标
- `nzg::Int`: 土壤层数
- `freedrain::Int`: 自由排水标志
- `dtll::Float64`: 时间步长 (s)
- `slz::Vector{Float64}`: 土壤层深度 (m) [nzg+1]
- `dz::Vector{Float64}`: 土壤层厚度 (m) [nzg]
- `soiltxt::Int`: 土壤类型
- `smoiwtd::Float64`: 地下水含水量
- `transp::Vector{Float64}`: 蒸腾提取量 (m) [nzg]
- `transpdeep::Float64`: 深层蒸腾 (m)
- `smoi::Vector{Float64}`: 土壤含水量 [nzg]
- `wtd::Float64`: 地下水位深度 (m)
- `precip::Float64`: 降水 (mm)
- `pet_s::Float64`: 土壤蒸发潜力 (mm)
- `fdepth::Float64`: 根系深度因子 (m)
- `qlat::Float64`: 侧向流 (m)
- `qrf::Float64`: 河流反馈流 (m)
- `flood::Float64`: 洪水深度 (m)
- `icefactor::Vector{Int8}`: 冰冻因子 [nzg]
- `smoieq::Vector{Float64}`: 平衡含水量 [nzg]
- `o18::Vector{Float64}`: 氧18同位素 [nzg]
- `precipo18::Float64`: 降水氧18
- `tempsfc::Float64`: 地表温度 (K)
- `qlato18::Float64`: 侧向流氧18
- `transpo18::Float64`: 蒸腾氧18

# 返回
- `Tuple`: (et_s, runoff, rech, flux, qrfcorrect, transpo18, updated_smoi, updated_o18)
"""
function soilfluxes(
    i::Int, j::Int, nzg::Int, freedrain::Int, dtll::Float64,
    slz::Vector{Float64}, dz::Vector{Float64}, soiltxt::Int,
    smoiwtd::Float64, transp::Vector{Float64}, transpdeep::Float64,
    θ::Vector{Float64}, wtd::Float64, precip::Float64, pet_s::Float64,
    fdepth::Float64, qlat::Float64, qrf::Float64, flood::Float64,
    icefactor::Vector{Int8}, θ_eq::Vector{Float64},
    o18::Vector{Float64}, precipo18::Float64, tempsfc::Float64, qlato18::Float64, transpo18::Float64)

    # 计算辅助变量
    vctr2 = 1.0 ./ dz  # 层厚度倒数
    vctr4 = 0.5 .* (slz[1:nzg] .+ slz[2:nzg+1])  # 层中心深度
    vctr5 = zeros(Float64, nzg)
    vctr6 = zeros(Float64, nzg)

    for k in 2:nzg
        vctr5[k] = vctr4[k] - vctr4[k-1]
        vctr6[k] = 1.0 / vctr5[k]
    end

    # 初始化变量
    kfmid = zeros(Float64, nzg)
    diffmid = zeros(Float64, nzg)
    vt3di = zeros(Float64, nzg + 1)
    gravflux = zeros(Float64, nzg + 1)
    capflux = zeros(Float64, nzg + 1)
    flux = zeros(Float64, nzg + 1)
    qlatflux = zeros(Float64, nzg + 2)

    rech = 0.0
    runoff = 0.0

    # 复制初始值
    o18ratio = o18 ./ θ
    θ_old = copy(θ)

    qgw = qlat - qrf

    # 顶部边界条件：入渗 + 土壤蒸发
    soil_params = get_soil_params(soiltxt)
    smoicp = soil_params.ρb

    if θ[nzg] <= smoicp
        pet_s_actual = 0.0
    else
        pet_s_actual = pet_s
    end

    vt3di[nzg+1] = (-precip + pet_s_actual) * 1.0e-3 - flood

    # 检查入渗能力
    if -vt3di[nzg+1] > soil_params.Ksat * dtll
        runoff = -vt3di[nzg+1] - soil_params.Ksat * dtll
        vt3di[nzg+1] = -soil_params.Ksat * dtll
    end

    # 确定地下水位位置
    if freedrain == 0
        iwtd = 1
        for k in 1:nzg
            if wtd < slz[k]
                iwtd = k
                break
            end
        end
    else
        iwtd = 0
    end

    # 计算中间层的导水率和扩散率
    for k in max(iwtd - 1, 2):nzg
        θmid = θ[k] + (θ[k] - θ[k-1]) * (slz[k] - vctr4[k]) * vctr6[k]

        # 获取土壤参数
        soil_params = get_soil_params(soiltxt)
        (; Ksat, θ_sat, ψsat, b) = soil_params

        _Ksat = soil_params.Ksat * max(min(exp((slz[k] + 1.5) / fdepth), 1.0), 0.1)
        _θsat = soil_params.θ_sat * max(min(exp((slz[k] + 1.5) / fdepth), 1.0), 0.1)
        _ψsat = soil_params.ψsat * min(max(exp(-(slz[k] + 1.5) / fdepth), 1.0), 10.0)

        θmid = min(θmid, _θsat)

        # 考虑冰冻因子
        f_ice = icefactor[k] == 0 ? 1.0 : 0.0

        kfmid[k] = f_ice * _Ksat * (θmid / _θsat)^(2.0 * soil_params.b + 3.0)
        diffmid[k] = -f_ice * (_Ksat * _ψsat * soil_params.b / _θsat) *
                     (θmid / _θsat)^(soil_params.b + 2.0)
    end

    # 计算三对角矩阵元素
    aa = zeros(Float64, nzg)
    bb = zeros(Float64, nzg)
    cc = zeros(Float64, nzg)
    rr = zeros(Float64, nzg)

    for k in max(iwtd, 3):nzg
        aa[k] = diffmid[k] * vctr6[k]
        cc[k] = diffmid[k+1] * vctr6[k+1] # TODO: bug here, 数组越界
        bb[k] = -(aa[k] + cc[k] + dz[k] / dtll)
        rr[k] = -θ[k] * dz[k] / dtll - kfmid[k+1] + kfmid[k] + transp[k] / dtll
    end

    # 顶部边界条件
    if iwtd - 1 == nzg
        aa[nzg] = 0.0
        cc[nzg] = 0.0
        bb[nzg] = -dz[nzg] / dtll
        rr[nzg] = vt3di[nzg+1] / dtll - θ[nzg] * dz[nzg] / dtll + transp[nzg] / dtll +
                  min(kfmid[nzg] + diffmid[nzg] * vctr6[nzg] * (θ[nzg] - θ[nzg-1]), 0.0)
    else
        aa[nzg] = diffmid[nzg] * vctr6[nzg]
        cc[nzg] = 0.0
        bb[nzg] = -aa[nzg] - dz[nzg] / dtll
        rr[nzg] = vt3di[nzg+1] / dtll - θ[nzg] * dz[nzg] / dtll + kfmid[nzg] + transp[nzg] / dtll
    end

    # 底部边界条件
    if freedrain != 1
        if iwtd <= 2
            aa[1] = 0.0
            cc[1] = diffmid[2] * vctr6[2]
            bb[1] = -(cc[1] + dz[1] / dtll)
            rr[1] = -θ[1] * dz[1] / dtll - kfmid[2] + transp[1] / dtll

            k = 2
            aa[k] = diffmid[k] * vctr6[k]
            cc[k] = diffmid[k+1] * vctr6[k+1]
            bb[k] = -(aa[k] + cc[k] + dz[k] / dtll)
            rr[k] = -θ[k] * dz[k] / dtll - kfmid[k+1] + kfmid[k] + transp[k] / dtll
        else
            # 其他层保持不变
            for k in 1:(iwtd-3)
                aa[k] = 0.0
                cc[k] = 0.0
                bb[k] = 1.0
                rr[k] = θ[k]
            end

            k = iwtd - 1  # 地下水位层
            aa[k] = 0.0
            cc[k] = diffmid[k+1] * vctr6[k+1]
            bb[k] = -(cc[k] + dz[k] / dtll)
            rr[k] = -θ[k] * dz[k] / dtll - kfmid[k+1] + transp[k] / dtll +
                    min(kfmid[k] + diffmid[k] * vctr6[k] * (θ[k] - θ[k-1]), 0.0)

            k = iwtd - 2
            aa[k] = 0.0
            cc[k] = 0.0
            bb[k] = -dz[k] / dtll
            rr[k] = -θ[k] * dz[k] / dtll +
                    max(-kfmid[k+1] - diffmid[k+1] * vctr6[k+1] * (θ[k+1] - θ[k]), 0.0)
        end
    else
        # 重力排水边界
        soil_params = get_soil_params(soiltxt)
        hydcon = soil_params.Ksat * max(min(exp((slz[1] + 1.5) / fdepth), 1.0), 0.1)
        smoisat = soil_params.θ_sat * max(min(exp((slz[1] + 1.5) / fdepth), 1.0), 0.1)

        kfmid[1] = hydcon * (θ[1] / smoisat)^(2.0 * soil_params.b + 3.0)

        aa[1] = 0.0
        cc[1] = diffmid[2] * vctr6[2]
        bb[1] = -(cc[1] + dz[1] / dtll)
        rr[1] = -θ[1] * dz[1] / dtll - kfmid[2] + kfmid[1] + transp[1] / dtll
    end

    # 求解三对角系统
    tridag!(aa, bb, cc, rr, θ)

    # 计算通量
    for k in max(iwtd, 3):nzg
        gravflux[k] = -kfmid[k] * dtll
        capflux[k] = -aa[k] * (θ[k] - θ[k-1]) * dtll
        vt3di[k] = capflux[k] + gravflux[k]
    end

    if iwtd <= 2
        capflux[1] = 0.0
        gravflux[1] = 0.0
        vt3di[1] = 0.0
        gravflux[2] = -kfmid[2] * dtll
        capflux[2] = -aa[2] * (θ[2] - θ[1]) * dtll
        vt3di[2] = capflux[2] + gravflux[2]
    else
        for k in 1:(iwtd-2)
            capflux[k] = 0.0
            gravflux[k] = 0.0
            vt3di[k] = 0.0
        end

        k = iwtd - 1
        gravflux[k] = -kfmid[k] * dtll
        capflux[k] = -diffmid[k] * vctr6[k] * (θ_old[k] - θ_old[k-1]) * dtll

        if capflux[k] > -gravflux[k]
            vt3di[k] = capflux[k] + gravflux[k]
        elseif capflux[k] > 0.0
            gravflux[k] = -capflux[k]
            vt3di[k] = 0.0
        else
            capflux[k] = 0.0
            gravflux[k] = 0.0
            vt3di[k] = 0.0
        end
    end

    if freedrain == 0
        vt3di[1] = 0.0
    else
        vt3di[1] = -kfmid[1] * dtll
    end

    # 重新计算土壤含水量
    θ .= θ_old
    for k in 1:nzg
        θ_old[k] = θ_old[k] + (vt3di[k] - vt3di[k+1] - transp[k]) * vctr2[k]
    end

    # 检查并修正土壤含水量边界
    for k in 1:nzg
        soil_params = get_soil_params(soiltxt)
        smoisat = soil_params.θ_sat * max(min(exp((vctr4[k] + 1.5) / fdepth), 1.0), 0.1)

        if θ_old[k] > smoisat
            dsmoi = max((θ_old[k] - smoisat) * dz[k], 0.0)
            if k < nzg
                θ_old[k+1] = θ_old[k+1] + dsmoi * vctr2[k+1]
                vt3di[k+1] = vt3di[k+1] + dsmoi
            else
                vt3di[k+1] = vt3di[k+1] + dsmoi
                runoff = runoff + dsmoi
            end
            θ_old[k] = smoisat

            if capflux[k+1] < 0.0
                gravflux[k+1] = gravflux[k+1] + capflux[k+1]
                capflux[k+1] = 0.0
            end
            gravflux[k+1] = gravflux[k+1] + dsmoi
            if gravflux[k+1] > 0.0
                capflux[k+1] = capflux[k+1] + gravflux[k+1]
                gravflux[k+1] = 0.0
            end
        end
    end

    # 处理最上层的干燥限制
    k = nzg
    soil_params = get_soil_params(soiltxt)
    smoicp = soil_params.ρb * max(min(exp((vctr4[k] + 1.5) / fdepth), 1.0), 0.1)

    et_s = pet_s
    if θ_old[k] < smoicp
        dsmoi = max((smoicp - θ_old[k]) * dz[k], 0.0)
        if vt3di[k+1] > dsmoi
            et_s = max(0.0, pet_s - dsmoi * 1.0e3)
            θ_old[k] = smoicp
            vt3di[k+1] = vt3di[k+1] - dsmoi
        else
            et_s = max(0.0, pet_s - max(vt3di[k+1], 0.0) * 1.0e3)
            vt3di[k+1] = min(vt3di[k+1], 0.0)
            θ_old[k] = θ_old[k] + max(vt3di[k+1], 0.0) / dz[k]

            dsmoi = max((smoicp - θ_old[k]) * dz[k], 0.0)
            θ_old[k-1] = θ_old[k-1] - dsmoi * vctr2[k-1]
            vt3di[k] = vt3di[k] + dsmoi
            θ_old[k] = smoicp
        end
    end

    # 继续处理下层
    for k in (nzg-1):-1:1
        soil_params = get_soil_params(soiltxt)
        smoicp = soil_params.ρb * max(min(exp((vctr4[k] + 1.5) / fdepth), 1.0), 0.1)

        if θ_old[k] < smoicp
            dsmoi = max((smoicp - θ_old[k]) * dz[k], 0.0)
            if k > 1
                θ_old[k-1] = θ_old[k-1] - dsmoi * vctr2[k-1]
            end
            vt3di[k] = vt3di[k] + dsmoi
            θ_old[k] = smoicp
        end
    end

    qrfcorrect = 0.0
    if vt3di[1] > 0.0
        qrfcorrect = -min(vt3di[1], max(qrf, 0.0))
    end

    # 累积通量
    flux .+= vt3di

    # 处理自由排水
    if freedrain == 1
        rech = vt3di[1]
        # smoiwtd = smoiwtd - vt3di[1]  # 这个变量需要在调用函数中处理
    end

    # 氧18同位素计算（简化版本）
    transpo18_out = 0.0
    for k in 1:nzg
        transpo18_out += 0.5 * (o18[k] / θ_old[k] + o18ratio[k]) * transp[k]
        if o18[k] < 0.0
            println("警告: O18小于零! i=$i, j=$j, k=$k")
        end
        if o18[k] > θ_old[k]
            println("警告: O18大于土壤含水量! i=$i, j=$j, k=$k")
        end
    end

    o18 .= max.(o18, 0.0)
    θ .= θ_old
    return et_s, runoff, rech, flux, qrfcorrect, transpo18_out, θ, o18
end


"""
三对角矩阵求解器（原地修改版本）

# 参数
- `a::Vector{Float64}`: 下对角线元素
- `b::Vector{Float64}`: 主对角线元素  
- `c::Vector{Float64}`: 上对角线元素
- `r::Vector{Float64}`: 右端向量
- `u::Vector{Float64}`: 解向量（输出）
"""
function tridag!(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64},
    r::Vector{Float64}, u::Vector{Float64})
    n = length(b)
    gam = zeros(Float64, n)

    if b[1] == 0.0
        error("tridag: 需要重写方程")
    end

    bet = b[1]
    u[1] = r[1] / bet

    for j in 2:n
        gam[j] = c[j-1] / bet
        bet = b[j] - a[j] * gam[j]
        if bet == 0.0
            error("tridag 失败")
        end
        u[j] = (r[j] - a[j] * u[j-1]) / bet
    end

    for j in (n-1):-1:1
        u[j] = u[j] - gam[j+1] * u[j+1]
    end
end
