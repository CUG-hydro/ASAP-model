# 水位计算的核心模块
# module WaterTableCalculations

# using ..Constants: π₄
# using ..SoilParameters: SoilType, get_soil_params
export wtable!, updatewtd!, lateral_flow!


"""
    wtable!(水位参数...)

主要的水位计算子程序，处理侧向流、深层补给和水位更新

参数：
- imax, jmax: 网格尺寸
- is, ie, js, je: 计算范围
- nzg: 土壤层数
- slz: 土壤层深度
- dz: 土壤层厚度
- area: 网格面积
- soiltxt: 土壤类型
- wtd: 水位深度
- bottomflux: 底部通量
- rech: 补给量
- qslat: 侧向流
- fdepth: 排水深度
- topo: 地形高程
- landmask: 陆地掩膜
- Δt: 时间步长
- smoi: 土壤湿度
- smoieq: 平衡土壤湿度
- smoiwtd: 水位处土壤湿度
- qsprings: 泉水流量
"""
function wtable!(
    imax::Int, jmax::Int, is::Int, ie::Int, js::Int, je::Int, nzg::Int,
    slz::V, dz::V, area::M,
    soiltxt::Array{Int,3}, wtd::M, bottomflux::M,
    rech::M, qslat::M, fdepth::M,
    topo::M, landmask::Matrix{Int}, Δt::T,
    smoi::Array{T,3}, smoieq::Array{T,3}, smoiwtd::M,
    qsprings::M
) where {T<:AbstractFloat,V<:Vector{T},M<:Matrix{T}}
    # 计算侧向流导水率
    κlat = zeros(size(wtd))

    for j in js:je
        for i in 1:imax
            nsoil = soiltxt[1, i, j]
            soil_params = get_soil_params(nsoil)
            κlat[i, j] = soil_params.K_sat * soil_params.K_latfactor
        end
    end

    # 计算侧向流
    qlat = zeros(size(wtd))
    lateral_flow!(imax, jmax, is, ie, js, je, wtd, qlat, fdepth, topo,
        landmask, Δt, area, κlat)

    # 累积侧向流 (转换为mm)
    qslat .+= qlat .* 1e3

    # 计算深层补给
    deeprech = zeros(size(wtd))

    for j in (js+1):(je-1)
        for i in 1:imax
            if landmask[i, j] == 1
                if wtd[i, j] < slz[1] - dz[1]
                    # 计算排水导水率
                    nsoil = soiltxt[1, i, j]
                    soil_params = get_soil_params(nsoil)

                    wgpmid = 0.5 * (smoiwtd[i, j] + soil_params.θ_sat)
                    kfup = soil_params.K_sat *
                           (wgpmid / soil_params.θ_sat)^(2.0 * soil_params.b + 3.0)

                    # 计算湿度势
                    vt3dbdw = soil_params.ψ *
                              (soil_params.θ_sat / smoiwtd[i, j])^soil_params.b

                    # 计算通量(=补给)
                    deeprech[i, j] = Δt * kfup *
                                     ((soil_params.ψ - vt3dbdw) / (slz[1] - wtd[i, j]) - 1.0)

                    # 更新水位处土壤湿度
                    newwgp = smoiwtd[i, j] + (deeprech[i, j] - bottomflux[i, j]) / (slz[1] - wtd[i, j])

                    if newwgp < soil_params.ρb
                        deeprech[i, j] += (soil_params.ρb - newwgp) * (slz[1] - wtd[i, j])
                        newwgp = soil_params.ρb
                    end

                    if newwgp > soil_params.θ_sat
                        deeprech[i, j] -= (soil_params.θ_sat - newwgp) * (slz[1] - wtd[i, j])
                        newwgp = soil_params.θ_sat
                    end

                    smoiwtd[i, j] = newwgp
                    rech[i, j] += deeprech[i, j] * 1e3
                end
            end
        end
    end

    # 重置底部通量
    bottomflux .= 0.0

    # 更新水位和土壤湿度
    for j in (js+1):(je-1)
        for i in 1:imax
            if landmask[i, j] == 1
                # 地下水平衡总量
                totwater = qlat[i, j] - deeprech[i, j]

                qspring = 0.0
                updatewtd!(nzg, slz, dz, wtd[i, j], qspring, totwater,
                    view(smoi, :, i, j), view(smoieq, :, i, j),
                    soiltxt[:, i, j], smoiwtd[i, j])

                qsprings[i, j] += qspring * 1e3
            end
        end
    end

    return nothing
end

"""
    lateral_flow!(网格参数...)

计算地下水侧向流

使用三角形差分格式计算8方向邻居的流量
"""
function lateral_flow!(
    imax::Int, jmax::Int, is::Int, ie::Int, js::Int, je::Int,
    wtd::M, qlat::M, fdepth::M,
    topo::M, landmask::Matrix{Int}, Δt::T,
    area::M, κlat::M
) where {T<:AbstractFloat,M<:Matrix{T}}
    # 流动角度因子
    fangle = sqrt(tan(π₄ / 32.0)) / (2.0 * sqrt(2.0))

    # 计算有效导水率和水头
    κcell = zeros(size(wtd))
    head = zeros(size(wtd))

    for j in max(js, 1):min(je, jmax)
        for i in max(is, 1):min(ie, imax)
            if fdepth[i, j] < 1e-6
                κcell[i, j] = 0.0
            elseif wtd[i, j] < -1.5
                κcell[i, j] = fdepth[i, j] * κlat[i, j] * exp((wtd[i, j] + 1.5) / fdepth[i, j])
            else
                κcell[i, j] = κlat[i, j] * (wtd[i, j] + 1.5 + fdepth[i, j])
            end

            head[i, j] = topo[i, j] + wtd[i, j]
        end
    end

    # 计算侧向流
    for j in (js+1):(je-1)
        for i in (is+1):(ie-1)
            if landmask[i, j] == 1
                q = 0.0

                # 8个方向的流量计算
                # 对角线方向 (除以sqrt(2))
                q += (κcell[i-1, j+1] + κcell[i, j]) * (head[i-1, j+1] - head[i, j]) / sqrt(2.0)
                q += (κcell[i-1, j-1] + κcell[i, j]) * (head[i-1, j-1] - head[i, j]) / sqrt(2.0)
                q += (κcell[i+1, j+1] + κcell[i, j]) * (head[i+1, j+1] - head[i, j]) / sqrt(2.0)
                q += (κcell[i+1, j-1] + κcell[i, j]) * (head[i+1, j-1] - head[i, j]) / sqrt(2.0)

                # 正交方向
                q += (κcell[i-1, j] + κcell[i, j]) * (head[i-1, j] - head[i, j])
                q += (κcell[i, j+1] + κcell[i, j]) * (head[i, j+1] - head[i, j])
                q += (κcell[i, j-1] + κcell[i, j]) * (head[i, j-1] - head[i, j])
                q += (κcell[i+1, j] + κcell[i, j]) * (head[i+1, j] - head[i, j])

                qlat[i, j] = fangle * q * Δt / area[i, j]
            end
        end
    end

    return nothing
end

"""
    updatewtd!(水位更新参数...)

更新水位深度和土壤湿度分布

处理正流量(水位上升)和负流量(水位下降)两种情况
"""
function updatewtd!(
    nzg::Int, slz::V, dz::V,
    wtd::T, qspring::T, totwater::T,
    smoi::AbstractVector{T}, smoieq::AbstractVector{T},
    soiltextures::Vector{Int}, smoiwtd::T
) where {T<:AbstractFloat,V<:Vector{T}}
    # 确定土壤类型分布
    soiltxt = similar(slz, Int)
    for k in 1:nzg
        if slz[k] < -0.3
            soiltxt[k] = soiltextures[1]
        else
            soiltxt[k] = soiltextures[2]
        end
    end

    iwtd = 1
    qspring = 0.0

    if totwater > 0.0  # 水位上升情况
        if wtd >= slz[1]
            # 找到水位所在层
            k = 2
            while k <= nzg && wtd < slz[k]
                k += 1
            end
            iwtd = k
            kwtd = iwtd - 1
            nsoil = soiltxt[kwtd]
            soil_params = get_soil_params(nsoil)

            # 该层最大可容纳水量
            maxwatup = dz[kwtd] * (soil_params.θ_sat - smoi[kwtd])

            if totwater <= maxwatup
                smoi[kwtd] += totwater / dz[kwtd]
                smoi[kwtd] = min(smoi[kwtd], soil_params.θ_sat)

                if smoi[kwtd] > smoieq[kwtd]
                    wtd = min((smoi[kwtd] * dz[kwtd] - smoieq[kwtd] * slz[iwtd] +
                               soil_params.θ_sat * slz[kwtd]) /
                              (soil_params.θ_sat - smoieq[kwtd]), slz[iwtd])
                end
                totwater = 0.0
            else
                # 水量足够饱和该层，继续向上层填充
                smoi[kwtd] = soil_params.θ_sat
                totwater -= maxwatup

                for k in iwtd:(nzg+1)
                    wtd = slz[k]
                    iwtd = k + 1
                    if k == nzg + 1
                        break
                    end

                    nsoil = soiltxt[k]
                    soil_params = get_soil_params(nsoil)
                    maxwatup = dz[k] * (soil_params.θ_sat - smoi[k])

                    if totwater <= maxwatup
                        smoi[k] += totwater / dz[k]
                        smoi[k] = min(smoi[k], soil_params.θ_sat)

                        if smoi[k] > smoieq[k]
                            wtd = min((smoi[k] * dz[k] - smoieq[k] * slz[iwtd] +
                                       soil_params.θ_sat * slz[k]) /
                                      (soil_params.θ_sat - smoieq[k]), slz[iwtd])
                        end
                        totwater = 0.0
                        break
                    else
                        smoi[k] = soil_params.θ_sat
                        totwater -= maxwatup
                    end
                end
            end

        elseif wtd >= slz[1] - dz[1]  # 水位在土壤模型底部以下
            nsoil = soiltxt[1]
            soil_params = get_soil_params(nsoil)
            maxwatup = (soil_params.θ_sat - smoiwtd) * dz[1]

            if totwater <= maxwatup
                smoieqwtd = soil_params.θ_sat *
                            (soil_params.ψ / (soil_params.ψ - dz[1]))^(1.0 / soil_params.b)
                smoieqwtd = max(smoieqwtd, soil_params.ρb)

                smoiwtd += totwater / dz[1]
                smoiwtd = min(smoiwtd, soil_params.θ_sat)

                if smoiwtd > smoieqwtd
                    wtd = min((smoiwtd * dz[1] - smoieqwtd * slz[1] +
                               soil_params.θ_sat * (slz[1] - dz[1])) /
                              (soil_params.θ_sat - smoieqwtd), slz[1])
                end
                totwater = 0.0
            else
                # 继续处理超出部分...
                smoiwtd = soil_params.θ_sat
                totwater -= maxwatup
                # 向上层传播...
            end

        else  # 深部水位
            nsoil = soiltxt[1]
            soil_params = get_soil_params(nsoil)
            maxwatup = (soil_params.θ_sat - smoiwtd) * (slz[1] - dz[1] - wtd)

            if totwater <= maxwatup
                wtd += totwater / (soil_params.θ_sat - smoiwtd)
                totwater = 0.0
            else
                # 水位上升到土壤层...
                totwater -= maxwatup
                wtd = slz[1] - dz[1]
                # 继续处理...
            end
        end

        # 地表出水
        qspring = totwater

    elseif totwater < 0.0  # 水位下降情况

        if wtd >= slz[1]  # 水位在已解析层内
            # 找到水位层
            k = 2
            while k <= nzg && wtd < slz[k]
                k += 1
            end
            iwtd = k

            # 从水位层开始向下排水
            for kwtd in (iwtd-1):-1:1
                nsoil = soiltxt[kwtd]
                soil_params = get_soil_params(nsoil)

                # 该层最大可产水量
                maxwatdw = dz[kwtd] * (smoi[kwtd] - smoieq[kwtd])

                if -totwater <= maxwatdw
                    smoi[kwtd] += totwater / dz[kwtd]

                    if smoi[kwtd] > smoieq[kwtd]
                        wtd = (smoi[kwtd] * dz[kwtd] - smoieq[kwtd] * slz[iwtd] +
                               soil_params.θ_sat * slz[kwtd]) /
                              (soil_params.θ_sat - smoieq[kwtd])
                    else
                        wtd = slz[kwtd]
                        iwtd -= 1
                    end
                    totwater = 0.0
                    break
                else
                    wtd = slz[kwtd]
                    iwtd -= 1
                    if maxwatdw >= 0.0
                        smoi[kwtd] = smoieq[kwtd]
                        totwater += maxwatdw
                    end
                end
            end

            # 如果到达底层仍有水量需要移除
            if iwtd == 1 && totwater < 0.0
                nsoil = soiltxt[1]
                soil_params = get_soil_params(nsoil)

                smoieqwtd = soil_params.θ_sat *
                            (soil_params.ψ / (soil_params.ψ - dz[1]))^(1.0 / soil_params.b)
                smoieqwtd = max(smoieqwtd, soil_params.ρb)

                maxwatdw = dz[1] * (smoiwtd - smoieqwtd)

                if -totwater <= maxwatdw
                    smoiwtd += totwater / dz[1]
                    wtd = max((smoiwtd * dz[1] - smoieqwtd * slz[1] +
                               soil_params.θ_sat * (slz[1] - dz[1])) /
                              (soil_params.θ_sat - smoieqwtd), slz[1] - dz[1])
                else
                    wtd = slz[1] - dz[1]
                    smoiwtd += totwater / dz[1]
                    # 进一步向下...
                    dzup = (smoieqwtd - smoiwtd) * dz[1] / (soil_params.θ_sat - smoieqwtd)
                    wtd -= dzup
                    smoiwtd = smoieqwtd
                end
            end

        else  # 水位已经在底部以下
            # 处理深部水位变化...
            nsoil = soiltxt[1]
            soil_params = get_soil_params(nsoil)

            wgpmid = soil_params.θ_sat *
                     (soil_params.ψ / (soil_params.ψ - (slz[1] - wtd)))^(1.0 / soil_params.b)
            wgpmid = max(wgpmid, soil_params.ρb)

            syielddw = soil_params.θ_sat - wgpmid
            wtdold = wtd
            wtd = wtdold + totwater / syielddw

            # 更新水位处湿度
            smoiwtd = (smoiwtd * (slz[1] - wtdold) + wgpmid * (wtdold - wtd)) / (slz[1] - wtd)
        end

        qspring = 0.0
    end

    return wtd, qspring, smoiwtd
end
