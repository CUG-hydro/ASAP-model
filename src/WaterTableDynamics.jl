"""
地下水位更新模块
计算浅层地下水位的变化
"""
module WaterTableDynamics

using ..SoilParameters

export updateshallowwtd, updatewtdqlat

"""
更新浅层地下水位

# 参数
- `i::Int`, `j::Int`: 网格坐标
- `nzg::Int`: 土壤层数
- `freedrain::Int`: 自由排水标志
- `slz::Vector{Float64}`: 土壤层深度 (m) [nzg+1]
- `dz::Vector{Float64}`: 土壤层厚度 (m) [nzg]
- `soiltxt::Int`: 土壤类型
- `smoieq::Vector{Float64}`: 平衡含水量 [nzg]
- `smoiwtd::Float64`: 地下水含水量
- `smoi::Vector{Float64}`: 土壤含水量 [nzg]
- `wtd::Float64`: 地下水位深度 (m)
- `fdepth::Float64`: 根系深度因子 (m)

# 返回
- `Tuple{Float64, Float64}`: (更新的地下水位, 补给量)
"""
function updateshallowwtd(i::Int, j::Int, nzg::Int, freedrain::Int, slz::Vector{Float64},
                         dz::Vector{Float64}, soiltxt::Int, smoieq::Vector{Float64},
                         smoiwtd::Float64, smoi::Vector{Float64}, wtd::Float64, fdepth::Float64)
    
    rech = 0.0
    vctr4 = 0.5 .* (slz[1:nzg] .+ slz[2:nzg+1])  # 层中心深度
    
    # 确定地下水位所在层
    iwtd = 1
    for k in 1:nzg
        if wtd < slz[k]
            iwtd = k
            break
        end
    end
    
    while true
        flag = 0
        kwtd = iwtd - 1
        
        if kwtd > 0  # 地下水位在已解析层中
            wtdold = wtd
            
            # 获取土壤参数
            soil_params = get_soil_params(soiltxt)
            
            # 检查上层是否接近饱和
            if kwtd > 1
                smoisat = soil_params.slmsts * max(min(exp((vctr4[kwtd-1] + 1.5) / fdepth), 1.0), 0.1)
                if wtd < slz[kwtd] + 0.01 && smoi[kwtd-1] < smoisat
                    flag = 1
                end
            end
            
            smoisat = soil_params.slmsts * max(min(exp((vctr4[kwtd] + 1.5) / fdepth), 1.0), 0.1)
            
            if smoi[kwtd] > smoieq[kwtd] && flag == 0
                if smoi[kwtd] == smoisat  # 地下水位上升到上层
                    wtd = slz[iwtd]
                    rech = (wtdold - wtd) * (smoisat - smoieq[kwtd])
                    iwtd = iwtd + 1
                    kwtd = kwtd + 1
                    
                    if kwtd <= nzg
                        if smoi[kwtd] > smoieq[kwtd]
                            wtdold = wtd
                            smoisat = soil_params.slmsts * max(min(exp((vctr4[kwtd] + 1.5) / fdepth), 1.0), 0.1)
                            wtd = min((smoi[kwtd] * dz[kwtd] - smoieq[kwtd] * slz[iwtd] + smoisat * slz[kwtd]) /
                                     (smoisat - smoieq[kwtd]), slz[iwtd])
                            rech = rech + (wtdold - wtd) * (smoisat - smoieq[kwtd])
                        end
                    else
                        break
                    end
                else  # 地下水位在层内
                    wtd = min((smoi[kwtd] * dz[kwtd] - smoieq[kwtd] * slz[iwtd] + smoisat * slz[kwtd]) /
                             (smoisat - smoieq[kwtd]), slz[iwtd])
                    rech = (wtdold - wtd) * (smoisat - smoieq[kwtd])
                    break
                end
                
            else  # 地下水位下降到下层
                wtd = slz[kwtd]
                rech = (wtdold - wtd) * (smoisat - smoieq[kwtd])
                kwtd = kwtd - 1
                iwtd = iwtd - 1
                
                # 调整到下层
                if kwtd >= 1
                    wtdold = wtd
                    smoisat = soil_params.slmsts * max(min(exp((vctr4[kwtd] + 1.5) / fdepth), 1.0), 0.1)
                    
                    if smoi[kwtd] > smoieq[kwtd]
                        wtd = min((smoi[kwtd] * dz[kwtd] - smoieq[kwtd] * slz[iwtd] + smoisat * slz[kwtd]) /
                                 (smoisat - smoieq[kwtd]), slz[iwtd])
                        rech = rech + (wtdold - wtd) * (smoisat - smoieq[kwtd])
                        break
                    else
                        wtd = slz[kwtd]
                        rech = rech + (wtdold - wtd) * (smoisat - smoieq[kwtd])
                    end
                else
                    break
                end
            end
        else
            break
        end
    end
    
    if wtd < slz[1]
        println("地下水位问题: wtd=$wtd, i=$i, j=$j")
    end
    
    return wtd, rech
end

"""
通过侧向流更新地下水位

# 参数
- `nzg::Int`: 土壤层数
- `slz::Vector{Float64}`: 土壤层深度 (m) [nzg+1]
- `dz::Vector{Float64}`: 土壤层厚度 (m) [nzg]
- `wtd::Float64`: 地下水位深度 (m)
- `qlat::Float64`: 侧向流 (m)
- `smoi::Vector{Float64}`: 土壤含水量 [nzg]
- `smoieq::Vector{Float64}`: 平衡含水量 [nzg]
- `soiltextures::Vector{Int}`: 土壤类型 [2]
- `smoiwtd::Float64`: 地下水含水量
- `qlatflux::Vector{Float64}`: 侧向流通量 [nzg+2]
- `fdepth::Float64`: 根系深度因子 (m)

# 返回
- `Tuple{Float64, Float64, Vector{Float64}, Vector{Float64}}`: 
  (泉水流量, 更新的地下水位, 更新的土壤含水量, 更新的侧向流通量)
"""
function updatewtdqlat(nzg::Int, slz::Vector{Float64}, dz::Vector{Float64}, wtd::Float64,
                      qlat::Float64, smoi::Vector{Float64}, smoieq::Vector{Float64},
                      soiltextures::Vector{Int}, smoiwtd::Float64, qlatflux::Vector{Float64},
                      fdepth::Float64)
    
    vctr4 = 0.5 .* (slz[1:nzg] .+ slz[2:nzg+1])  # 层中心深度
    
    # 确定各层土壤类型
    soiltxt = zeros(Int, nzg)
    for k in 1:nzg
        soiltxt[k] = slz[k] < -0.3 ? soiltextures[1] : soiltextures[2]
    end
    
    qspring = 0.0
    totwater = qlat
    
    iwtd = 1
    
    # 情况1: totwater > 0 (地下水位上升)
    if totwater > 0.0
        # 确定当前地下水位层
        for k in 2:nzg
            if wtd < slz[k]
                iwtd = k
                break
            end
        end
        
        kwtd = iwtd - 1
        soil_params = get_soil_params(soiltxt[kwtd])
        smoisat = soil_params.slmsts * max(min(exp((vctr4[kwtd] + 1.5) / fdepth), 1.0), 0.1)
        
        # 该层可容纳的最大水量
        maxwatup = dz[kwtd] * (smoisat - smoi[kwtd])
        
        if totwater <= maxwatup
            smoi[kwtd] = smoi[kwtd] + totwater / dz[kwtd]
            qlatflux[kwtd] = qlatflux[kwtd] + totwater
            smoi[kwtd] = min(smoi[kwtd], smoisat)
            
            if smoi[kwtd] > smoieq[kwtd]
                wtd = min((smoi[kwtd] * dz[kwtd] - smoieq[kwtd] * slz[iwtd] + smoisat * slz[kwtd]) /
                         (smoisat - smoieq[kwtd]), slz[iwtd])
            end
            totwater = 0.0
        else  # 水量足够饱和该层
            smoi[kwtd] = smoisat
            qlatflux[kwtd] = qlatflux[kwtd] + maxwatup
            totwater = totwater - maxwatup
            
            k1 = iwtd
            for k in k1:(nzg+1)
                wtd = slz[k]
                iwtd = k + 1
                if k == nzg + 1
                    break
                end
                
                soil_params = get_soil_params(soiltxt[k])
                smoisat = soil_params.slmsts * max(min(exp((vctr4[k] + 1.5) / fdepth), 1.0), 0.1)
                maxwatup = dz[k] * (smoisat - smoi[k])
                
                if totwater <= maxwatup
                    smoi[k] = smoi[k] + totwater / dz[k]
                    qlatflux[k] = qlatflux[k] + totwater
                    smoi[k] = min(smoi[k], smoisat)
                    
                    if smoi[k] > smoieq[k]
                        wtd = min((smoi[k] * dz[k] - smoieq[k] * slz[iwtd] + smoisat * slz[k]) /
                                 (smoisat - smoieq[k]), slz[iwtd])
                    end
                    totwater = 0.0
                    break
                else
                    smoi[k] = smoisat
                    qlatflux[k] = qlatflux[k] + maxwatup
                    totwater = totwater - maxwatup
                end
            end
        end
        
        # 地表泉水
        qspring = totwater
        
    # 情况2: totwater < 0 (地下水位下降)
    elseif totwater < 0.0
        # 确定当前地下水位层
        for k in 2:nzg
            if wtd < slz[k]
                iwtd = k
                break
            end
        end
        
        k1 = iwtd - 1
        for kwtd in k1:-1:1
            soil_params = get_soil_params(soiltxt[kwtd])
            smoisat = soil_params.slmsts * max(min(exp((vctr4[kwtd] + 1.5) / fdepth), 1.0), 0.1)
            
            # 该层可提供的最大水量
            maxwatdw = dz[kwtd] * (smoi[kwtd] - smoieq[kwtd])
            
            if -totwater <= maxwatdw
                smoi[kwtd] = smoi[kwtd] + totwater / dz[kwtd]
                qlatflux[kwtd] = qlatflux[kwtd] + totwater
                
                if smoi[kwtd] > smoieq[kwtd]
                    wtd = (smoi[kwtd] * dz[kwtd] - smoieq[kwtd] * slz[iwtd] + smoisat * slz[kwtd]) /
                          (smoisat - smoieq[kwtd])
                else
                    wtd = slz[kwtd]
                    iwtd = iwtd - 1
                end
                totwater = 0.0
                break
            else
                wtd = slz[kwtd]
                iwtd = iwtd - 1
                if maxwatdw >= 0.0
                    smoi[kwtd] = smoieq[kwtd]
                    qlatflux[kwtd] = qlatflux[kwtd] + maxwatdw
                    totwater = totwater + maxwatdw
                end
            end
        end
        
        if iwtd == 1 && totwater < 0.0
            soil_params = get_soil_params(soiltxt[1])
            smoisat = soil_params.slmsts * max(min(exp((vctr4[1] + 1.5) / fdepth), 1.0), 0.1)
            
            smoi[1] = smoi[1] + totwater / dz[1]
            qlatflux[1] = qlatflux[1] + totwater
            wtd = max((smoi[1] * dz[1] - smoieq[1] * slz[2] + smoisat * slz[1]) /
                     (smoisat - smoieq[1]), slz[1])
        end
        
        qspring = 0.0
    end
    
    return qspring, wtd, smoi, qlatflux
end

end # module WaterTableDynamics
