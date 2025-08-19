"""
    updatewtd!(水位更新参数...)

更新水位深度和土壤湿度分布

处理正流量(水位上升)和负流量(水位下降)两种情况
"""
function updatewtd!(
  nzg::Int, slz::V, dz::V,
  wtd::T, qspring::T, totwater::T,
  θ::AbstractVector{T}, θ_eq::AbstractVector{T},
  soiltextures::Vector{Int}, θ_wtd::T
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
      soil = get_soil_params(nsoil)
      (; θ_sat) = soil

      # 该层最大可容纳水量
      maxwatup = dz[kwtd] * (θ_sat - θ[kwtd])

      if totwater <= maxwatup
        θ[kwtd] += totwater / dz[kwtd]
        θ[kwtd] = min(θ[kwtd], θ_sat)

        if θ[kwtd] > θ_eq[kwtd]
          wtd = min((θ[kwtd] * dz[kwtd] - θ_eq[kwtd] * slz[iwtd] +
                     θ_sat * slz[kwtd]) /
                    (θ_sat - θ_eq[kwtd]), slz[iwtd])
        end
        totwater = 0.0
      else
        # 水量足够饱和该层，继续向上层填充
        θ[kwtd] = soil.θ_sat
        totwater -= maxwatup

        for k in iwtd:(nzg+1)
          wtd = slz[k]
          iwtd = k + 1
          if k == nzg + 1
            break
          end

          nsoil = soiltxt[k]
          soil = get_soil_params(nsoil)
          maxwatup = dz[k] * (soil.θ_sat - θ[k])

          if totwater <= maxwatup
            θ[k] += totwater / dz[k]
            θ[k] = min(θ[k], soil.θ_sat)

            if θ[k] > θ_eq[k]
              wtd = min((θ[k] * dz[k] - θ_eq[k] * slz[iwtd] +
                         soil.θ_sat * slz[k]) /
                        (soil.θ_sat - θ_eq[k]), slz[iwtd])
            end
            totwater = 0.0
            break
          else
            θ[k] = soil.θ_sat
            totwater -= maxwatup
          end
        end
      end

    elseif wtd >= slz[1] - dz[1]  # 水位在土壤模型底部以下
      nsoil = soiltxt[1]
      soil = get_soil_params(nsoil)
      maxwatup = (soil.θ_sat - θ_wtd) * dz[1]

      if totwater <= maxwatup
        smoieqwtd = soil.θ_sat *
                    (soil.ψsat / (soil.ψsat - dz[1]))^(1.0 / soil.b)
        smoieqwtd = max(smoieqwtd, soil.θ_cp)

        θ_wtd += totwater / dz[1]
        θ_wtd = min(θ_wtd, soil.θ_sat)

        if θ_wtd > smoieqwtd
          wtd = min((θ_wtd * dz[1] - smoieqwtd * slz[1] +
                     soil.θ_sat * (slz[1] - dz[1])) /
                    (soil.θ_sat - smoieqwtd), slz[1])
        end
        totwater = 0.0
      else
        # 继续处理超出部分...
        θ_wtd = soil.θ_sat
        totwater -= maxwatup
        # 向上层传播...
      end

    else  # 深部水位
      nsoil = soiltxt[1]
      soil = get_soil_params(nsoil)
      maxwatup = (soil.θ_sat - θ_wtd) * (slz[1] - dz[1] - wtd)

      if totwater <= maxwatup
        wtd += totwater / (soil.θ_sat - θ_wtd)
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
        soil = get_soil_params(nsoil)

        # 该层最大可产水量
        maxwatdw = dz[kwtd] * (θ[kwtd] - θ_eq[kwtd])

        if -totwater <= maxwatdw
          θ[kwtd] += totwater / dz[kwtd]

          if θ[kwtd] > θ_eq[kwtd]
            wtd = (θ[kwtd] * dz[kwtd] - θ_eq[kwtd] * slz[iwtd] +
                   soil.θ_sat * slz[kwtd]) /
                  (soil.θ_sat - θ_eq[kwtd])
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
            θ[kwtd] = θ_eq[kwtd]
            totwater += maxwatdw
          end
        end
      end

      # 如果到达底层仍有水量需要移除
      if iwtd == 1 && totwater < 0.0
        nsoil = soiltxt[1]
        soil = get_soil_params(nsoil)

        smoieqwtd = soil.θ_sat *
                    (soil.ψsat / (soil.ψsat - dz[1]))^(1.0 / soil.b)
        smoieqwtd = max(smoieqwtd, soil.θ_cp)

        maxwatdw = dz[1] * (θ_wtd - smoieqwtd)

        if -totwater <= maxwatdw
          θ_wtd += totwater / dz[1]
          wtd = max((θ_wtd * dz[1] - smoieqwtd * slz[1] +
                     soil.θ_sat * (slz[1] - dz[1])) /
                    (soil.θ_sat - smoieqwtd), slz[1] - dz[1])
        else
          wtd = slz[1] - dz[1]
          θ_wtd += totwater / dz[1]
          # 进一步向下...
          dzup = (smoieqwtd - θ_wtd) * dz[1] / (soil.θ_sat - smoieqwtd)
          wtd -= dzup
          θ_wtd = smoieqwtd
        end
      end

    else  # 水位已经在底部以下
      # 处理深部水位变化...
      nsoil = soiltxt[1]
      soil = get_soil_params(nsoil)

      wgpmid = soil.θ_sat *
               (soil.ψsat / (soil.ψsat - (slz[1] - wtd)))^(1.0 / soil.b)
      wgpmid = max(wgpmid, soil.θ_cp)

      syielddw = soil.θ_sat - wgpmid
      wtdold = wtd
      wtd = wtdold + totwater / syielddw

      # 更新水位处湿度
      θ_wtd = (θ_wtd * (slz[1] - wtdold) + wgpmid * (wtdold - wtd)) / (slz[1] - wtd)
    end
    qspring = 0.0
  end
  return wtd, qspring, θ_wtd
end
