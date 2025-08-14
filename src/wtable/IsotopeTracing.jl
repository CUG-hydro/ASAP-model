# 同位素追踪模块
export lateral_isotope!, updatedeepwtable!

"""
    lateral_isotope!(同位素侧向流计算参数...)

计算地下水侧向流的同位素传输

参数：
- 网格参数: imax, jmax, is, ie, js, je, nzg
- 空间参数: soiltxt, wtd, fdepth, topo, landmask, area, lats, dxy
- 同位素参数: o18, o18wtd
- 流量参数: qlat, qlato18, qlatin, qlatout, qlatino18, qlatouto18
- 累积参数: qlatinsum, qlatoutsum, qlatino18sum, qlatouto18sum
- 时间参数: Δt
"""
function lateral_isotope!(
  imax::Int, jmax::Int, is::Int, ie::Int, js::Int, je::Int, nzg::Int,
  soiltxt::Array{Int,3}, wtd::M, qlat::M,
  fdepth::M, topo::M, landmask::Matrix{Int},
  Δt::T, area::M, lats::M, dxy::T,
  slz::Vector{T}, o18::Array{T,3}, smoi::Array{T,3},
  qlato18::M, qlatin::M, qlatout::M,
  qlatino18::M, qlatouto18::M,
  qlatinsum::M, qlatoutsum::M,
  qlatino18sum::M, qlatouto18sum::M
) where {T<:AbstractFloat,M<:Matrix{T}}
  # 计算侧向流导水率
  κlat = zeros(size(wtd))
  o18wtd = zeros(size(wtd))

  for j in max(js, 1):min(je, jmax)
    for i in max(is, 1):min(ie, imax)
      nsoil = soiltxt[1, i, j]
      soil_params = get_soil_params(nsoil)
      κlat[i, j] = soil_params.Ksat * soil_params.K_latfactor

      # 确定水位所在层的同位素浓度
      k = 1
      while k <= nzg && wtd[i, j] < slz[k]
        k += 1
      end
      kwtd = max(k - 1, 1)
      o18wtd[i, j] = o18[kwtd, i, j] / smoi[kwtd, i, j]
    end
  end

  # 计算带同位素的侧向流
  lateralflow_with_isotope!(
    imax, jmax, is, ie, js, je, wtd, qlat, fdepth, topo, landmask,
    Δt, area, κlat, lats, dxy, o18wtd, qlato18, qlatin, qlatout,
    qlatino18, qlatouto18
  )

  # 累积流量 (转换为mm)
  qlatinsum .+= qlatin .* 1e3
  qlatoutsum .+= qlatout .* 1e3
  qlatino18sum .+= qlatino18 .* 1e3
  qlatouto18sum .+= qlatouto18 .* 1e3

  return nothing
end


"""
    lateralflow_with_isotope!(带同位素的侧向流计算)

考虑同位素传输的四方向侧向流计算
"""
function lateralflow_with_isotope!(
  imax::Int, jmax::Int, is::Int, ie::Int, js::Int, je::Int,
  wtd::M, qlat::M, fdepth::M,
  topo::M, landmask::Matrix{Int}, Δt::T,
  area::M, κlat::M, lats::M,
  dxy::T, o18wtd::M, qlato18::M,
  qlatin::M, qlatout::M, qlatino18::M,
  qlatouto18::M
) where {T<:AbstractFloat,M<:Matrix{T}}
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

  # 初始化输出数组
  qlat .= 0.0
  qlato18 .= 0.0
  qlatin .= 0.0
  qlatout .= 0.0
  qlatino18 .= 0.0
  qlatouto18 .= 0.0

  # 计算四方向侧向流
  for j in max(js + 1, 2):min(je - 1, jmax - 1)
    for i in max(is + 1, 2):min(ie - 1, imax - 1)
      if landmask[i, j] == 1
        q = 0.0
        qo18 = 0.0
        qin = 0.0
        qout = 0.0
        qino18 = 0.0
        qouto18 = 0.0

        # 北方向流量
        qn = (κcell[i, j+1] + κcell[i, j]) * (head[i, j+1] - head[i, j]) *
             cos(π2r * (lats[i, j] + 0.5 * dxy))

        if qn > 0.0
          qno18 = qn * o18wtd[i, j+1]
          qin += qn
          qino18 += qno18
        else
          qno18 = qn * o18wtd[i, j]
          qout += qn
          qouto18 += qno18
        end

        # 南方向流量
        qs = (κcell[i, j-1] + κcell[i, j]) * (head[i, j-1] - head[i, j]) *
             cos(π2r * (lats[i, j] - 0.5 * dxy))

        if qs > 0.0
          qso18 = qs * o18wtd[i, j-1]
          qin += qs
          qino18 += qso18
        else
          qso18 = qs * o18wtd[i, j]
          qout += qs
          qouto18 += qso18
        end

        # 西方向流量
        qw = (κcell[i-1, j] + κcell[i, j]) * (head[i-1, j] - head[i, j]) /
             cos(π2r * lats[i, j])

        if qw > 0.0
          qwo18 = qw * o18wtd[i-1, j]
          qin += qw
          qino18 += qwo18
        else
          qwo18 = qw * o18wtd[i, j]
          qout += qw
          qouto18 += qwo18
        end

        # 东方向流量
        qe = (κcell[i+1, j] + κcell[i, j]) * (head[i+1, j] - head[i, j]) /
             cos(π2r * lats[i, j])

        if qe > 0.0
          qeo18 = qe * o18wtd[i+1, j]
          qin += qe
          qino18 += qeo18
        else
          qeo18 = qe * o18wtd[i, j]
          qout += qe
          qouto18 += qeo18
        end

        # 总流量
        q = qn + qs + qw + qe
        qlat[i, j] = 0.5 * q * Δt / area[i, j]

        # 同位素流量
        qo18 = qno18 + qso18 + qwo18 + qeo18
        qlato18[i, j] = 0.5 * qo18 * Δt / area[i, j]

        # 入流和出流
        qlatin[i, j] = 0.5 * qin * Δt / area[i, j]
        qlatout[i, j] = -0.5 * qout * Δt / area[i, j]

        qlatino18[i, j] = 0.5 * qino18 * Δt / area[i, j]
        qlatouto18[i, j] = -0.5 * qouto18 * Δt / area[i, j]
      end
    end
  end

  return nothing
end

"""
    updatedeepwtable!(深层水位更新参数...)

更新深层地下水位和计算补给量
"""
function updatedeepwtable!(
  imax::Int, jmax::Int, js::Int, je::Int, nzg::Int,
  slz::Vector{T}, dz::Vector{T}, soiltxt::Array{Int,3},
  wtd::M, bottomflux::M, rech::M,
  qslat::M, qlat::M, landmask::Matrix{Int},
  Δt::T, smoi::Array{T,3}, smoieq::Array{T,3},
  smoiwtd::M, qsprings::M
) where {T<:AbstractFloat,M<:Matrix{T}}
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
          kfup = soil_params.Ksat *
                 (wgpmid / soil_params.θ_sat)^(2.0 * soil_params.b + 3.0)

          # 计算湿度势
          vt3dbdw = soil_params.ψsat *
                    (soil_params.θ_sat / smoiwtd[i, j])^soil_params.b

          # 计算通量(=补给)
          deeprech[i, j] = Δt * kfup *
                           ((soil_params.ψsat - vt3dbdw) / (slz[1] - wtd[i, j]) - 1.0)

          # 更新水位处土壤湿度
          newwgp = smoiwtd[i, j] + (deeprech[i, j] - bottomflux[i, j]) / (slz[1] - wtd[i, j])

          if newwgp < soil_params.θ_cp
            deeprech[i, j] += (soil_params.θ_cp - newwgp) * (slz[1] - wtd[i, j])
            newwgp = soil_params.θ_cp
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

  # 更新水位
  for j in (js+1):(je-1)
    for i in 1:imax
      if landmask[i, j] == 1
        # 地下水平衡总量
        totwater = qlat[i, j] - qslat[i, j] - deeprech[i, j]

        qspring = 0.0
        wtd_new, qspring_new, smoiwtd_new = updatewtd_simple(
          nzg, slz, dz, wtd[i, j], qspring, totwater,
          view(smoi, :, i, j), view(smoieq, :, i, j),
          soiltxt[:, i, j], smoiwtd[i, j]
        )

        wtd[i, j] = wtd_new
        smoiwtd[i, j] = smoiwtd_new
        qsprings[i, j] += qspring_new * 1e3
      end
    end
  end

  return nothing
end

"""
    updatewtd_simple(简化的水位更新函数)

简化版本的水位更新，用于深层水位计算
"""
function updatewtd_simple(
  nzg::Int, slz::V, dz::V,
  wtd::T, qspring::T, totwater::T,
  smoi::AbstractVector{T}, smoieq::AbstractVector{T},
  soiltextures::Vector{Int}, smoiwtd::T
) where {T<:AbstractFloat,V<:Vector{T}}
  # 这是一个简化版本，主要处理深层水位变化
  # 完整实现请参考WaterTableCalculations模块中的updatewtd!函数

  new_wtd = wtd
  new_qspring = 0.0
  new_smoiwtd = smoiwtd

  if totwater != 0.0
    # 简单的线性调整
    if totwater > 0.0
      # 水位上升
      nsoil = soiltextures[1]
      soil_params = get_soil_params(nsoil)

      if wtd < slz[1] - dz[1]
        # 深部水位上升
        maxwatup = (soil_params.θ_sat - smoiwtd) * (slz[1] - dz[1] - wtd)
        if totwater <= maxwatup
          new_wtd = wtd + totwater / (soil_params.θ_sat - smoiwtd)
        else
          new_wtd = slz[1] - dz[1]
          new_qspring = totwater - maxwatup
        end
      else
        new_qspring = totwater
      end
    else
      # 水位下降
      nsoil = soiltextures[1]
      soil_params = get_soil_params(nsoil)

      if wtd >= slz[1] - dz[1]
        # 在土壤层附近
        maxwatdw = dz[1] * (smoiwtd - soil_params.θ_cp)
        if -totwater <= maxwatdw
          new_smoiwtd = smoiwtd + totwater / dz[1]
        else
          new_wtd = wtd + totwater / (soil_params.θ_sat - soil_params.θ_cp)
        end
      else
        # 深部水位下降
        wgpmid = soil_params.θ_sat *
                 (soil_params.ψsat / (soil_params.ψsat - (slz[1] - wtd)))^(1.0 / soil_params.b)
        wgpmid = max(wgpmid, soil_params.θ_cp)

        syielddw = soil_params.θ_sat - wgpmid
        new_wtd = wtd + totwater / syielddw
      end
    end
  end
  return new_wtd, new_qspring, new_smoiwtd
end
