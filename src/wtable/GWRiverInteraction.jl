# 河流-地下水相互作用模块

export gw2river!, rivers_kw_flood!, rivers_dw_flood!, flooding!, moveqrf!


"""
    gw2river!(地下水-河流交换参数...)

计算地下水与河流之间的水量交换

参数：
- 网格参数  : imax, jmax, is, ie, js, je, nzg
- 土壤层参数: slz
- 时间参数  : Δt
- 土壤类型  : soiltxt
- 空间参数  : landmask, wtd, maxdepth, riverdepth, width, length, area, fdepth
- 输出     : qrf (河流-地下水通量)
"""
function gw2river!(
  imax::Int, jmax::Int, is::Int, ie::Int, js::Int, je::Int, nzg::Int,
  slz::V, Δt::T, soiltxt::Array{Int,3},
  landmask::Matrix{Int}, wtd::M, maxdepth::M,
  riverdepth::M, width::M, length::M,
  area::M, fdepth::M, qrf::M
) where {T<:AbstractFloat,V<:Vector{T},M<:Matrix{T}}
  qrf .= 0.0

  for j in (js+1):(je-1)
    for i in max(is + 1, 2):min(ie - 1, imax - 1)
      landmask[i, j] == 0 || width[i, j] == 0.0 && continue

      rdepth = max(riverdepth[i, j], 0.0)        # 河深
      riversurface = -(maxdepth[i, j] - rdepth)  # 河面高程，相对于fullbank
      riversurface >= 0.0 && continue

      nsoil = soiltxt[2, i, j]
      soil = get_soil_params(nsoil)

      # Miguez-Macho, 2007, Eq. 2d, 注意这里未做积分, 认为dl = dh
      # 这里做了很多简化, 理论情况只有losing river适用该公式
      K_rb = soil.Ksat * clamp(exp((-maxdepth[i, j] + 1.5) / fdepth[i, j]), 0.1, 1.0) # river bed
      rcond = width[i, j] * length[i, j] * K_rb

      if wtd[i, j] > riversurface # 地下水位高于河流水面，向河流排水
        qrf[i, j] = rcond * (wtd[i, j] - riversurface) * (Δt / area[i, j])
        # 限制突然下降，最大50mm/day
        qrf[i, j] = min(qrf[i, j], Δt * 0.05 / 86400.0)

      elseif wtd[i, j] > -maxdepth[i, j]  # 水位与河流连接但低于河流水面, 反向river补给GW
        soilwatercap = -rcond * (wtd[i, j] - riversurface) * (Δt / area[i, j]) # 负号为了将其转换为正的数值

        soilwatercap = min(soilwatercap, Δt * 0.05 / 86400.0)
        qrf[i, j] = -max(min(soilwatercap, riverdepth[i, j]), 0.0) *
                    min(width[i, j] * length[i, j] / area[i, j], 1.0)
      else
        # 水位在河床以下，断开连接，仅渗透
        qrf[i, j] = -max(min(soil.Ksat * Δt, rdepth), 0.0) *
                    min(width[i, j] * length[i, j] / area[i, j], 1.0)
      end
    end
  end

  return nothing
end

"""
    rivers_kw_flood!(运动波河流洪水路由参数...)

使用运动波方程进行河流洪水路由计算
"""
function rivers_kw_flood!(
  imax::Int, jmax::Int, is::Int, ie::Int, js::Int, je::Int,
  Δt::T, dtlr::T, fd::Matrix{Int}, bfd::Matrix{Int},
  qnew::M, qs::M, qrf::M,
  delsfcwat::M, slope::M, depth::M,
  width::M, length::M, maxdepth::M,
  area::M, riverarea::M, floodarea::M,
  riverchannel::M, qmean::M, floodheight::M,
  topo::M
) where {T<:AbstractFloat,M<:Matrix{T}}
  # 计算外部输入流量
  qext = zeros(size(qnew))

  for j in max(js + 1, 2):min(je - 1, jmax - 1)
    for i in max(is + 1, 2):min(ie - 1, imax - 1)
      if fd[i, j] != 0
        qext[i, j] = (qrf[i, j] + qs[i, j] + delsfcwat[i, j]) / Δt * area[i, j]
      end
    end
  end

  # 保存当前流量
  q = copy(qnew)

  # 计算入流
  qin = zeros(size(q))

  for j in js:je
    for i in is:ie
      if fd[i, j] > 0
        i1, j1 = flowdir(fd, i, j)
        if i1 > is && i1 < ie && j1 > js && j1 < je
          qin[i1, j1] += q[i, j]
        end
      end
    end
  end

  # 更新河流深度和洪水高度
  for j in max(js + 1, 2):min(je - 1, jmax - 1)
    for i in max(is + 1, 2):min(ie - 1, imax - 1)
      if fd[i, j] != 0
        # 计算流入该网格的总流量变化
        dsnew = qin[i, j] - q[i, j]

        # 添加特定河流分流(原始Fortran中的硬编码分流)
        dsnew = apply_specific_diversions(i, j, q, dsnew)

        # 新的河流存储量
        snew = depth[i, j] * riverarea[i, j] + floodheight[i, j] * floodarea[i, j] +
               (dsnew + qext[i, j]) * dtlr

        # 在河道和洪泛区之间重新分配水量
        if snew >= riverchannel[i, j]
          floodheight[i, j] = (snew - riverchannel[i, j]) / max(area[i, j], riverarea[i, j])
          depth[i, j] = floodheight[i, j] + maxdepth[i, j]
        else
          floodheight[i, j] = 0.0
          if riverarea[i, j] > 0.0
            depth[i, j] = snew / riverarea[i, j]
          else
            depth[i, j] = 0.0
          end
        end
      end
    end
  end

  # 计算新的流量
  for j in (js+1):(je-1)
    for i in (is+1):(ie-1)
      if fd[i, j] != 0
        if width[i, j] * depth[i, j] > 1e-9 && fd[i, j] > 0
          # 使用曼宁公式计算流速
          aa = depth[i, j] * width[i, j] / (2.0 * depth[i, j] + width[i, j])
          flowwidth = width[i, j]

          if floodheight[i, j] > 0.05
            # 计算瞬时坡度
            i1, j1 = flowdir(fd, i, j)
            waterelevij = topo[i, j] - maxdepth[i, j] + depth[i, j]
            waterelevi1j1 = topo[i1, j1] - maxdepth[i1, j1] + max(depth[i1, j1], 0.0)
            slopefor = (waterelevij - waterelevi1j1) / (0.5 * (length[i, j] + length[i1, j1]))

            if bfd[i, j] > 0
              i2, j2 = flowdir(bfd, i, j)
              waterelevi2j2 = topo[i2, j2] - maxdepth[i2, j2] + max(depth[i2, j2], 0.0)
              slopeback = (waterelevi2j2 - waterelevij) / (0.5 * (length[i2, j2] + length[i, j]))
              slope_inst = 0.5 * (slopefor + slopeback)
            else
              slope_inst = slopefor
            end

            slope_inst = 0.25 * slope_inst + 0.75 * slope[i, j]
            if slope_inst < 0.0
              slope_inst = slope[i, j]
            end

            speed = (aa^(2.0 / 3.0)) * sqrt(slope_inst) / 0.03
            speed = max(min(speed, length[i, j] / dtlr), 0.01)
          else
            slope_inst = slope[i, j]
            speed = (aa^(2.0 / 3.0)) * sqrt(slope_inst) / 0.03
            speed = max(min(speed, length[i, j] / dtlr), 0.01)
          end
        else
          speed = 0.0
        end

        # 计算新流量
        qnew[i, j] = speed * depth[i, j] * width[i, j]
      else
        qnew[i, j] = 0.0
      end
    end
  end

  # 累积平均流量
  qmean .+= qnew .* dtlr
  return nothing
end


"""
    rivers_dw_flood!(扩散波河流洪水路由参数...)

使用扩散波方程进行河流洪水路由计算
"""
function rivers_dw_flood!(
  imax::Int, js::Int, je::Int, Δt::T, dtlr::T,
  fd::Matrix{Int}, bfd::Matrix{Int}, qnew::M,
  qs::M, qrf::M, delsfcwat::M,
  slope::M, depth::M, width::M,
  length::M, maxdepth::M, area::M,
  riverarea::M, floodarea::M, riverchannel::M,
  qmean::M, floodheight::M, topo::M
) where {T<:AbstractFloat,M<:Matrix{T}}
  # 计算外部输入流量
  qext = zeros(size(qnew))

  for j in (js+1):(je-1)
    for i in 2:(imax-1)
      if fd[i, j] != 0
        qext[i, j] = (qrf[i, j] + qs[i, j] + delsfcwat[i, j]) / Δt * area[i, j]
      end
    end
  end

  # 保存当前流量
  q = copy(qnew)

  # 计算入流
  qin = zeros(size(q))

  for j in js:je
    for i in 1:imax
      if fd[i, j] > 0
        i1, j1 = flowdir(fd, i, j)
        if i1 > 1 && i1 < imax && j1 > js && j1 < je
          qin[i1, j1] += q[i, j]
        end
      end
    end
  end

  # 更新河流深度和洪水高度 (类似KW方法)
  for j in (js+1):(je-1)
    for i in 2:(imax-1)
      if fd[i, j] != 0
        dsnew = qin[i, j] - q[i, j]
        dsnew = apply_specific_diversions(i, j, q, dsnew)

        snew = depth[i, j] * riverarea[i, j] + floodheight[i, j] * floodarea[i, j] +
               (dsnew + qext[i, j]) * dtlr

        if snew >= riverchannel[i, j]
          floodheight[i, j] = (snew - riverchannel[i, j]) / max(area[i, j], riverarea[i, j])
          depth[i, j] = floodheight[i, j] + maxdepth[i, j]
        else
          floodheight[i, j] = 0.0
          if riverarea[i, j] > 0.0
            depth[i, j] = snew / riverarea[i, j]
          else
            depth[i, j] = 0.0
          end
        end
      end
    end
  end

  # 使用扩散波方程计算新流量
  for j in (js+1):(je-1)
    for i in 2:(imax-1)
      if fd[i, j] > 0
        i1, j1 = flowdir(fd, i, j)

        if floodheight[i, j] > 0.05 || depth[i1, j1] > maxdepth[i1, j1] + 0.05
          # 洪水条件下使用扩散波
          if width[i, j] > 0
            aa = depth[i, j] * width[i, j] / (2.0 * depth[i, j] + width[i, j])
          else
            aa = depth[i, j]
          end

          waterelevij = topo[i, j] - maxdepth[i, j] + depth[i, j]
          waterelevi1j1 = topo[i1, j1] - maxdepth[i1, j1] + max(depth[i1, j1], 0.0)
          slopefor = (waterelevi1j1 - waterelevij) / (0.5 * (length[i, j] + length[i1, j1]))

          if bfd[i, j] > 0
            i2, j2 = flowdir(bfd, i, j)
            waterelevi2j2 = topo[i2, j2] - maxdepth[i2, j2] + max(depth[i2, j2], 0.0)
            slopeback = (waterelevij - waterelevi2j2) / (0.5 * (length[i2, j2] + length[i, j]))
            slopeinst = 0.5 * (slopefor + slopeback)
          else
            slopeinst = slopefor
          end

          # 扩散波方程
          qnew[i, j] = (q[i, j] - g0 * depth[i, j] * dtlr * slopeinst) /
                       (1.0 + g0 * dtlr * 0.03^2 * q[i, j] /
                              (aa^(4.0 / 3.0) * depth[i, j]))

          if width[i, j] == 0.0
            flowwidth = sqrt(area[i, j])
          else
            flowwidth = width[i, j]
          end

          qnew[i, j] *= flowwidth
        else
          # 非洪水条件下使用运动波
          aa = depth[i, j] * width[i, j] / (2.0 * depth[i, j] + width[i, j])
          speed = (aa^(2.0 / 3.0)) * sqrt(slope[i, j]) / 0.03
          speed = max(min(speed, length[i, j] / dtlr), 0.01)
          qnew[i, j] = speed * depth[i, j] * width[i, j]
        end
      else
        qnew[i, j] = 0.0
      end
    end
  end

  # 累积平均流量
  qmean .+= qnew .* dtlr

  return nothing
end


"""
    flooding!(洪水漫流参数...)

计算洪泛区的水位扩散
"""
function flooding!(
  imax::Int, jmax::Int, is::Int, ie::Int, js::Int, je::Int, Δt::T,
  fd::Matrix{Int}, bfd::Matrix{Int}, topo::M, area::M,
  riverwidth::M, riverlength::M, riverdepth::M,
  floodheight::M, delsfcwat::M
) where {T<:AbstractFloat,M<:Matrix{T}}
  ntsplit = 1  # 时间分割数
  dflood = zeros(size(floodheight))

  for n in 1:ntsplit
    for j in max(js + 1, 2):min(je - 1, jmax - 1)
      for i in max(is + 1, 2):min(ie - 1, imax - 1)
        if fd[i, j] == 0
          continue
        end

        if floodheight[i, j] > 0.05
          # 找到最低高程的邻居
          dhmax = 0.0
          ilow = i
          jlow = j

          for jj in (j-1):(j+1)
            for ii in (i-1):(i+1)
              if ii == i && jj == j
                continue
              end

              dh = floodheight[i, j] + topo[i, j] - (floodheight[ii, jj] + topo[ii, jj])
              if ii != i && jj != j
                dh /= sqrt(2.0)  # 对角线距离调整
              end

              if dh > dhmax
                ilow = ii
                jlow = jj
                dhmax = dh
              end
            end
          end

          # 向最低邻居漫流
          if dhmax > 0.0
            i1, j1 = flowdir(fd, i, j)

            dtotal = floodheight[i, j] + floodheight[ilow, jlow]
            dij = max(floodheight[i, j] - max(0.5 * (topo[ilow, jlow] - topo[i, j] + dtotal), 0.0), 0.0)

            # 如果是主河道方向，考虑河道流量
            if ilow == i1 && jlow == j1
              dij = max(dij - (riverwidth[i, j] * floodheight[i, j] * riverlength[i, j]) / area[i, j], 0.0)
            end

            # 考虑地表水变化
            if delsfcwat[i, j] < 0.0
              dij = max(min(dij, floodheight[i, j] + delsfcwat[i, j]), 0.0)
            end

            dflood[i, j] -= dij
            dflood[ilow, jlow] += dij * area[i, j] / area[ilow, jlow]
          end
        end
      end
    end

    # 更新地表水变化
    delsfcwat .+= dflood
    dflood .= 0.0
  end
  return nothing
end


"""
    moveqrf!(移动河流-地下水通量参数...)

将小河流网格的通量移动到下游
"""
function moveqrf!(
  imax::Int, js::Int, je::Int, fd::Matrix{Int},
  qrf::M, area::M, width::M
) where {T<:AbstractFloat,M<:Matrix{T}}
  qrfextra = zeros(size(qrf))

  for j in (js+1):(je-1)
    for i in 2:(imax-1)
      if fd[i, j] > 0
        if width[i, j] < 1.0
          iout, jout = flowdir(fd, i, j)
          qrfextra[iout, jout] += qrf[i, j] * area[i, j] / area[iout, jout]
          qrf[i, j] = 0.0
        end
      end
    end
  end
  # 应用额外通量
  qrf .+= qrfextra
  return nothing
end


"""
    apply_specific_diversions(特定分流应用...)

应用原始代码中硬编码的特定河流分流
"""
function apply_specific_diversions(i::Int, j::Int, q::M, dsnew::T) where {T<:AbstractFloat,M<:Matrix{T}}
  # Taquari河流分流
  (i == 4498 && j == 4535) && (dsnew += q[4499, 4534] / 4.0)
  (i == 4498 && j == 4534) && (dsnew -= q[4499, 4534] / 4.0)
  (i == 4464 && j == 4536) && (dsnew += q[4465, 4535] / 2.0)
  (i == 4465 && j == 4534) && (dsnew -= q[4465, 4535] / 2.0)
  (i == 4346 && j == 4560) && (dsnew += q[4346, 4561] / 3.0)
  (i == 4345 && j == 4561) && (dsnew -= q[4346, 4561] / 3.0)
  (i == 4444 && j == 4551) && (dsnew += q[4444, 4552] / 3.0)
  (i == 4443 && j == 4553) && (dsnew -= q[4444, 4552] / 3.0)
  (i == 4350 && j == 4497) && (dsnew += q[4352, 4496] / 3.0)
  (i == 4351 && j == 4496) && (dsnew -= q[4352, 4496] / 3.0)

  # Sao Lourenco河流分流
  (i == 4439 && j == 4772) && (dsnew += q[4440, 4773] / 2.0)
  (i == 4440 && j == 4772) && (dsnew -= q[4440, 4773] / 2.0)
  (i == 4400 && j == 4685) && (dsnew += q[4401, 4685] / 2.0)
  (i == 4400 && j == 4684) && (dsnew -= q[4401, 4685] / 2.0)
  (i == 4418 && j == 4688) && (dsnew += q[4418, 4689] / 5.0)
  (i == 4417 && j == 4689) && (dsnew -= q[4418, 4689] / 5.0)
  (i == 4367 && j == 4698) && (dsnew += q[4368, 4699] / 6.0)
  (i == 4368 && j == 4698) && (dsnew -= q[4368, 4699] / 6.0)
  (i == 4363 && j == 4667) && (dsnew += q[4364, 4668] / 6.0)
  (i == 4364 && j == 4667) && (dsnew -= q[4364, 4668] / 6.0)
  (i == 4475 && j == 4718) && (dsnew += q[4475, 4717] / 6.0)
  (i == 4474 && j == 4717) && (dsnew -= q[4475, 4717] / 6.0)
  return dsnew
end
