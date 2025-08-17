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
        fd[i, j] == 0 && continue

        if floodheight[i, j] > 0.05
          # 找到最低高程的邻居
          dhmax = 0.0
          ilow = i
          jlow = j

          for jj in (j-1):(j+1), ii in (i-1):(i+1)
            ii == i && jj == j && continue

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

  for j in (js+1):(je-1), i in 2:(imax-1)
    fd[i, j] == 0 && continue
    if width[i, j] < 1.0
      iout, jout = flowdir(fd, i, j)
      qrfextra[iout, jout] += qrf[i, j] * area[i, j] / area[iout, jout]
      qrf[i, j] = 0.0
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
