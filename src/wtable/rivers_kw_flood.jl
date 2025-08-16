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

  q = copy(qnew)       # 保存当前流量
  qin = zeros(size(q)) # 计算入流

  for j in js:je, i in is:ie
    if fd[i, j] > 0
      i1, j1 = flowdir(fd, i, j)
      if i1 > is && i1 < ie && j1 > js && j1 < je
        qin[i1, j1] += q[i, j]
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
  for j in (js+1):(je-1), i in (is+1):(ie-1)
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

  qmean .+= qnew .* dtlr # 累积平均流量
  return nothing
end
