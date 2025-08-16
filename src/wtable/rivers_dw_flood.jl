
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

  for j in (js+1):(je-1), i in 2:(imax-1)
    if fd[i, j] != 0
      qext[i, j] = (qrf[i, j] + qs[i, j] + delsfcwat[i, j]) / Δt * area[i, j]
    end
  end

  q = copy(qnew)         # 保存当前流量
  qin = zeros(size(q))   # 计算入流

  for j in js:je, i in 1:imax
    if fd[i, j] > 0
      i1, j1 = flowdir(fd, i, j)
      if i1 > 1 && i1 < imax && j1 > js && j1 < je
        qin[i1, j1] += q[i, j]
      end
    end
  end

  # 更新河流深度和洪水高度 (类似KW方法)
  for j in (js+1):(je-1), i in 2:(imax-1)
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

  # 使用扩散波方程计算新流量
  for j in (js+1):(je-1), i in 2:(imax-1)
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

  qmean .+= qnew .* dtlr # 累积平均流量
  return nothing
end
