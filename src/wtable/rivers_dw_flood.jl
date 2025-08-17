const g0 = 9.81  # 重力加速度 (m/s²)

"""
    rivers_dw_flood!(扩散波河流洪水路由参数...)

使用扩散波方程进行河流洪水路由计算
"""
function rivers_dw_flood!(
  imax::Int, js::Int, je::Int, Δt::T, δt::T,
  fd::Matrix{Int}, bfd::Matrix{Int}, qnew::M,
  qs::M, qrf::M, delsfcwat::M,
  slope::M, depth::M, width::M,
  length::M, maxdepth::M, area::M,
  riverarea::M, floodarea::M, riverchannel::M,
  qmean::M, floodheight::M, topo::M
) where {T<:AbstractFloat,M<:Matrix{T}}

  qext = zeros(size(qnew)) # 计算外部输入流量

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
             (dsnew + qext[i, j]) * δt

      if snew >= riverchannel[i, j]
        floodheight[i, j] = (snew - riverchannel[i, j]) / max(area[i, j], riverarea[i, j])
        depth[i, j] = floodheight[i, j] + maxdepth[i, j]
      else
        floodheight[i, j] = 0.0
        depth[i, j] = riverarea[i, j] > 0.0 ? snew / riverarea[i, j] : 0.0
      end
    end
  end

  # 使用扩散波方程计算新流量
  for j in (js+1):(je-1), i in 2:(imax-1)
    qnew[i, j] = 0.0
    fd[i, j] <= 0 && continue

    i1, j1 = flowdir(fd, i, j)
    A = depth[i, j] * width[i, j]

    if floodheight[i, j] > 0.05 || depth[i1, j1] > maxdepth[i1, j1] + 0.05
      # 洪水条件下使用扩散波
      R = width[i, j] > 0 ? A / (2.0 * depth[i, j] + width[i, j]) : depth[i, j]

      z0 = topo[i, j] - maxdepth[i, j] + depth[i, j]
      z1 = topo[i1, j1] - maxdepth[i1, j1] + max(depth[i1, j1], 0.0)
      slope_for = (z1 - z0) / (0.5 * (length[i, j] + length[i1, j1]))

      if bfd[i, j] > 0
        i2, j2 = flowdir(bfd, i, j)
        z2 = topo[i2, j2] - maxdepth[i2, j2] + max(depth[i2, j2], 0.0)
        slope_back = (z0 - z2) / (0.5 * (length[i2, j2] + length[i, j]))
        slope_inst = 0.5 * (slope_for + slope_back)
      else
        slope_inst = slope_for
      end

      # 扩散波方程
      qnew[i, j] = (q[i, j] - g0 * depth[i, j] * δt * slope_inst) /
                   (1.0 + g0 * δt * 0.03^2 * q[i, j] /
                          (R^(4.0 / 3.0) * depth[i, j]))

      flowwidth = width[i, j] == 0.0 ? sqrt(area[i, j]) : width[i, j]
      qnew[i, j] *= flowwidth
    else
      # 非洪水条件下使用运动波
      R = A / (2.0 * depth[i, j] + width[i, j])
      speed = (R^(2.0 / 3.0)) * sqrt(slope[i, j]) / 0.03
      speed = clamp(speed, 0.01, length[i, j] / δt)
      qnew[i, j] = speed * A
    end
  end

  qmean .+= qnew .* δt # 累积平均流量
  return nothing
end
