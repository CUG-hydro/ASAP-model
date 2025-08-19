"""
    rivers_kw_flood!(运动波河流洪水路由参数...)

使用运动波方程进行河流洪水路由计算

# Arguments
- `qrf`: runoff–fast/overland, 地表径流
- `qs` : runoff–slow/subsurface, 地下径流
- `delsfcwat`: Δ surface water storage to channel
"""
function rivers_kw_flood!(
  imax::Int, jmax::Int, is::Int, ie::Int, js::Int, je::Int,
  Δt::T, δt::T, 
  fd::Matrix{Int}, bfd::Matrix{Int},
  qnew::M, qs::M, qrf::M,
  delsfcwat::M, depth::M,
  riverarea::M, floodarea::M, floodheight::M,
  width::M, length::M, maxdepth::M, slope::M, area::M, topo::M, riverchannel::M, # 静态变量
  qmean::M,
) where {T<:AbstractFloat,M<:Matrix{T}}
  
  qext = zeros(size(qnew)) # 计算外部输入流量
  for j in max(js + 1, 2):min(je - 1, jmax - 1)
    for i in max(is + 1, 2):min(ie - 1, imax - 1)
      if fd[i, j] != 0
        qext[i, j] = (qrf[i, j] + qs[i, j] + delsfcwat[i, j]) / Δt * area[i, j]
      end
    end
  end

  q = deepcopy(qnew)   # 保存当前流量
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
      fd[i, j] == 0 && continue

      dsnew = qin[i, j] - q[i, j]  # 计算流入该网格的总流量变化, [m3 s-1]
      dsnew = apply_specific_diversions(i, j, q, dsnew) # 添加特定河流分流(原始Fortran中的硬编码分流)

      # 新的河流存储量
      snew = depth[i, j] * riverarea[i, j] + floodheight[i, j] * floodarea[i, j] +
             (dsnew + qext[i, j]) * δt

      # 在河道和洪泛区之间重新分配水量
      if snew >= riverchannel[i, j]
        floodheight[i, j] = (snew - riverchannel[i, j]) / max(area[i, j], riverarea[i, j])
        depth[i, j] = floodheight[i, j] + maxdepth[i, j]
      else
        floodheight[i, j] = 0.0
        depth[i, j] = riverarea[i, j] > 0.0 ? snew / riverarea[i, j] : 0.0
      end
    end
  end


  # 计算曼宁公式计算新的流量
  for j in (js+1):(je-1), i in (is+1):(ie-1)
    qnew[i, j] = 0.0
    A = width[i, j] * depth[i, j] # 过水断面
    A < 1e-9 || fd[i, j] <= 0 && continue

    if floodheight[i, j] > 0.05
      i1, j1 = flowdir(fd, i, j)
      z0 = topo[i, j] - maxdepth[i, j] + depth[i, j]                  # 当前网格水面高程
      z1 = topo[i1, j1] - maxdepth[i1, j1] + max(depth[i1, j1], 0.0)  # 下游网格水面高程

      slopefor = (z0 - z1) / (0.5 * (length[i, j] + length[i1, j1]))

      if bfd[i, j] > 0
        i2, j2 = flowdir(bfd, i, j)                                   # 上游网格水面高程
        z2 = topo[i2, j2] - maxdepth[i2, j2] + max(depth[i2, j2], 0.0)
        slope_back = (z2 - z0) / (0.5 * (length[i2, j2] + length[i, j]))
        slope_inst = 0.5 * (slopefor + slope_back)
      else
        slope_inst = slopefor
      end
      slope_inst = 0.25 * slope_inst + 0.75 * slope[i, j]
      slope_inst < 0.0 && (slope_inst = slope[i, j])
    else
      slope_inst = slope[i, j]
    end
    R = A / (2.0 * depth[i, j] + width[i, j]) # 水力半径
    speed = (R^(2.0 / 3.0)) * sqrt(slope_inst) / 0.03
    speed = clamp(speed, 0.01, length[i, j] / δt)
    qnew[i, j] = speed * A # 计算新流量
  end

  qmean .+= qnew .* δt # 累积平均流量
  return nothing
end
