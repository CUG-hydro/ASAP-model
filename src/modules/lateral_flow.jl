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
  # 流动角度因子, 除2是为了应对(κcell[i-1, j+1] + κcell[i, j])
  fangle = sqrt(tan(4π / 32.0)) / (2.0 * sqrt(2.0)) 

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

        qlat[i, j] = fangle * q * Δt / area[i, j] # 转化为深度
      end
    end
  end
  return nothing
end
