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
