"""
    flowdir(流向参数...)

根据流向编码确定下游网格位置

流向编码：
- 1: 东 (→)
- 2: 东南 (↘)
- 4: 南 (↓)
- 8: 西南 (↙)
- 16: 西 (←)
- 32: 西北 (↖)
- 64: 北 (↑)
- 128: 东北 (↗)
"""
function flowdir(fd::Matrix{Int}, ii::Int, jj::Int)
  # 确定j方向
  j = if fd[ii, jj] in [2, 4, 8]
    jj - 1
  elseif fd[ii, jj] in [1, 16]
    jj
  elseif fd[ii, jj] in [32, 64, 128]
    jj + 1
  else
    0
  end

  # 确定i方向
  i = if fd[ii, jj] in [128, 1, 2]
    ii + 1
  elseif fd[ii, jj] in [4, 64]
    ii
  elseif fd[ii, jj] in [8, 16, 32]
    ii - 1
  else
    0
  end
  return i, j
end
