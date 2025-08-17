using DataFrames

# 预定义的层厚度（40层）
const DZ2 = [
  0.1, 0.1, 0.1, 0.1, 0.1, 0.2, 0.2, 0.2, 0.2, 0.2,
  0.3, 0.3, 0.3, 0.3, 0.4, 0.4, 0.4, 0.5, 0.5, 0.6,
  0.7, 0.7, 0.8, 0.9, 1.0, 1.0, 1.2, 1.2, 1.5, 1.5,
  2.0, 2.0, 3.0, 6.0, 11.0, 20.0, 50.0, 100.0, 250.0, 540.0]

export initializesoildepth_clm, initializesoildepth


"""
CLM方法初始化土壤层深度

# 参数
- `nzg::Int`: 土壤层数

# 返回
- `slz`: [N, 0], 土壤上边界，向下为负
- `dz` : [N, 1], 土壤深度，正值

# Note
该算法存在问题，slz计算到了N+1层，但dz没有N+1层
"""
function initializesoildepth_clm(nzg::Int)
  slz = zeros(Float64, nzg + 1)
  dz = zeros(Float64, nzg)
  vctr4 = zeros(Float64, nzg) # 土壤的下边界

  # 计算节点深度（指数分布）
  for k in 1:nzg
    vctr4[k] = 0.025 * (exp(0.5 * (k - 0.5)) - 1.0)
  end

  # 计算层厚度
  for k in 2:(nzg-1)
    dz[k] = 0.5 * (vctr4[k+1] - vctr4[k-1])
  end
  dz[1] = 0.5 * (vctr4[1] + vctr4[2]) # vctr4[0] = -vctr4[1]
  dz[nzg] = vctr4[nzg] - vctr4[nzg-1]

  # 计算层中心深度, 与下文计算的对象不同
  slz[1] = 0.5vctr4[1]
  for k in 2:nzg
    slz[k] = 0.5 * (vctr4[k-1] + vctr4[k]) # 中心
  end
  slz[nzg+1] = vctr4[nzg] + 0.5 * dz[nzg] # 计算到了N+1层

  # 反转深度（从地表向下为负值）
  (; slz=reverse(slz), dz=reverse(dz))
end


"""
固定层厚方法初始化土壤层深度

# 参数
- `nzg::Int`: 土壤层数

# 返回  
- `slz`: [N, 0], 土壤下边界
- `dz` : [N, 1], 土壤层厚度
"""
function initializesoildepth(nzg::Int)
  slz = zeros(Float64, nzg + 1)
  dz = zeros(Float64, nzg)

  if nzg > length(DZ2)
    error("土壤层数不能超过 $(length(DZ2))")
  end
  slz[nzg+1] = 0.0  # 地表深度

  # 从地表向下计算各层深度
  for k in nzg:-1:1
    dz[k] = DZ2[nzg-k+1]
    slz[k] = slz[k+1] - dz[k]
  end
  return slz, dz
end
