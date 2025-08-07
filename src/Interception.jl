"""
截留模块
计算植被截留和截留蒸发
"""
module Interception

export interception

"""
计算植被截留和截留蒸发

# 参数
- `minpprate::Float64`: 最小降水率阈值 (mm/timestep)
- `precip::Float64`: 降水量 (mm)
- `lai::Float64`: 叶面积指数
- `intercepstore::Float64`: 截留水量 (mm)
- `pet_i::Float64`: 截留蒸发潜力 (mm)

# 返回
- `Tuple{Float64, Float64, Float64}`: (穿透降水量, 截留蒸发量, 更新的截留水量)
"""
function interception(minpprate::Float64, precip::Float64, lai::Float64, 
                     intercepstore::Float64, pet_i::Float64)
    
    # 最大截留量，与叶面积指数成正比
    intercepmax = 0.2 * lai
    
    # 截留容量不足量
    deficit = intercepmax - intercepstore
    
    ppdrip = 0.0      # 穿透降水
    et_i = 0.0        # 截留蒸发
    new_intercepstore = intercepstore  # 新的截留水量
    
    if precip > deficit
        # 降水量超过截留容量不足量
        if precip < minpprate
            # 降水量小于阈值，有截留蒸发损失
            et_i = min(intercepmax, pet_i)
        else
            # 降水量大于阈值，无截留蒸发损失
            et_i = 0.0
        end
        
        new_intercepstore = intercepmax - et_i
        ppdrip = precip - deficit
        
    else
        # 降水量不超过截留容量不足量
        if precip < minpprate
            # 降水量小于阈值，有截留蒸发损失
            et_i = min(intercepstore + precip, pet_i)
        else
            # 降水量大于阈值，无截留蒸发损失  
            et_i = 0.0
        end
        
        new_intercepstore = intercepstore + precip - et_i
        ppdrip = 0.0
    end
    
    return ppdrip, et_i, new_intercepstore
end

end # module Interception
