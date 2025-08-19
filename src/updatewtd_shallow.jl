# z: 负值
function cal_factor(z::T, fdepth::T) where {T<:Real}
  clamp(exp((z + 1.5) / fdepth), 0.1, 1.0)
end


"""
地下水位更新模块
计算浅层地下水位的变化

# 参数
- `nzg::Int`: 土壤层数
- `freedrain::Int`: 自由排水标志
- `slz::Vector{Float64}`: 土壤层深度 (m) [nzg+1]
- `dz::Vector{Float64}`: 土壤层厚度 (m) [nzg]
- `soiltxt::Int`: 土壤类型
- `smoieq::Vector{Float64}`: 平衡含水量 [nzg]
- `smoiwtd::Float64`: 地下水含水量
- `smoi::Vector{Float64}`: 土壤含水量 [nzg]
- `wtd::Float64`: 地下水位深度 (m)
- `fdepth::Float64`: 根系深度因子 (m)

# 返回
- `Tuple{Float64, Float64}`: (更新的地下水位, 补给量)
"""
function updatewtd_shallow(nzg::Int, freedrain::Int, z₋ₕ::Vector{Float64},
  dz::Vector{Float64}, soiltxt::Int, θ_eq::Vector{Float64},
  θ_wtd::Float64, θ::Vector{Float64}, wtd::Float64, fdepth::Float64)

  rech = 0.0
  z = 0.5 .* (z₋ₕ[1:nzg] .+ z₋ₕ[2:nzg+1])  # 层中心深度
  jwt = find_jwt(wtd, z₋ₕ)                 # 地下水水位所在的下一层

  while true
    flag = 0
    kwt = jwt - 1
    kwt <= 0 && break # 地下水位在已解析层中

    wtd_old = wtd
    soil = get_soil_params(soiltxt)

    # 检查上层是否接近饱和
    if kwt > 1
      θ_sat = soil.θ_sat * cal_factor(z[kwt-1], fdepth)
      if wtd < z₋ₕ[kwt] + 0.01 && θ[kwt-1] < θ_sat
        flag = 1 # 这种情况是错误的？需要进行修正
      end
    end
    θ_sat = soil.θ_sat * cal_factor(z[kwt], fdepth)

    if θ[kwt] > θ_eq[kwt] && flag == 0 # 水位上升
      if θ[kwt] == θ_sat  # 上升一层
        wtd = z₋ₕ[jwt]
        rech = (wtd_old - wtd) * (θ_sat - θ_eq[kwt]) # release, 正释放了多少水, 若水位上升，则吸收水
        jwt = jwt + 1
        kwt = kwt + 1

        kwt > nzg && break
        # 已知这一层的土壤含水量，推求zwt到哪了
        if θ[kwt] > θ_eq[kwt]
          wtd_old = wtd
          θ_sat = soil.θ_sat * cal_factor(z[kwt], fdepth)
          wtd = min( (θ[kwt] * dz[kwt] - θ_eq[kwt] * z₋ₕ[jwt] + θ_sat * z₋ₕ[kwt]) /
                    (θ_sat - θ_eq[kwt]), z₋ₕ[jwt] ) # 更新之后的地下水水位
          rech = rech + (wtd_old - wtd) * (θ_sat - θ_eq[kwt])
        end
      else  # 地下水位在层内
        wtd = min((θ[kwt] * dz[kwt] - θ_eq[kwt] * z₋ₕ[jwt] + θ_sat * z₋ₕ[kwt]) /
                  (θ_sat - θ_eq[kwt]), z₋ₕ[jwt])
        rech = (wtd_old - wtd) * (θ_sat - θ_eq[kwt])
        break
      end

    else  # 地下水位下降到下层
      wtd = z₋ₕ[kwt]
      rech = (wtd_old - wtd) * (θ_sat - θ_eq[kwt])
      kwt = kwt - 1
      jwt = jwt - 1

      # 调整到下层
      if kwt >= 1
        wtd_old = wtd
        θ_sat = soil.θ_sat * cal_factor(z[kwt], fdepth)

        if θ[kwt] > θ_eq[kwt]
          wtd = min((θ[kwt] * dz[kwt] - θ_eq[kwt] * z₋ₕ[jwt] + θ_sat * z₋ₕ[kwt]) /
                    (θ_sat - θ_eq[kwt]), z₋ₕ[jwt])
          rech = rech + (wtd_old - wtd) * (θ_sat - θ_eq[kwt])
          break
        else
          wtd = z₋ₕ[kwt]
          rech = rech + (wtd_old - wtd) * (θ_sat - θ_eq[kwt])
        end
      else
        break
      end
    end

  end
  wtd < z₋ₕ[1] && (println("地下水位问题: wtd=$wtd")) # 地下水水位过深
  return wtd, rech
end
