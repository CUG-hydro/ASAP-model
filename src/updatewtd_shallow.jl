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
    kwtd = jwt - 1

    if kwtd > 0  # 地下水位在已解析层中
      wtd_old = wtd
      soil = get_soil_params(soiltxt)

      # 检查上层是否接近饱和
      if kwtd > 1
        θ_sat = soil.θ_sat * cal_factor(z[kwtd - 1], fdepth)
        if wtd < z₋ₕ[kwtd] + 0.01 && θ[kwtd-1] < θ_sat
          flag = 1
        end
      end
      θ_sat = soil.θ_sat * cal_factor(z[kwtd], fdepth)

      if θ[kwtd] > θ_eq[kwtd] && flag == 0
        if θ[kwtd] == θ_sat  # 地下水位上升到上层
          wtd = z₋ₕ[jwt]
          rech = (wtd_old - wtd) * (θ_sat - θ_eq[kwtd])
          jwt = jwt + 1
          kwtd = kwtd + 1

          if kwtd <= nzg
            if θ[kwtd] > θ_eq[kwtd]
              wtd_old = wtd
              θ_sat = soil.θ_sat * cal_factor(z[kwtd], fdepth)
              wtd = min((θ[kwtd] * dz[kwtd] - θ_eq[kwtd] * z₋ₕ[jwt] + θ_sat * z₋ₕ[kwtd]) /
                        (θ_sat - θ_eq[kwtd]), z₋ₕ[jwt])
              rech = rech + (wtd_old - wtd) * (θ_sat - θ_eq[kwtd])
            end
          else
            break
          end
        else  # 地下水位在层内
          wtd = min((θ[kwtd] * dz[kwtd] - θ_eq[kwtd] * z₋ₕ[jwt] + θ_sat * z₋ₕ[kwtd]) /
                    (θ_sat - θ_eq[kwtd]), z₋ₕ[jwt])
          rech = (wtd_old - wtd) * (θ_sat - θ_eq[kwtd])
          break
        end

      else  # 地下水位下降到下层
        wtd = z₋ₕ[kwtd]
        rech = (wtd_old - wtd) * (θ_sat - θ_eq[kwtd])
        kwtd = kwtd - 1
        jwt = jwt - 1

        # 调整到下层
        if kwtd >= 1
          wtd_old = wtd
          θ_sat = soil.θ_sat * cal_factor(z[kwtd], fdepth)

          if θ[kwtd] > θ_eq[kwtd]
            wtd = min((θ[kwtd] * dz[kwtd] - θ_eq[kwtd] * z₋ₕ[jwt] + θ_sat * z₋ₕ[kwtd]) /
                      (θ_sat - θ_eq[kwtd]), z₋ₕ[jwt])
            rech = rech + (wtd_old - wtd) * (θ_sat - θ_eq[kwtd])
            break
          else
            wtd = z₋ₕ[kwtd]
            rech = rech + (wtd_old - wtd) * (θ_sat - θ_eq[kwtd])
          end
        else
          break
        end
      end
    else
      break
    end
  end
  wtd < z₋ₕ[1] && (println("地下水位问题: wtd=$wtd")) # 地下水水位过深
  return wtd, rech
end
