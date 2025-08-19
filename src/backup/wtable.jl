"""
    wtable!(水位参数...)

主要的水位计算子程序，处理侧向流、深层补给和水位更新

参数：
- imax, jmax : 网格尺寸
- is, ie, js, je: 计算范围
- nzg        : 土壤层数
- slz        : 土壤层深度
- dz         : 土壤层厚度
- area       : 网格面积
- soiltxt    : 土壤类型
- wtd        : 水位深度
- bottomflux : 底部通量
- rech       : 补给量
- qslat      : 侧向流
- fdepth     : 排水深度
- topo       : 地形高程
- landmask   : 陆地掩膜
- Δt         : 时间步长
- smoi       : 土壤湿度
- smoieq     : 平衡土壤湿度
- smoiwtd    : 水位处土壤湿度
- qsprings   : 泉水流量
"""
function wtable!(
  imax::Int, jmax::Int, is::Int, ie::Int, js::Int, je::Int,
  nzg::Int, z₋ₕ::V, dz::V, area::M,
  soiltxt::Array{Int,3}, wtd::M, bottomflux::M,
  rech::M, qslat::M, fdepth::M,
  topo::M, landmask::Matrix{Int}, Δt::T,
  θ::Array{T,3}, θ_eq::Array{T,3}, θ_wtd::M,
  qsprings::M
) where {T<:AbstractFloat,V<:Vector{T},M<:Matrix{T}}

  # 计算侧向流导水率
  κlat = zeros(size(wtd))

  for j in js:je
    for i in 1:imax
      nsoil = soiltxt[1, i, j]
      soil = get_soil_params(nsoil)
      κlat[i, j] = soil.Ksat * soil.K_latfactor
    end
  end

  # 计算侧向流
  qlat = zeros(size(wtd))
  lateral_flow!(imax, jmax, is, ie, js, je, wtd, qlat, fdepth, topo,
    landmask, Δt, area, κlat)

  # 累积侧向流 (转换为mm)
  qslat .+= qlat .* 1e3

  # 计算深层补给
  deeprech = zeros(size(wtd))

  for j in (js+1):(je-1)
    for i in 1:imax
      if landmask[i, j] == 1
        if wtd[i, j] < z₋ₕ[1] - dz[1]
          # 计算排水导水率
          nsoil = soiltxt[1, i, j]
          soil = get_soil_params(nsoil)

          θ_mid = 0.5 * (θ_wtd[i, j] + soil.θ_sat)
          K = soil.Ksat * (θ_mid / soil.θ_sat)^(2.0 * soil.b + 3.0)

          # 计算湿度势
          ψ = soil.ψsat * (soil.θ_sat / θ_wtd[i, j])^soil.b

          # 计算通量(=补给), 
          deeprech[i, j] = Δt * K * ((soil.ψsat - ψ) / (z₋ₕ[1] - wtd[i, j]) - 1.0)

          # 更新水位处土壤湿度
          _θ_wtd = θ_wtd[i, j] + (deeprech[i, j] - bottomflux[i, j]) / (z₋ₕ[1] - wtd[i, j])

          if _θ_wtd < soil.θ_cp
            deeprech[i, j] += (soil.θ_cp - _θ_wtd) * (z₋ₕ[1] - wtd[i, j])
            _θ_wtd = soil.θ_cp
          end

          if _θ_wtd > soil.θ_sat
            deeprech[i, j] -= (soil.θ_sat - _θ_wtd) * (z₋ₕ[1] - wtd[i, j])
            _θ_wtd = soil.θ_sat
          end

          θ_wtd[i, j] = _θ_wtd
          rech[i, j] += deeprech[i, j] * 1e3
        end
      end
    end
  end

  # 重置底部通量
  bottomflux .= 0.0

  # 更新水位和土壤湿度
  for j in (js+1):(je-1)
    for i in 1:imax
      if landmask[i, j] == 1
        # 地下水平衡总量
        totwater = qlat[i, j] - deeprech[i, j]

        qspring = 0.0
        updatewtd!(nzg, z₋ₕ, dz, wtd[i, j], qspring, totwater,
          view(θ, :, i, j), view(θ_eq, :, i, j),
          soiltxt[:, i, j], θ_wtd[i, j])

        qsprings[i, j] += qspring * 1e3
      end
    end
  end

  return nothing
end
