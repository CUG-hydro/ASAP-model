"""
主要的根系深度计算函数

使用泛型类型参数提供灵活性：
- T: 数值类型 (通常为 Float64 或 Float32)
- V: 向量类型 Vector{T} 
- M: 矩阵类型 Matrix{T}
- A3: 三维数组类型 Array{T,3}

# 参数
- `freedrain::Int`: 自由排水标志
- `is::Int`, `ie::Int`, `js::Int`, `je::Int`: 网格范围
- `nzg::Int`: 土壤层数
- `slz::V`: 土壤层深度 (m) [nzg+1]
- `dz::V`: 土壤层厚度 (m) [nzg]
- `Δt::Float64`: 时间步长 (s)
- `landmask::Matrix{Int}`: 陆地掩码 [is:ie, js:je]
- `veg::M`: 植被类型 [is:ie, js:je]
- `hveg::M`: 植被高度 (m) [is:ie, js:je]
- `soiltxt::Array{Int,3}`: 土壤类型 [2, is:ie, js:je]
- `wind::M`: 风速 (m/s) [is:ie, js:je]
- `temp::M`: 温度 (K) [is:ie, js:je]
- `qair::M`: 比湿 [is:ie, js:je]
- `press::M`: 大气压力 (Pa) [is:ie, js:je]
- `netrad::M`: 净辐射 (W/m²) [is:ie, js:je]
- `rshort::M`: 短波辐射 (W/m²) [is:ie, js:je]
- `lai::M`: 叶面积指数 [is:ie, js:je]
- `precip::M`: 降水 (mm) [is:ie, js:je]
- `qsrun::M`: 地表径流 (m) [is:ie, js:je]
- `smoi::A3`: 土壤含水量 [nzg, is:ie, js:je]
- `smoieq::A3`: 平衡含水量 [nzg, is:ie, js:je]
- `smoiwtd::M`: 地下水含水量 [is:ie, js:je]
- `wtd::M`: 地下水位深度 (m) [is:ie, js:je]
- `waterdeficit::M`: 水分不足量 (mm) [is:ie, js:je]
- `watext::A3`: 水分提取量 (mm) [nzg, is:ie, js:je]
- `watextdeep::M`: 深层水分提取 (mm) [is:ie, js:je]
- `rech::M`: 补给量 (mm) [is:ie, js:je]
- `deeprech::M`: 深层补给 (mm) [is:ie, js:je]
- `et_s::M`: 土壤蒸发 (mm) [is:ie, js:je]
- `et_i::M`: 截留蒸发 (mm) [is:ie, js:je]
- `et_c::M`: 冠层蒸腾 (mm) [is:ie, js:je]
- `intercepstore::M`: 截留存储 (mm) [is:ie, js:je]
- `ppacum::M`: 累计降水 (mm) [is:ie, js:je]
- `pppendepth::M`: 降水渗透深度 (m) [is:ie, js:je]
- `pppendepthold::Matrix{Int8}`: 渗透深度记录 [is:ie, js:je]
- `qlat::M`: 侧向流 (m) [is:ie, js:je]
- `qlatsum::M`: 累计侧向流 (m) [is:ie, js:je]
- `qsprings::M`: 泉水流 (m) [is:ie, js:je]
- `inactivedays::Array{Int,3}`: 非活跃天数 [0:nzg+1, is:ie, js:je]
- `maxinactivedays::Int`: 最大非活跃天数
- `fieldcp::M`: 田间持水量 [nzg, nstyp]
- `fdepth::M`: 根系深度因子 (m) [is:ie, js:je]
- `steps::Float64`: 时间步数
- `floodheight::M`: 洪水高度 (m) [is:ie, js:je]
- `qrf::M`: 河流反馈流 (m) [is:ie, js:je]
- `delsfcwat::M`: 地表水变化 (m) [is:ie, js:je]
- `icefactor::Array{Int8,3}`: 冰冻因子 [is:ie, js:je, 26:40]
- `wtdflux::M`: 地下水通量 (mm) [is:ie, js:je]
- `et_s_daily::M`: 日土壤蒸发 (mm) [is:ie, js:je]
- `et_c_daily::M`: 日冠层蒸腾 (mm) [is:ie, js:je]
- `transptop::M`: 表层蒸腾 (mm) [is:ie, js:je]
- `infilk::Matrix{Int8}`: 入渗层指标 [is:ie, js:je]
- `infilflux::A3`: 入渗通量 (mm) [nzg, is:ie, js:je]
- `infilfluxday::A3`: 日入渗通量 (mm) [nzg, is:ie, js:je]
- `infilcounter::Array{Int16,3}`: 入渗计数器 [nzg, is:ie, js:je]
- `hour::Int`: 小时
- `o18::A3`: 氧18同位素 [nzg, is:ie, js:je]
- `o18ratiopp::M`: 降水氧18比值 [is:ie, js:je]
- `tempsfc::M`: 地表温度 (K) [is:ie, js:je]
- `qlato18::M`: 侧向流氧18 [is:ie, js:je]
- `transpo18::M`: 蒸腾氧18 [is:ie, js:je]
- `upflux::A3`: 上升通量 (mm) [nzg, is:ie, js:je]

# 返回
- 更新所有输入的状态变量

# 示例
```julia
# 使用 Float64 类型
slz_f64 = Vector{Float64}([...])
smoi_f64 = Array{Float64,3}([...])
rootdepth_main(freedrain, is, ie, js, je, nzg, slz_f64, ..., smoi_f64, ...)

# 使用 Float32 类型提高性能
slz_f32 = Vector{Float32}([...])
smoi_f32 = Array{Float32,3}([...])
rootdepth_main(freedrain, is, ie, js, je, nzg, slz_f32, ..., smoi_f32, ...)
```
"""
function rootdepth_main(
  freedrain::Int, is::Int, ie::Int, js::Int, je::Int, nzg::Int,
  z₋ₕ::V, dz::V, Δt::Float64,
  landmask::Matrix{Int}, veg::M, hveg::M,
  soiltxt::Array{Int,3}, wind::M, temp::M,
  qair::M, press::M, netrad::M,
  rshort::M, lai::M, precip::M,
  qsrun::M, θ::A3, θ_eq::A3,
  θ_wtd::M, wtd::M, waterdeficit::M,
  watext::A3,
  watextdeep::M, rech::M,
  deeprech::M,
  et_s::M, et_i::M, et_c::M,
  intercepstore::M, ppacum::M,
  pInfiltDepth::M, pInfiltDepthK_old::Matrix{Int8}, qlat::M,
  qlatsum::M, qsprings::M, inactivedays::Array{Int,3},
  maxinactivedays::Int, fdepth::M,
  steps::Float64, floodheight::M, qrf::M,
  delsfcwat::M, icefactor::Array{Int8,3}, wtdflux::M,
  et_s_daily::M, et_c_daily::M, transptop::M,
  infilk::Matrix{Int8}, infilflux::A3, infilfluxday::A3,
  infilcounter::Array{Int16,3}, hour::Int,
  o18::A3, o18ratiopp::M, tempsfc::M,
  qlato18::M, transpo18::M, upflux::A3) where {
  T<:Real,V<:Vector{T},M<:Matrix{T},A3<:Array{T,3}}


  # 初始化冰冻因子数组
  icefac = zeros(Int8, nzg)

  # 主循环
  for j in (js+1):(je-1)
    for i in (is+1):(ie-1)
      landmask[i, j] == 0 && continue
      
      # 复制冰冻因子
      icefac[26:40] .= icefactor[i, j, 26:40]

      # 检查是否有洪水
      floodflag = floodheight[i, j] > 0.05 ? 1 : 0

      # 计算 Shuttleworth-Wallace 蒸散发
      Δ, γ, λ, ra_a, ra_c, rs_c, R_a, R_s, petfactor_s, petfactor_c, petstep_w, petstep_i =
        potevap_shutteworth_wallace(Δt, temp[i, j], netrad[i, j], rshort[i, j],
          press[i, j], qair[i, j], wind[i, j], lai[i, j],
          veg[i, j], hveg[i, j], floodflag)

      # 累计土壤蒸发
      et_s[i, j] += petstep_w
      if floodflag == 1 && round(Int, veg[i, j]) <= 1
        delsfcwat[i, j] -= petstep_w * 1.0e-3
      end

      # 如果是水体或裸地，跳过植被相关计算
      round(Int, veg[i, j]) <= 1 && continue
      
      # 植被截留计算
      ppdrip, etstep_i, new_intercepstore = interception(precip[i, j], lai[i, j],
        intercepstore[i, j], petstep_i)

      # 更新截留蒸发和存储
      et_i[i, j] += etstep_i
      intercepstore[i, j] = new_intercepstore

      # 分步计算
      ppdrip_step = ppdrip / steps  # mm
      floodstep = floodheight[i, j] / steps  # m
      qlatstep = qlat[i, j] / steps  # m
      qrfstep = qrf[i, j] / steps  # m
      qlato18step = qlato18[i, j] / steps  # m

      # 初始化通量数组
      flux = zeros(Float64, nzg + 1)
      qlatflux = zeros(Float64, nzg + 2)

      wtd_old = wtd[i, j]

      # 时间子循环
      for itime in 1:round(Int, steps)
        # 水分提取计算
        pet_s, pet_c, watdef, dθ, dsmoideep = extraction(nzg, z₋ₕ, dz, Δt / steps,
          soiltxt[1, i, j], 
          wtd[i, j], θ[:, i, j], θ_wtd[i, j], 
          Δ, γ, λ, lai[i, j],
          ra_a, ra_c, rs_c, R_a, R_s, petfactor_s, petfactor_c, 
          inactivedays[:, i, j], maxinactivedays, hveg[i, j], fdepth[i, j], icefac)

        # 更新蒸腾和水分不足量
        et_c[i, j] += pet_c - watdef * 1.0e3
        waterdeficit[i, j] += watdef * 1.0e3
        watext[:, i, j] .+= dθ .* 1.0e3
        transptop[i, j] += dθ[nzg] * 1.0e3
        et_c_daily[i, j] += pet_c - watdef * 1.0e3

        # 土壤水流计算
        # updated_o18, transpo18step
        et_s_step, runoff, rechstep, flux_step, qrfcorrect, updated_smoi =
          soilfluxes(nzg, Δt / steps, z₋ₕ, dz, soiltxt[1, i, j],
            θ_wtd[i, j],
            dθ, dsmoideep, θ[:, i, j], wtd[i, j], ppdrip_step, pet_s,
            fdepth[i, j], qlatstep, qrfstep, floodstep, icefac
            ; freedrain
          )
        # θ_eq[:, i, j], o18[:, i, j], o18ratiopp[i, j], tempsfc[i, j],
        # qlato18step, transpo18step

        # 更新状态变量
        delsfcwat[i, j] -= max(floodstep - runoff, 0.0)  # m
        qsrun[i, j] += max(runoff - floodstep, 0.0)  # m
        rech[i, j] += rechstep * 1.0e3
        et_s[i, j] += et_s_step
        transpo18[i, j] += transpo18step * 1.0e3
        et_s_daily[i, j] += et_s_step
        ppacum[i, j] += ppdrip_step
        qrf[i, j] += qrfcorrect

        # 更新土壤含水量和同位素
        θ[:, i, j] .= updated_smoi
        o18[:, i, j] .= updated_o18

        # 更新浅层地下水位
        wtd[i, j], rech_additional = updateshallowwtd(i, j, nzg, freedrain, z₋ₕ, dz,
          soiltxt[1, i, j], θ_eq[:, i, j], θ_wtd[i, j], θ[:, i, j], wtd[i, j], fdepth[i, j])

        rech[i, j] += rech_additional * 1.0e3
        qlatsum[i, j] += qlatstep

        # 累计通量
        flux .+= flux_step
      end

      # 保存入渗信息
      for k in nzg:-1:1
        if flux[k+1] < 0.0
          infilflux[k, i, j] += flux[k+1] * 1.0e3
        else
          upflux[k, i, j] += flux[k+1] * 1.0e3
        end
        infilfluxday[k, i, j] += flux[k+1] * 1.0e3
      end

      # 入渗计数器
      if hour == 0
        for k in nzg:-1:1
          if infilfluxday[k, i, j] < -0.01
            infilcounter[k, i, j] += 1
          end
          infilfluxday[k, i, j] = 0.0
        end
      end

      # 确定入渗深度
      _kInfilt = nzg + 1
      _pInfiltDepth = 0.0
      flux[nzg+1] = -1.0

      for k in nzg:-1:0
        if k <= nzg - 2
          if pInfiltDepthK_old[i, j] >= k + 3 # 下渗深度，不可能一下增加三层
            break
          end
        end
        
        # 0.012 [cm/h] to [m/s]
        if flux[k+1] < -0.333e-5  # 时间步长从3h改为1h后相应调整阈值
          if k == 0
            if -flux[1] > -qlatflux[1] && _pInfiltDepth > z₋ₕ[1]
              _pInfiltDepth = z₋ₕ[1]
              _kInfilt = 1
            end
          elseif -flux[k+1] + flux[k] > -qlatflux[k] + dθ[k] && _pInfiltDepth > z₋ₕ[k+1]
            _pInfiltDepth = z₋ₕ[k+1]
            _kInfilt = k + 1
          end
        end
      end

      pInfiltDepthK_old[i, j] = _kInfilt
      pInfiltDepth[i, j] > _pInfiltDepth && (pInfiltDepth[i, j] = _pInfiltDepth)
      
      if z₋ₕ[max(_kInfilt - 1, 1)] <= wtd_old
        wtdflux[i, j] -= flux[_kInfilt] * 1.0e3
      end
      infilk[i, j] > _kInfilt && (infilk[i, j] = _kInfilt)
    end
  end
  return nothing  # 所有变量都是原地更新
end
