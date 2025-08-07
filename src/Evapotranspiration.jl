"""
蒸散发计算模块
包含多种蒸散发计算方法
"""
module Evapotranspiration

export potevap_priestly_taylor, potevap_penman_monteith, potevap_shutteworth_wallace


const CP = 1013.0  # J/kg/K
const VK = 0.41    # von Karman 常数
const RD = 287.0   # 干空气气体常数

# 植被叶阻力参数
const RL = [150.0, 150.0, 500.0, 500.0, 175.0, 240.0, 110.0, 100.0, 250.0, 150.0,
    80.0, 225.0, 225.0, 250.0, 180.0, 180.0, 240.0, 500.0, 240.0, 500.0,
    175.0, 250.0, 250.0, 175.0, 225.0, 150.0, 110.0, 180.0, 250.0, 250.0]

# 零位移高度
const ZDIS = [0.1, 0.1, 0.1, 15.0, 20.0, 15.0, 20.0, 0.2, 1.0, 0.1, 0.5, 0.1, 1.0, 1.0,
    20.0, 0.7, 0.7, 1.0, 10.2, 20.7, 9.2, 7.2, 6.5, 7.4, 3.6, 1.4, 0.2, 0.2, 0.2, 0.2, 1.1]


"""
Priestley-Taylor 法计算潜在蒸散发

# 参数
- `i::Int`, `j::Int`: 网格坐标
- `tempk::Float64`: 温度 (K)
- `rad::Float64`: 净辐射 (W/m²)
- `presshp::Float64`: 大气压力 (hPa)

# 返回
- `Float64`: 潜在蒸散发 (mm)
"""
function potevap_priestly_taylor(i::Int, j::Int, tempk::Float64, rad::Float64, presshp::Float64)
    const CP = 1013.0e-6  # 比热容

    tempc = tempk - 273.15  # 转换为摄氏度
    presskp = presshp * 0.1  # 转换为 kPa
    rad_mj = rad * 24.0 * 3600.0 * 1.0e-6  # 转换为 MJ/day/m²

    alpha = 1.26
    delta = 0.2 * (0.00738 * tempc + 0.8072)^7 - 0.000116
    lambda = 2.501 - 0.002361 * tempc
    gamma = (CP * presskp) / (0.622 * lambda)

    pet = alpha * rad_mj * delta / (delta + gamma)
    pet = pet / lambda

    return pet
end

"""
Penman-Monteith 法计算潜在蒸散发

# 参数
- `i::Int`, `j::Int`: 网格坐标  
- `tempk::Float64`: 温度 (K)
- `rad::Float64`: 净辐射 (W/m²)
- `rshort::Float64`: 短波辐射 (W/m²)
- `press::Float64`: 大气压力 (Pa)
- `qair::Float64`: 比湿
- `wind::Float64`: 风速 (m/s)
- `lai::Float64`: 叶面积指数
- `veg::Float64`: 植被类型
- `hveg::Float64`: 植被高度 (m)

# 返回
- `Float64`: 潜在蒸散发 (mm/3h)
"""
function potevap_penman_monteith(i::Int, j::Int, tempk::Float64, rad::Float64, rshort::Float64,
    press::Float64, qair::Float64, wind::Float64, lai::Float64,
    veg::Float64, hveg::Float64)

    tempc = tempk - 273.15
    pressesat = 610.8 * exp(17.27 * tempc / (tempc + 237.3))  # Pa
    pressvap = qair * press / (0.622 + qair)  # Pa
    delta = 4098.0 * pressesat / (tempc + 237.3)^2  # Pa/K
    vpd = pressesat - pressvap
    lambda = (2.501 - 0.002361 * tempc) * 1.0e6  # J/kg
    gamma = (CP * press) / (0.622 * lambda)
    dens = press / (RD * tempk * (1.0 + 0.608 * qair))

    # 空气动力学阻力计算
    hdisp = 0.7 * hveg
    zm = max(10.0, hveg)
    z0m = 0.1 * hveg
    zh = max(2.0, hveg)
    z0h = 0.1 * z0m

    ra = log((zm - hdisp) / z0m) * log((zh - hdisp) / z0h) / (VK^2 * wind)

    # 体表阻力计算
    frad = min(1.0, (0.004 * rshort + 0.05) / (0.81 * (1.0 + 0.004 * rshort)))
    fswp = 1.0

    g_d = hdisp > 2.0 ? 0.0003 : 0.0
    fvpd = exp(-g_d * vpd)

    slai = 0.5 * lai
    veg_int = round(Int, veg)
    veg_int = max(1, min(veg_int, length(RL)))

    if slai * frad * fswp * fvpd == 0.0
        rs = 5000.0
    else
        rs = min(RL[veg_int] / (slai * frad * fswp * fvpd), 5000.0)
    end

    # 计算蒸散发
    pet = (delta * rad + dens * CP * vpd / ra) / (delta + gamma * (1.0 + rs / ra))
    pet = 3.0 * 3600.0 * pet / lambda

    return pet
end


"""
Shuttleworth-Wallace 双源法计算蒸散发组分

# 参数
- `i::Int`, `j::Int`: 网格坐标
- `deltat::Float64`: 时间步长 (s)
- `tempk::Float64`: 温度 (K)
- `rad::Float64`: 净辐射 (W/m²)
- `rshort::Float64`: 短波辐射 (W/m²)
- `press::Float64`: 大气压力 (Pa)
- `qair::Float64`: 比湿
- `wind::Float64`: 风速 (m/s)
- `lai::Float64`: 叶面积指数
- `veg::Float64`: 植被类型
- `hhveg::Float64`: 植被高度 (m)
- `floodflag::Int`: 洪水标志

# 返回
- `Tuple`: (delta, gamma, lambda, ra_a, ra_c, rs_c, R_a, R_s, pet_s, pet_c, pet_w, pet_i)
"""
function potevap_shutteworth_wallace(i::Int, j::Int, deltat::Float64, tempk::Float64, rad::Float64,
    rshort::Float64, press::Float64, qair::Float64, wind::Float64,
    lai::Float64, veg::Float64, hhveg::Float64, floodflag::Int)
    const CP = 1013.0  # J/kg/K
    const VK = 0.41    # von Karman 常数  
    const RD = 287.0   # 干空气气体常数

    # 植被叶阻力参数
    const RL = [150.0, 150.0, 500.0, 500.0, 175.0, 240.0, 110.0, 100.0, 250.0, 150.0,
        80.0, 225.0, 225.0, 250.0, 180.0, 180.0, 240.0, 500.0, 240.0, 500.0,
        175.0, 250.0, 250.0, 175.0, 225.0, 150.0, 110.0, 180.0, 250.0, 250.0]

    # 生物物理参数 [地表粗糙度, 最大叶宽度]
    const BIOPARMS = [
        0.001 0.0;    # 0  海洋
        0.001 0.0;    # 1  湖泊、河流
        0.001 0.0;    # 2  冰川
        0.02 0.001;   # 3  常绿针叶林
        0.02 0.001;   # 4  落叶针叶林
        0.02 0.08;    # 5  落叶阔叶林
        0.02 0.05;    # 6  常绿阔叶林
        0.01 0.01;    # 7  短草
        0.01 0.01;    # 8  高草
        0.001 0.01;   # 9  沙漠
        0.01 0.01;    # 10 半沙漠
        0.01 0.01;    # 11 苔原
        0.02 0.01;    # 12 常绿灌木
        0.02 0.01;    # 13 落叶灌木
        0.02 0.04;    # 14 混交林
        0.005 0.01;   # 15 农田
        0.005 0.01;   # 16 灌溉农田
        0.01 0.01;    # 17 沼泽
        0.01 0.001;   # 18 常绿针叶林
        0.02 0.05;    # 19 常绿阔叶林
        0.02 0.001;   # 20 落叶针叶林
        0.02 0.08;    # 21 落叶阔叶林
        0.01 0.01;    # 22 混合覆盖
        0.02 0.04;    # 23 林地
        0.02 0.01;    # 24 草木地
        0.02 0.01;    # 25 密闭灌丛
        0.02 0.01;    # 26 开放灌丛
        0.01 0.01;    # 27 草地
        0.005 0.01;   # 28 农田
        0.001 0.01;   # 29 裸地
        0.02 0.0      # 30 城市
    ]

    tempc = tempk - 273.15
    pressesat = 610.8 * exp(17.27 * tempc / (tempc + 237.3))  # Pa
    pressvap = qair * press / (0.622 + qair)  # Pa
    delta = 4098.0 * pressesat / (tempc + 237.3)^2  # Pa/K
    vpd = pressesat - pressvap
    lambda = (2.501 - 0.002361 * tempc) * 1.0e6  # J/kg
    gamma = (CP * press) / (0.622 * lambda)
    dens = press / (RD * tempk * (1.0 + 0.608 * qair))

    # 确保植被高度不为零
    hveg = max(hhveg, 0.1)

    if round(Int, veg) <= 1
        # 水体或裸地
        pet_w = (delta * rad + gamma * 6.43 * (1.0 + 0.536 * wind) * vpd / (24.0 * 3600.0)) / (delta + gamma)
        pet_w = max(deltat * pet_w / lambda, 0.0)

        return (delta, gamma, lambda, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, pet_w, 0.0)
    else
        pet_w = 0.0

        # 地面净辐射
        Rn_s = rad * exp(-0.5 * lai)

        # 阻力计算
        za = hveg + 2.0  # 参考高度

        # 冠层粗糙度
        if hveg <= 1.0
            z0c = 0.13 * hveg
        elseif hveg > 1.0 && hveg < 10.0
            z0c = 0.139 * hveg - 0.009 * hveg^2
        else
            z0c = 0.05 * hveg
        end

        # 平均阻力系数
        c_d = hveg == 0.0 ? 1.4e-3 : (-1.0 + exp(0.909 - 3.03 * z0c / hveg))^4 / 4.0

        # 零位移高度
        d0 = lai >= 4.0 ? max(hveg - z0c / 0.3, 0.0) : 1.1 * hveg * log(1.0 + (c_d * lai)^0.25)

        # 地表粗糙度
        veg_int = round(Int, veg)
        veg_int = max(1, min(veg_int, size(BIOPARMS, 1)))
        z0g = floodflag == 0 ? BIOPARMS[veg_int, 1] : BIOPARMS[1, 1]

        # 冠层粗糙度
        z0 = min(0.3 * (hveg - d0), z0g + 0.3 * hveg * (c_d * lai)^0.5)
        z0 = max(z0, z0g)

        # 摩擦速度
        ustar = VK * wind / log(10.0 / z0)

        # 涡扩散系数
        K_h = VK * ustar * (hveg - d0)

        # 植被衰减常数
        if hveg <= 1.0
            n = 2.5
        elseif hveg > 1.0 && hveg < 10.0
            n = 2.306 + 0.194 * hveg
        else
            n = 4.25
        end

        # 首选参数
        Z0 = 0.13 * hveg
        dp = 0.63 * hveg

        # 空气动力学阻力
        ra_a = log((za - d0) / (hveg - d0)) / (VK * ustar) +
               hveg * (exp(n * (1.0 - (Z0 + dp) / hveg)) - 1.0) / (n * K_h)

        # 土壤到冠层阻力
        ra_s = hveg * exp(n) * (exp(-n * z0g / hveg) - exp(-n * (Z0 + dp) / hveg)) / (n * K_h)

        # 冠层边界层阻力
        uc = ustar * log((hveg - d0) / z0) / VK

        # 叶宽度
        wleaf = veg_int in [4, 5, 13, 20, 21] ?
                BIOPARMS[veg_int, 2] * (1.0 - exp(-0.6 * lai)) :
                BIOPARMS[veg_int, 2]

        rb = 100.0 * (wleaf / uc)^0.5 / ((1.0 - exp(-n / 2.0)) * n)

        ra_c = lai > 0.1 ? rb * 0.5 / lai : 0.0

        # 冠层气孔阻力
        frad = min(1.0, (0.004 * rshort + 0.05) / (0.81 * (1.0 + 0.004 * rshort)))
        fswp = 1.0
        g_d = d0 > 2.0 ? 0.0003 : 0.0
        fvpd = exp(-g_d * vpd)
        slai = 0.5 * lai

        if slai * frad * fswp * fvpd == 0.0
            rs_c = 5000.0
        else
            rs_c = min(RL[veg_int] / (slai * frad * fswp * fvpd), 5000.0)
        end

        # 组合阻力
        R_a = (delta + gamma) * ra_a
        R_c = (delta + gamma) * ra_c + gamma * rs_c
        R_s = (delta + gamma) * ra_s

        # 蒸散发计算
        pet_c = (delta * rad + (dens * CP * vpd - delta * ra_c * Rn_s) / (ra_a + ra_c))
        pet_s = (delta * rad + (dens * CP * vpd - delta * ra_s * (rad - Rn_s)) / (ra_a + ra_s))

        # 截留蒸发（气孔阻力为0）
        C_c = lai < 0.001 ? 0.0 : 1.0 / (1.0 + R_a * R_c / (R_s * (R_c + R_a)))
        pet_i = C_c * (delta * rad + (dens * CP * vpd - delta * ra_c * Rn_s) / (ra_a + ra_c)) / (delta + gamma)
        pet_i = max(deltat * pet_i / lambda, 0.0)

        return (delta, gamma, lambda, ra_a, ra_c, rs_c, R_a, R_s, pet_s, pet_c, pet_w, pet_i)
    end
end

end # module Evapotranspiration
