"""
土壤参数模块
定义各种土壤类型的参数和常数
"""
module SoilParameters

export SoilType, get_soil_params, init_soil_param, khyd

# 常数定义
const NVTYP = 30  # 植被类型数量
const NSTYP = 13  # 土壤类型数量

# 土壤参数结构体
struct SoilType
    slmsts::Float64   # 饱和含水量
    soilcp::Float64   # 土壤容重
    slbs::Float64     # 土壤b参数
    slcons::Float64   # 饱和导水率
    slpots::Float64   # 饱和基质势
    slwilt::Float64   # 凋萎点含水量
    klatfactor::Float64  # 侧向流因子
end

# 土壤参数数据
const SLMSTS = [0.395, 0.410, 0.435, 0.485, 0.451, 0.420, 0.477, 0.476, 0.426, 0.492, 0.482, 0.863, 0.476]
const SOILCP = [0.050, 0.052, 0.092, 0.170, 0.125, 0.148, 0.195, 0.235, 0.202, 0.257, 0.268, 0.195, 0.235]
const SLBS = [4.05, 4.38, 4.9, 5.3, 5.39, 7.12, 7.75, 8.52, 10.4, 10.4, 11.4, 7.75, 8.52]
const SLCONS = [0.000176, 0.0001563, 0.00003467, 0.0000072, 0.00000695, 0.0000063, 
                0.0000017, 0.00000245, 0.000002167, 0.000001033, 0.000001283, 0.0000080, 0.000005787]
const SLPOTS = [-0.121, -0.090, -0.218, -0.786, -0.478, -0.299, -0.356, -0.630, -0.153, -0.490, -0.405, -0.356, -0.630]
const KLATFACTOR = [2.0, 3.0, 4.0, 10.0, 12.0, 14.0, 20.0, 24.0, 28.0, 40.0, 48.0, 48.0, 48.0]

"""
获取指定土壤类型的参数

# 参数
- `soil_type::Int`: 土壤类型索引 (1-13)

# 返回
- `SoilType`: 土壤参数结构体
"""
function get_soil_params(soil_type::Int)
    if soil_type < 1 || soil_type > NSTYP
        error("土壤类型索引必须在 1-$NSTYP 范围内")
    end
    
    return SoilType(
        SLMSTS[soil_type],
        SOILCP[soil_type],
        SLBS[soil_type],
        SLCONS[soil_type],
        SLPOTS[soil_type],
        0.0,  # slwilt 将在初始化时计算
        KLATFACTOR[soil_type]
    )
end

"""
初始化土壤参数，计算凋萎点

# 参数
- `nzg::Int`: 土壤层数

# 返回
- `Array{Float64,2}`: 田间持水量数组 (nzg × nstyp)
"""
function init_soil_param(nzg::Int)
    POTWILT_LOCAL = -153.0  # 凋萎点基质势
    
    fieldcp = zeros(Float64, nzg, NSTYP)
    slwilt = zeros(Float64, NSTYP)
    
    # 计算各土壤类型的凋萎点含水量
    for nsoil in 1:NSTYP
        slwilt[nsoil] = SLMSTS[nsoil] * (SLPOTS[nsoil] / POTWILT_LOCAL)^(1.0 / SLBS[nsoil])
    end
    
    return fieldcp, slwilt
end

"""
计算土壤导水率

# 参数
- `smoi::Float64`: 土壤含水量
- `nsoil::Int`: 土壤类型索引

# 返回
- `Float64`: 导水率
"""
function khyd(smoi::Float64, nsoil::Int)
    if nsoil < 1 || nsoil > NSTYP
        error("土壤类型索引必须在 1-$NSTYP 范围内")
    end
    
    return SLCONS[nsoil] * (smoi / SLMSTS[nsoil])^(2.0 * SLBS[nsoil] + 3.0)
end

end # module SoilParameters
