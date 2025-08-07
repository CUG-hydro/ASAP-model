    WaterTable模块

实现地下水位动态模拟的综合功能，包括：

## 主要组件
- `Constants`: 物理常数定义
- `WaterTableCalculations`: 核心水位计算
- `GroundwaterRiverInteraction`: 河流-地下水相互作用
- `IsotopeTracing`: 同位素传输追踪

## 主要功能
1. **水位计算** (`wtable!`)
   - 地下水侧向流计算
   - 深层补给估算
   - 水位深度更新
   - 土壤湿度分布调整

2. **河流相互作用** (`gw2river!`)
   - 地下水向河流排水
   - 河流向地下水渗透
   - 河床连通性判断

3. **洪水路由** (`rivers_kw_flood!`, `rivers_dw_flood!`)
   - 运动波方程路由
   - 扩散波方程路由
   - 洪泛区漫流计算

4. **同位素追踪** (`lateral_isotope!`)
   - 侧向流同位素传输
   - 浓度梯度计算
   - 质量守恒追踪

## 变量命名约定
使用希腊字母表示物理量：
- `θ`: 体积含水量 (theta)
- `ψ`: 基质势 (psi) 
- `ρ`: 密度相关参数 (rho)
- `κ`: 导水率 (kappa)
- `Δt`: 时间步长 (Delta t)
- `α`, `β`, `γ`, `δ`, `λ`: 其他参数

## 使用示例
```julia
using ASAP.WaterTable

# 初始化参数
nzg = 10
slz = collect(range(-0.1, -5.0, length=nzg+1))
dz = -diff(slz)

# 创建网格数组
imax, jmax = 100, 100
wtd = zeros(imax, jmax)
smoi = zeros(nzg, imax, jmax)
# ... 其他数组初始化

# 主计算循环
Δt = 3600.0  # 1小时时间步长
wtable!(imax, jmax, 1, imax, 1, jmax, nzg,
        slz, dz, area, soiltxt, wtd, bottomflux,
        rech, qslat, fdepth, topo, landmask, Δt,
        smoi, smoieq, smoiwtd, qsprings)
```

## 算法说明
### 侧向流计算
使用达西定律的有限差分近似：
```
q = -K * ∇h * A
```
其中 K 为导水率，h 为水头，A 为截面积。

### 水位更新
基于质量守恒原理：
```
∂θ/∂t = ∇·(K∇h) + S
```
其中 S 为源汇项。

### 河流交换
考虑三种情况：
1. 地下水位 > 河流水位：排水
2. 地下水位 < 河流水位但连通：补给  
3. 地下水位 < 河床：渗透

## 注意事项
- 所有长度单位为米(m)，时间单位为秒(s)
- 流量单位为 m³/s，转换为mm时乘以1000
- 数组索引遵循Julia的1-based约定
- 边界条件需要在调用前正确设置
