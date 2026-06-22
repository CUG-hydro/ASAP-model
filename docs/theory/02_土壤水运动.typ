// docs/theory/02_土壤水运动.typ
// ============================
// 章：土壤水运动（1D Richards 方程）
// 状态：骨架（待摄取；现有 docs/土壤水运动_ASAP.typ 保持讲稿风格不动）
// 来源：src/SoilFluxes.jl::soilfluxes
// Fortran 参考：fortran/soilfluxes.f90::SOILFLUXES
#import "@preview/modern-cug-report:0.1.3": *
#show: doc => template(doc, size: 12pt, footer: "CUG水文气象学2026", header: "")
#import "../preamble.typ": *

// 公式编号（被 @eq_2 / @eq_4 等引用的公式需要编号）
// supplement: none → 引用显示为「式 (3)」而非「式 Equation 3」
#set math.equation(numbering: "(1)", supplement: none)

= 1 土壤水运动

== 1.1 物理过程

土壤水运动由 Richards 方程控制，变量为体积含水量 $theta$，参数包括
Campbell 导水率 $K(theta)$ 与水力扩散度 $D(theta)$。

$ pdv(theta, t) = pdv(, z) [K(theta) + D(theta) pdv(theta, z)] - S $ <eq:richards_1d>

*数值离散与求解*

- Crank-Nicolson 时间离散（L+ L- 加权）
- 中心差分空间离散
- Thomas 三对角求解（`src/SoilFluxes.jl` 中 `tridag!` 实现）
- 顶部入渗能力截断 + 底部自由/受限排水分支
- 含水量限幅 $[theta_"cp", theta_"sat"]$

*代码实现*

- 主求解器：`src/SoilFluxes.jl:L42-L342`（`soilfluxes` 全函数）
- Thomas 算法：`src/SoilFluxes.jl:LXXX-LXXX`（`tridag!`，待定位行号）
- 公式对照：`fortran/soilfluxes.f90:LXXX-LXXX`
- 测试断言：`test/test_soil_fluxes.jl`


*局限与已知问题*

- 同位素追踪段（`src/SoilFluxes.jl:L347-L351`）按 CLAUDE.md §6.1 暂缓启用。
- 现有讲稿式文档 `docs/土壤水运动_ASAP.typ` 提供更详细的公式推导，本章
  后续章节应与之交叉引用。

== 1.3 参数与单位

#figure(
  caption: [土壤水运动模型的参数与单位],
  kind: table,
  three-line-table(
    fill: (_, y) => if calc.odd(y) { rgb("EAF2F5") },
  )[
    | 符号            | 含义                       | 单位   |
    | --------------- | -------------------------- | ------ |
    | $theta$         | 体积含水量                 | m³/m³  |
    | $K(theta)$      | Campbell 导水率            | m/s    |
    | $D(theta)$      | 水力扩散度 $K pdv(psi, theta)$ | m²/s |
    | $psi$           | 基质势                     | m      |
    | $Q$             | 通量（向下为负）           | mm/s   |
    | $theta_"sat"$   | 饱和含水量                 | m³/m³  |
    | $theta_"cp"$    | 残余/凋萎含水量下限        | m³/m³  |
    | $theta_"fc"$    | 田间持水量                 | m³/m³  |
  ],
) <tbl:params>


== 1.4 核心理论

#let code = false

定义$z$和$Q$向下为负，向上为正。第1层为土壤底部；第N层为土壤表层。
// 注意代码中的

$
  Q_(k+1/2) =
  underbrace(- K_(k+1/2), "重力步长通量") -
  underbrace(D_(k+1/2)(theta_(k+1)^(n+1) - theta_k^(n+1)) /(Delta z_(k+1/2) ), "毛管步长通量")
$

$
  Q_(k-1/2) =
  underbrace(- K_(k-1/2), "重力步长通量") -
  underbrace(D_(k-1/2)(theta_(k)^(n+1) - theta_(k-1)^(n+1)) /(Delta z_(k-1/2) ), "毛管步长通量")
$ <eq_2>

#let mprev = markhl.with(color: red)
#let mcurr = markhl.with(color: yellow)
#let mnext = markhl.with(color: blue)
// #let theta-prev = $markhl(theta_{k-1}^{n+1}, color: #purple)$
// #let theta-curr = $markul(theta_{k}^{n+1})$
// #let theta-next = $markhl(theta_{k+1}^{n+1}, color: #blue)$
//
// $ pdv(theta, t) = - pdv(Q, z) -S =  pdv(, z)[D(theta) pdv(theta, z) - K(theta)] - S $

$
  (Delta z_k)/(Delta t)( theta_k^(n+1) - theta_k^n ) = - Q_(k+1/2) + Q_(k-1/2) - S_k
$ <eq_4>

$
  (Delta z_k)/(Delta t)( mcurr(theta_k^(n+1)) - theta_k^n ) =
  (D_(k+1/2))/ ( Delta z_(k+1/2))( mnext(theta_(k+1)^(n+1)) - mcurr(theta_k^(n+1)) )
  - (D_(k-1/2))/(Delta z_(k-1/2))( mcurr(theta_k^(n+1)) - mprev(theta_(k-1)^(n+1)) )
  + (K_(k+1/2) - K_(k-1/2))
  - S_k
$ <eq_5>

其中，上标$n$代表时间，Q单位为$"mm"thin s^(-1)$。合并整理成$a theta_(k-1)^(n+1) + b theta_k^(n+1) + c theta_(k+1)^(n+1) = d$的形式：

// $
//   underbrace((-D_(k-1/2))/(Delta z_(k-1/2)), A_k) mprev(theta_(k-1)^(n+1)) +
//   underbrace( ((Delta z_k)/(Delta t) + (D_(k-1/2))/(Delta z_(k-1/2)) + (D_(k+1/2))/(Delta z_(k+1/2))), B_k) mcurr(theta_k^(n+1)) +
//   underbrace((-D_(k+1/2))/(Delta z_(k+1/2)), C_k) mnext(theta_(k+1)^(n+1))
//   = (Delta z_k)/(Delta t) theta_k^n
//   + (K_(k+1/2) - K_(k-1/2)) - S_k
// $
// 等号两边同时乘-1即可得到代码中的公式形式：

$
  underbrace((D_(k-1/2))/(Delta z_(k-1/2)), A_k) mprev(theta_(k-1)^(n+1))
  -underbrace(((Delta z_k)/(Delta t) + (D_(k-1/2))/(Delta z_(k-1/2)) + (D_(k+1/2))/(Delta z_(k+1/2))), B_k) mcurr(theta_k^(n+1)) +
  underbrace((D_(k+1/2))/(Delta z_(k+1/2)), C_k) mnext(theta_(k+1)^(n+1))
  = - (Delta z_k)/(Delta t) theta_k^n - K_(k+1/2) + K_(k-1/2) + S_k
$

#v(-0.6em)

- 第$3~N-1$层，上式与代码中的公式形式完全一致。
  #v(-0.6em)
  #if code {
    ```julia
    for k in max(iwtd, 3):(nzg-1)
      aa[k] = D_mid[k] / Δz₊ₕ[k]
      cc[k] = D_mid[k+1] / Δz₊ₕ[k+1]
      bb[k] = -(aa[k] + cc[k] + Δz[k] / dt)
      rr[k] = -θ[k] * Δz[k] / dt - K_mid[k+1] + K_mid[k] + transp[k] / dt
    end
    ```
  }

- *上边界层$k = N$*

  对于第k=N层，已知$Q_(k+1/2) = I$，将式 @eq_2 带入式 @eq_4，重新整理可得：

  *(1) 若地下水未没过地表: *
  $
    (Delta z_k)/(Delta t)( mcurr(theta_k^(n+1)) - theta_k^n ) =
    - Q_(k+1/2) - K_(k-1/2) - D_(k-1/2) / (Delta z_(k-1/2) ) (mcurr(theta_(k)^(n+1)) - mprev(theta_(k-1)^(n+1)))
    - S_k
  $

  $
    underbrace((D_(k-1/2))/(Delta z_(k-1/2)), A_k) mprev(theta_(k-1)^(n+1))
    -underbrace(((Delta z_k)/(Delta t) + (D_(k-1/2))/(Delta z_(k-1/2)) ), B_k) mcurr(theta_k^(n+1)) +
    underbrace(0, C_k) mnext(theta_(k+1)^(n+1))
    = Q_(k+1/2) - (Delta z_k)/(Delta t) theta_k^n + K_(k-1/2) + S_k
  $

  #if code {
    ```julia
    aa[nzg] = D_mid[nzg] / Δz₊ₕ[nzg]
    cc[nzg] = 0.0
    bb[nzg] = -aa[nzg] - Δz[nzg] / dt
    rr[nzg] = Q[nzg+1] / dt - θ[nzg] * Δz[nzg] / dt + K_mid[nzg] + transp[nzg] / dt
    ```
  }

  *(2) 若地下水没过地表: * $Q_(k-1/2) = K_(k-1/2)$
  $
    (Delta z_k)/(Delta t)( mcurr(theta_k^(n+1)) - theta_k^n ) = - Q_(k+1/2) -
    K_(k-1/2) - S_k \
    - (Delta z_k)/(Delta t) mcurr(theta_k^(n+1)) = Q_(k+1/2) - (Delta z_k)/(Delta t) theta_k^n + Q_(k-1/2) + S_k
  $

  #if code {
    ```julia
    aa[nzg] = 0.0
    cc[nzg] = 0.0
    bb[nzg] = -Δz[nzg] / dt
    rr[nzg] = Q[nzg+1] / dt - θ[nzg] * Δz[nzg] / dt + transp[nzg] / dt +
        min(K_mid[nzg] + D_mid[nzg] / Δz₊ₕ[nzg] * (θ[nzg] - θ[nzg-1]), 0.0)
    ```
  }

- *下边界层$k = 1$*，已知$Q_(k-1/2)$

  *若为重力排水*，$Q_(k-1/2) = -K_(k-1/2)$
  $
    (Delta z_k)/(Delta t)( mcurr(theta_k^(n+1)) - theta_k^n ) =
    (D_(k+1/2))/ ( Delta z_(k+1/2))( mnext(theta_(k+1)^(n+1)) - mcurr(theta_k^(n+1)) ) + K_(k+1/2) + Q_(k-1/2) - S_k
  $

  $
    underbrace(- ((Delta z_k)/(Delta t) + (D_(k+1/2))/ ( Delta z_(k+1/2)) ), B_k) mcurr(theta_k^(n+1)) +
    underbrace((D_(k+1/2))/ ( Delta z_(k+1/2)), C_k) mnext(theta_(k+1)^(n+1)) =
    - (Delta z_k)/(Delta t) theta_k^n - K_(k+1/2) - Q_(k-1/2) + S_k
  $

  #if code {
    ```julia
    K_mid[1] = K * (θ[1] / smoisat)^(2.0 * soil_params.b + 3.0)
    aa[1] = 0.0
    cc[1] = D_mid[2] / Δz₊ₕ[2]
    bb[1] = -(cc[1] + Δz[1] / dt)
    rr[1] = -θ[1] * Δz[1] / dt - K_mid[2] + K_mid[1] + transp[1] / dt
    ```
  }

*若非重力排水*，地下水淹没部分Q=0。
// #highlight[该情景下，模型是如何处理的？请补全文档说明]
#include "../Figures/Figure_Soil_ASAP.typ"

- 第$1 ~ "jwt"-3$层（图中的1\~2层），$Q_"in" = 0, Q_"out" = 0$；

- 第$"jwt" - 2$层（图中第k=3层），$Q_(k-1) = 0$、$Q_k$正常计算，限制$Q_k$的方向为向下排水；

- 第$"jwt" - 1$层（图中第k=4层），$Q_(k-1)$、$Q_k$均正常计算，限制$Q_(k-1)$的方向为向下排水。
