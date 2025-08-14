#import "@local/modern-cug-report:0.1.2": *
#import "@preview/mannot:0.3.0": *

// #import "@preview/modern-cug-report:0.1.0": *
// #import "../lib.typ": *
#counter(heading).update(2)
#show: doc => template(doc, footer: "CUG水文气象学2025", header: "土壤水运动")

#set par(spacing: 1.24em + 0.0em, leading: 1.24em)
#show math.equation: set text(size: 10.2pt)

*坐标系定义*：$z$和$Q$向下为负，向上为正。注意代码中的第1层为土壤底部；第N层为土壤表层。

$
  Q_(k+1/2) =
  underbrace(- K_(k+1/2) , "重力步长通量") -
  underbrace(D_(k+1/2)(theta_(k+1)^(n+1) - theta_k^(n+1)) /(Delta z_(k+1/2) ) , "毛管步长通量")
$

$
  Q_(k-1/2) =
  underbrace(- K_(k-1/2) , "重力步长通量") -
  underbrace(D_(k-1/2)(theta_(k)^(n+1) - theta_{k-1}^(n+1)) /(Delta z_(k-1/2) ) , "毛管步长通量")
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
  (Delta z_k)/(Delta t)( {theta_k^(n+1)} - theta_k^n ) = - Q_{k+1/2} + Q_{k-1/2} - S_k
$ <eq_4>

$
  (Delta z_k)/(Delta t)( mcurr(theta_k^(n+1)) - theta_k^n ) =
  (D_(k+1/2))/ ( Delta z_(k+1/2))( mnext(theta_(k+1)^(n+1)) - mcurr(theta_k^(n+1)) )
  - (D_(k-1/2))/(Delta z_(k-1/2))( mcurr(theta_k^(n+1)) - mprev(theta_(k-1)^(n+1)) )
  + (K_(k+1/2) - K_(k-1/2))
  - S_k
$ <eq_5>

其中，上标$n$代表时间，Q单位为$"mm"thin s^{-1}$。合并整理成$a theta_(k-1)^(n+1) + b theta_k^(n+1) + c theta_(k+1)^(n+1) = d$的形式：

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
  -underbrace( ((Delta z_k)/(Delta t) + (D_(k-1/2))/(Delta z_(k-1/2)) + (D_(k+1/2))/(Delta z_(k+1/2))), B_k) mcurr(theta_k^(n+1)) +
  underbrace((D_(k+1/2))/(Delta z_(k+1/2)), C_k) mnext(theta_(k+1)^(n+1))
  = - (Delta z_k)/(Delta t) theta_k^n  - K_(k+1/2) + K_(k-1/2) + S_k
$

#v(-0.6em)

- 第$3~N-1$层，上式与代码中的公式形式完全一致。
  #v(-0.6em)
  ```julia
  for k in max(iwtd, 3):(nzg-1)
    aa[k] = D_mid[k] / Δz₊ₕ[k]
    cc[k] = D_mid[k+1] / Δz₊ₕ[k+1]
    bb[k] = -(aa[k] + cc[k] + Δz[k] / dt)
    rr[k] = -θ[k] * Δz[k] / dt - K_mid[k+1] + K_mid[k] + transp[k] / dt
  end
  ```

- *上边界层$k = N$*

  对于第k=N层，已知$Q_k+1/2 = I$，将式#[@eq_2]带入式#[@eq_4]，重新整理可得：

  *(1) 若地下水未没过地表: *
  $
    (Delta z_k)/(Delta t)( mcurr({theta_k^(n+1)}) - theta_k^n ) = 
      - Q_{k+1/2} - K_(k-1/2) - D_(k-1/2) / (Delta z_(k-1/2) ) (mcurr(theta_(k)^(n+1)) - mprev(theta_{k-1}^(n+1))) 
      - S_k
  $

  $
    underbrace((D_(k-1/2))/(Delta z_(k-1/2)), A_k) mprev(theta_(k-1)^(n+1)) 
    -underbrace( ((Delta z_k)/(Delta t) + (D_(k-1/2))/(Delta z_(k-1/2)) ), B_k) mcurr(theta_k^(n+1)) +
    underbrace(0, C_k) mnext(theta_(k+1)^(n+1))
    = Q_{k+1/2} - (Delta z_k)/(Delta t) theta_k^n + K_(k-1/2) + S_k
  $

  ```julia
  aa[nzg] = D_mid[nzg] / Δz₊ₕ[nzg]
  cc[nzg] = 0.0
  bb[nzg] = -aa[nzg] - Δz[nzg] / dt
  rr[nzg] = Q[nzg+1] / dt - θ[nzg] * Δz[nzg] / dt + K_mid[nzg] + transp[nzg] / dt
  ```

  *(2) 若地下水没过地表: * $Q_{k-1/2} = K_{k-1/2}$
  $
    (Delta z_k)/(Delta t)( mcurr({theta_k^(n+1)}) - theta_k^n ) = - Q_{k+1/2} -
      K_(k-1/2) - S_k \ 
    - delta(z_k) / delta(t) mcurr({theta_k^(n+1)}) = Q_{k+1/2} - (Delta z_k)/(Delta t) theta_k^n + Q_{k-1/2} + S_k
  $

  ```julia
  aa[nzg] = 0.0
  cc[nzg] = 0.0
  bb[nzg] = -Δz[nzg] / dt
  rr[nzg] = Q[nzg+1] / dt - θ[nzg] * Δz[nzg] / dt + transp[nzg] / dt +
      min(K_mid[nzg] + D_mid[nzg] / Δz₊ₕ[nzg] * (θ[nzg] - θ[nzg-1]), 0.0)
  ```

- *下边界层$k = 1$*，已知$Q_{k-1/2}$
  
  *若为重力排水*，$Q_{k-1/2} = -K_{k-1/2}$

  $
  (Delta z_k)/(Delta t)( mcurr(theta_k^(n+1)) - theta_k^n ) = 
    (D_(k+1/2))/ ( Delta z_(k+1/2))( mnext(theta_(k+1)^(n+1)) - mcurr(theta_k^(n+1)) ) + K_{k+1/2} + Q_{k-1/2} - S_k
  $
  
  $
    underbrace(- (delta(z_k) / delta(t) +  (D_(k+1/2))/ ( Delta z_(k+1/2)) ), B_k)  mcurr(theta_k^(n+1)) + 
    underbrace((D_(k+1/2))/ ( Delta z_(k+1/2)), C_k)                                mnext(theta_(k+1)^(n+1)) = 
     - (Delta z_k)/(Delta t) theta_k^n  - K_(k+1/2) - Q_(k-1/2) + S_k
  $
  
  ```julia
  K_mid[1] = K * (θ[1] / smoisat)^(2.0 * soil_params.b + 3.0)
  aa[1] = 0.0
  cc[1] = D_mid[2] / Δz₊ₕ[2]
  bb[1] = -(cc[1] + Δz[1] / dt)
  rr[1] = -θ[1] * Δz[1] / dt - K_mid[2] + K_mid[1] + transp[1] / dt
  ```

  *若非重力排水*，

  #highlight[该情景下，模型是如何处理的？请补全文档说明]
