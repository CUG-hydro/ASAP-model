#import "@local/modern-cug-report:0.1.2": *
// #import "@preview/mannot:0.3.0": *

#show: doc => template(doc, footer: "CUG水文气象学2025", header: "河道汇流")

基于圣维南方程组动量守恒公式，忽略侧向动量项$q$与风$F_f$、局部损失项$F_e$：

#box-red[
  *运动波*：陡峭山区，惯性项与附加比降都可忽略 
  #v(-0.4em)
  *扩散波*：河漫滩，惯性项不能忽略
]

$
  pdv(Q, t) + pdv(, x)(Q^2/A) + g A (pdv(h, x) + S_f) = 0
$

其中 $h = y + z$, $y$是水深、$z$是河床高度，河床比降$pdv(z, x) = -S_0$。

$
  pdv(h, x) = pdv(y, x) - S_0, 
$


== 1 等价形式

=== 1.1 q的形式

$q$为单宽流量（$q = Q/B$, $A = B y$），因此$q = V y$

$
  pdv(q B, t) + pdv(, x)({q^2 B^2}/ (B y)) + g B y (pdv(h, x) + S_f) = 0
$ <eq_q>

#mitex(`$$
\underbrace{\frac{\partial q}{\partial t}}_{\text{局部惯性}}
+\underbrace{\frac{\partial}{\partial x}\!\left(\frac{q^2}{y}\right)}_{\text{对流惯性}}
+\underbrace{g y \frac{\partial h}{\partial x}}_{\text{压强/水面坡度}}
+\underbrace{g y S_f}_{\text{摩阻}} = 0
$$`)

=== 1.2 V的形式

$
  pdv(Q, t) = pdv(A V, t) = A pdv(V, t) + V pdv(A, t)
$

#mitex(`$$
\frac{\partial}{\partial x}\!\left(\frac{Q^2}{A}\right)
=\frac{\partial}{\partial x}(A V^2)
=\underbrace{V^2 \frac{\partial A}{\partial x}}_{\text{产出自 }A}
+\underbrace{2AV \frac{\partial V}{\partial x}}_{\text{产出自 }V}\tag{2}
$$`)

结合连续性方程，

#mitex(`$$
\frac{\partial A}{\partial t} + \frac{\partial Q}{\partial x} = 0, \frac{\partial A}{\partial t}
=-\Bigl(V \frac{\partial A}{\partial x}+A \frac{\partial V}{\partial x}\Bigr)\tag{3}
$$`)

综上可解得：
#mitex(`$$
\frac{\partial V}{\partial t}
+V \frac{\partial V}{\partial x}
+g\!\left(\frac{\partial h}{\partial x}+S_f\right)=0
$$`)

Chow 2003, Eq. 9.1.37。

== 3 扩散波

#h(2em)
式#[@eq_q]忽略对流惯性项，

$
  pdv(q, t) + g y pdv(h, x) + g y S_f = 0 \ 
  
  pdv(q, t) + g y (pdv(h, x) + S_f) = 0
  // pdv(q, t) + g y (pdv(y, x) - S_0) + g y S_f = 0 \ 
  
  // pdv(q, t) + g y pdv(y, x) - g y (S_0 - S_f) = 0 \ 
$

根据曼宁公式$ V = sqrt(S_f) / n R^(2/3), S_f = n^2 / R^{4/3} q^2 / y^2$

将

$
  (q^{n+1} - q^n) / delta(t) = - g y^n (S_h^{n} + S_f^{n+1})
$
其中，上标$n$，$n+1$代表时刻。

#mitex(`$$
S_f^{ n+1}\ \approx\ \gamma^n q^{n+1},
\qquad
\gamma^n
:= \frac{n^2 q^n}{R^{4/3} (y^n)^2}\quad\text{（冻结系数线性化）}
$$`)

#mitex(`$$
\boxed{\
q^{n+1}
= \frac{ q^n - g y^n \Delta t S_h^n }
       { 1 + g \Delta t y^n \gamma^n }
= \frac{ q^n - g y^n \Delta t S_h^n }
       { 1 + g \Delta t \dfrac{n^2 q^n}{R^{4/3} y^n} }}
\tag{U}
$$`)

```julia
qnew = ( q - g0*depth*δt*slope_inst ) /
       ( 1 + g0*δt*0.03^2*q / ( R^(4/3)*depth ) )
```
