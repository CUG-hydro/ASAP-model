#import "@local/modern-cug-report:0.1.2": *
#show: doc => template(doc, footer: "CUG水文气象学2025", header: "下渗")

#beamer-block[王文焰2003，黄土中Green-Ampt 入渗模型的改进与验证]

$
  q = i = K(theta_s) [((H + 0) - (-0.5 L + S_m )) / (0.5 L)]
$

土壤基质势$S_m$为负，土壤越干燥、基质势越小。

$
  q d t  = (theta_s - theta_l) (d L)/2 + pi/4 (theta_s - theta_l) (d L) /2 \
  q d t = (pi + 4) / 8 (theta_s - theta_l) d L
$

第二项为椭圆面积公式。

$
  dv(L, t) = K(theta_s) [((H + 0) - (-0.5 L + S_m )) / (0.5 L)] 8 / (pi + 4) 1 / (theta_s - theta_l) \ 
  dv(L, t) = (16 K(theta_s) ) / {(pi + 4) (theta_s - theta_l)} { H + 0.5 L - S_m}/{L}
$

let $c = (16 K(theta_s) ) / {(pi + 4) (theta_s - theta_l)} $

$
  d L {L} / {H + 0.5L - S_m}= c thin d t \
  
  { [2(2H + L - 2S_m) - 4H - 4S_m ]} / {2H + L - 2S_m} d L = c thin d t \
  (2 - 4 (H + S_m) / (2H + L - 2S_m) ) d L = c thin d t
$

$
  2 L - 4 (H + S_m)  ln({2H + L - 2S_m} / {2H - 2S_m}) = c t
$

$
  t = {(4 + pi) (theta_s- theta_l)} / (16 K(theta_s)) [2 L - 4 (H + S_m)  ln({2H + L - 2S_m} / {2H - 2S_m})]
$

===

// 基于Julia语言的
// 考虑地下水垂向与侧向交互影响的简易水文模型研发与问题探讨

#v(2em)
// #pagebreak()
// 地下水在陆面水文循环过程中的重要性逐渐被广泛认知。但如何考虑地下水与土壤水在垂向交互，以及地下水侧向流动

// 土壤水与
// 
// 的影响，
// 
// 大尺度
// 
// 基于Julia语言的考虑侧向流动与垂向
// 水文过程的
