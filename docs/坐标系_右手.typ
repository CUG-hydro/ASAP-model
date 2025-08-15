#import "@preview/cetz:0.4.1"
#cetz.canvas({
  import cetz.draw: *

  // 基础倾斜视角（便于观察3D）：α绕x，β绕y
  let α = 50deg
  let β = 0deg
  // 题意：沿 z 轴顺时针 30°，按右手系取负角
  let γ = 40deg

  // 预计算三角函数
  let cα = calc.cos(α); let sα = calc.sin(α)
  let cβ = calc.cos(β); let sβ = calc.sin(β)
  let cγ = calc.cos(γ); let sγ = calc.sin(γ)

  // 组合旋转 R = Rz(γ) · Ry(β) · Rx(α)；投到2D：M = P·R
  let m11 =  cβ * cγ
  let m12 = -cβ * sγ
  let m13 =  sβ
  let m21 =  sα * sβ * cγ + cα * sγ
  let m22 = -sα * sβ * sγ + cα * cγ
  let m23 = -sα * cβ

  // 应用 4×4 变换矩阵（前两行给 x', y'；第三行不用；第四行齐次常数）
  set-transform((
    (m11, m12, m13, 0),
    (m21, m22, m23, 0),
    (0,   0,   0,   0),
    (0,   0,   0,   1),
  ))

  // 画坐标轴与网格
  let L = 6.0
  line((0,0,0), (L,0,0), mark: (end: ">"), name: "x")
  line((0,0,0), (0,L,0), mark: (end: ">"), name: "y")
  line((0,0,0), (0,0,L), mark: (end: ">"), name: "z")
  content("x.end", [$x$], anchor: "west")
  content("y.end", [$y$], anchor: "south")
  content("z.end", [$z$], anchor: "north-east")
  on-xy({ grid((0,0), (L,L), step: 1.0, stroke: gray + 0.2pt) })
  on-xz({ grid((0,0), (L,L), step: 1.0, stroke: black + 0.2pt) })
  on-yz({ grid((0,0), (L,L), step: 1.0, stroke: gray + 0.2pt) })
})
