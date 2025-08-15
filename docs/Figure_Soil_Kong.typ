#import "@preview/cetz:0.4.1"
#import "main_cetz.typ": *

#show "{": ""
#show "}": ""

#let v_arrow(height, x0, y0, text: "Q0", arrow: "straight", anchor: "west", ..args) = {
  line(
    (x0, y0),
    (x0, y0 - height),
    mark: (end: "stealth", fill: black),
    name: "arrow",
    ..args,
  )
  content("arrow.mid", label(outset: -1pt, inset: 3pt)[#text], anchor: anchor)
}

#let depths = (1, 1, 1, 1, 1, 1, 1, 1, 1).map(x => x * 2.5)
#let zs = depth-center(depths).rev()
#let zh = depth-edge(depths).rev()

// #zs, #zh
#let z(i) = depths.at(i-1) / 2

#let z1 = zs.at(0)
#let z2 = zs.at(1)
#let z3 = zs.at(2)
#let z4 = zs.at(3)
#let z5 = zs.at(4)

#let z0_h = zh.at(0)
#let z1_h = zh.at(1)
#let z2_h = zh.at(2)
#let z3_h = zh.at(3)
#let z4_h = zh.at(4)
#let z5_h = zh.at(5)
#let z6_h = zh.at(6)
#let z7_h = zh.at(7)

#let dz = depths.at(0)

#let width = 8
#let x_mid = 0.42 * width

// #(zs: zs, zh: zh)
// #(z0_h: z0_h, z1_h: z0_h, z2_h: z2_h, z3_h: z3_h)
// #(_z2, _z1, z3)
// #v(-3.5em)

#figure(
  cetz.canvas(length: 8.5mm, {
    set-style(stroke: (thickness: 0.8pt))
    // ---------------- (a) i = 1 ----------------
    // content((0, z0_h + 0.6), [*(a)  i = 1*], anchor: "north-west")
    // rect((0, z1_h), (width, z2_h), fill: luma(95%), stroke: none, inset: 15pt) // 中心土壤层

    let ngz = 5
    
    // 土壤边界
    for i in range(0, ngz+1) {
      let k = i // + 1
      let _zh = zh.at(i)
      hline(width, _zh)
      if i > 0 {
        content((width, _zh), label(inset: 3pt, outset: 0pt)[$z_{#k+1/2} ss K_{#k+1/2}$], anchor: "west")
      }
    }

    // 土壤中心
    for i in range(0, ngz) {
      let k = i + 1
      let _z = zs.at(i)
      hline(width, _z, stroke: (dash: "dashed"))
      content((x_mid, _z), label[$theta_#k ss psi_#k ss K_#k ss z_#k$], anchor: "center")
    }

    // z_{k+1/2}, K_{k+1/2}
    content((x_mid, z0_h), label[$theta_0 ss psi_0 ss K_0$], anchor: "center")

    content((width, z0_h), label(inset: 3pt, outset: 0pt)[$z_0$], anchor: "west") 
    content((width * 0.89, z0_h + 0.05 * dz), label[地表], anchor: "south")
    // content((width * 0.89, z5_h + 0.06 * dz), label[深层], anchor: "south")

    // 绘制地下水水位
    // let z_wt = (z3_h + z5_h * 1.18) / 2
    let jwt = 4
    let z_wt = z0_h - dz * (jwt - 1) - 0.25 * dz 
    rect((0, z5_h), (width, z_wt), fill: blue.transparentize(90%), stroke: none, inset: 15pt)
    // rect((0, z0_h), (width, z_wt), fill: blue.transparentize(90%), stroke: none, inset: 15pt)
    hline(width, z_wt, stroke: (dash: "dashed", paint: blue))
    content((width, z_wt), label(inset: 3pt, outset: 0pt, color: blue)[$z_"wtb"$], anchor: "west")

    v_arrow(0.5 * dz, 0.72 * width, z0_h, text: $Q_0$)
    v_arrow(dz, 0.75 * width, z1, text: $Q_1$)
    v_arrow(dz, 0.78 * width, z2, text: $Q_2$)
    v_arrow(dz, 0.83 * width, z3, text: $Q_{3}$)
    v_arrow(dz, 0.86 * width, z4, text: $Q_{4}$)
    v_arrow(0.5 * dz, 0.89 * width, z5, text: $Q_{5}$)

    mv_arrow(dz, -0.02 * width, z0_h, text: $Delta z_1$)
    mv_arrow(dz, -0.02 * width, z1_h, text: $Delta z_2$)
    mv_arrow(dz, -0.02 * width, z2_h, text: $Delta z_3$)
    mv_arrow(dz, -0.02 * width, z3_h, text: $Delta z_4$)
    mv_arrow(dz, -0.02 * width, z4_h, text: $Delta z_5$)

    let dy_h = 0.1
    mv_arrow(0.5 * dz, 0.04 * width, z0_h, text: $Delta z_{1/2}$, anchor: "west")
    mv_arrow(dz, 0.04 * width, z1, text: $Delta z_{1+1/2}$, anchor: "west")
    mv_arrow(dz, 0.04 * width, z2, text: $Delta z_{2+1/2}$, anchor: "west")
    mv_arrow(dz, 0.04 * width, z3, text: $Delta z_{3+1/2}$, anchor: "west")
    mv_arrow(dz, 0.04 * width, z4, text: $Delta z_{4+1/2}$, anchor: "west")
    mv_arrow(0.5 * dz, 0.04 * width, z5, text: $Delta z_{5+1/2}$, anchor: "west")
  }),
  caption: [土壤层结构示意图。（该情景地下水水位在第4层）
  ],
) <fig_soil>
