#import "@preview/cetz:0.4.1"
#import cetz.draw: *

#let ss = h(0.4em)

#let depth-center(depths) = {
  let sum = 0
  let result = ()

  for item in depths {
    result.push(sum + item / 2)
    sum += item
  }
  result
}

#let depth-edge(depths) = {
  let result = (0,)
  let sum = 0

  for item in depths {
    sum += item
    result.push(sum)
  }
  result
}


#let label(doc, outset: 2pt, color: black, ..args) = {
  box(fill: white, outset: outset, ..args)[#text(fill: color)[#doc]]
}

#let hline(width, y0, ..args) = line((0, y0), (width, y0), ..args)

// arrow: stealth, triangle, straight, barbed
// vertical arrow
#let v_arrow(height, x0, y0, name: "Q0", arrow: "straight", anchor: "west", ..args) = {
  line(
    (x0, y0),
    (x0, y0 - height),
    mark: (start: "stealth", fill: black),
    name: name,
    ..args,
  )
  content(name + ".mid", label(outset: -1pt, inset: 3pt)[#name], anchor: anchor)
}

// measurement arrow
#let mv_arrow(height, x0, y0, text: "Z1", arrow: "stealth", anchor: "east", ..args) = {
  let name = "mv"
  line(
    (x0, y0),
    (x0, y0 - height),
    mark: (end: arrow, start: arrow, fill: black, scale: 0.7),
    name: name,
    ..args,
  )
  content(name + ".mid", label(outset: -1pt, inset: 3pt)[#text], anchor: anchor)
}
