#import "@preview/ctheorems:1.1.0": *


// double dot (second derivative)
#let ddot(cont) = $dot.double(cont)$

// imaginary unit 
#let ii = $upright(i)$

// exponential function
#let expp(cont) = $upright(e)^(cont)$

// Theorems ----------------------------
#let exercice = thmbox(
  "exercise", 
  "Exercise", 
  inset: 0em,
  base: none,
  breakable: true, // allow page break
).with(numbering: "1")

#let solution = thmplain(
  "solution",
  "Solution",
  base: "exercise",
  inset:  (top: 0.5em, left: 0.5em, right: 0.5em, bottom: 0.5em),
  breakable: true, // allow page break
  fill: luma(240),
).with(numbering: none)
// --------------------------------------