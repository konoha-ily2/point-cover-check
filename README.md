# point-cover-check

An independent, exact re-check of the computational step of point-weight
lower-bound certificates for packing unit squares in a square:
**every closed unit square inside the side-*k* container captures point mass ≥ 1**.
It has been run on the three point-only certificates published by
Hosono (wand125):

| *n* | container side | points | total mass *M* | bound it gives |
|---|---|---|---|---|
| 21 | 5 | 4,604 | 20.99999200 | s(21) ≥ 5 |
| 45 | 7 | 12,645 | 44.99999100 | s(45) ≥ 7 |
| 61 | 8 | 15,193 | 60.99998780 | s(61) ≥ 8 |

All three pass for all poses (all sub-boxes finished, no violation, no box left
at the split limit).

> **Status.** Computer-assisted check written with an AI coding assistant
> (Claude Code) and tested only by the tests in this repository.
> It has **not** been reviewed by anyone else. It is not a proof-assistant proof.

## What this checks, and what it does not

A lower bound s(*n*) ≥ *k* from a point-weight certificate has three parts:

1. **Data**: the weights are non-negative, the total mass is *M* < *n*, and the
   weights are invariant under the symmetry group D4 of the container.
2. **Capture**: every closed unit square (any position, any angle) contained in
   the closed square [0, *k*]² captures weight ≥ 1.
3. **Scaling argument**: if *n* unit squares with disjoint interiors fit in a
   square of side *L* < *k*, enlarge by *k*/*L* and take the concentric closed unit
   squares strictly inside the enlarged ones; they are pairwise disjoint, so
   *n* ≤ *M*, a contradiction.

This repository checks **(1) and (2)** in exact integer and rational
arithmetic. It does **not** check (3); see Hosono's READMEs (and, for *n* = 21,
his Lean files).

The claim is therefore limited to:

> The capture step (2) of each of the three certificates below
> (SHA256 as listed) holds for all poses, as checked by code written
> independently of Hosono's verifier and of Evan Daniel's `zmx2`.

- The certificate data are Hosono's. The *n* = 21 file is derived from Evan
  Daniel's earlier point certificate (see Hosono's attribution). Nothing here
  improves them.
- The values themselves were already published: s(21) = 5 by Evan Daniel
  (mixed point-and-segment proof) and by Hosono (point-only certificate with a
  rational replay and a Lean reduction); s(45) = 7 and s(60) = 8 (hence
  s(61) = 8) by Evan Daniel with mixed measures; the point-only *n* = 45 and
  *n* = 61 certificates by Hosono, checked with Daniel's `zmx2`.
  This is a re-check, not a new result.

## Certificates

The certificates are **not redistributed here**. `fetch_certificate.py` downloads
them from Hosono's repository <https://github.com/wand125/square-packing-bounds>
(MIT License) at pinned commits and checks their hashes:

| *n* | commit | file | SHA256 |
|---|---|---|---|
| 21 | `39d8ecc74d651b54ec977c331c8f2015b442a6c4` | `point_n21_L5/certificates/n21-capture-one.txt` | `9f631fbae420e0376ea636a232b7680e937deebe205dc035cfb9c43a3d74b456` |
| 45 | `678aa1d029423177792f51d3e7e42895b749241e` | `point_n45_L7/cover.txt` | `f7d706aa07506c351d9c4b76082cc391ae3f1fd20f6759496bd7c601bb0ca192` |
| 61 | `678aa1d029423177792f51d3e7e42895b749241e` | `point_n61_L8/cover.txt` | `bcb66c7910ef844ce8d39c423d7a5a419bd530a344689331f7665ba6971d0374` |

All three use the same text format (Hosono's `point_n21_L5/FORMAT.md`):
container side as a fraction, coordinate denominator, weight denominator,
number of entries, then one `X Y w` line per point.

| *n* | side | coordinate denominator | weight denominator | total mass |
|---|---|---|---|---|
| 21 | 5/1 | 1000 | 999948000000 | 2624862500021/124993500000 |
| 45 | 7/1 | 2000 | 2⁴⁹ | 12666371418707823/2⁴⁸ |
| 61 | 8/1 | 2000 | 2⁴⁹ | 8584985072679551/2⁴⁷ |

The *n* = 21 file is Hosono's normalised one (weights divided by his capture
threshold 249987/250000). For all three, the threshold to check is exactly 1.
Hosono's repository has since moved to `wand125/square-packing`; the pinned
commits above stay readable in the original repository.

## Method (summary)

**Poses.** Let θ be the angle of the square and t = tan(θ/2), so that
cos θ = (1−t²)/(1+t²) and sin θ = 2t/(1+t²). The horizontal half-extent of the
rotated square is h(t) = (cos θ + sin θ)/2. Write the centre as
(c_x, c_y) = (h + s, h + s_y); then s, s_y ≥ 0 are exactly the conditions that
the square does not cross the left and bottom walls.

**Symmetry reduction.** Because the weights are D4-invariant, it suffices to
check c_x ≤ L/2, c_y ≤ L/2 and θ ∈ [0, π/4], i.e. t ∈ [0, tan(π/8)].
The checker uses t ∈ [0, 7/16] and s, s_y ∈ [0, s_max] (in coordinate units),
with s_max the least power of two ≥ (L − 1)/2 (2048 for *n* = 21, 8192 for
*n* = 45 and 61). This covers the region; every pose in it lies inside the
container (checked).

**Containment as polynomials.** A point (X, Y) is in the closed square iff four
polynomials in t of degree 4, with integer coefficients that are affine in
a = X − s and b = Y − s_y, are all ≤ 0 (the four sides, multiplied by (1+t²)²).
For t ∈ [0, 1) the sign of the coefficient of a and of b in each polynomial is
fixed, so over a box the maximum and minimum are attained at known corners in
(a, b).

**Bounds.** On a t-interval [t₀, t₀ + w], substitute t = t₀ + w·r with r ∈ [0, 1]
and write the polynomial as d₀ + r(d₁ + r(d₂ + …)). An upper bound is obtained
from the inside out using r·x ≤ max(0, x) (lower bound with min). All of this
is done in exact integers: box endpoints are dyadic rationals
(2⁻²⁰ in position, 2⁻⁴⁰ in t).

**Branch and bound.** For each box, a point is *in* (all four upper bounds ≤ 0),
*out* (some lower bound > 0) or undecided. A box is accepted when the weight of
points that are *in* for the whole box is ≥ 1. A box is a **violation** when even
the in-or-undecided weight is < 1 (then every pose in it captures < 1).
Otherwise one of the following exact stages is tried before splitting:

- **Opposite sides (bipartite stage).** For two points p, q at the same pose,
  (bottom constraint of p) + (top constraint of q) does not depend on s, s_y:
  it is a polynomial in t and the difference q − p only. If it is ≤ 0 on the
  box's t-interval, p cannot leave through the bottom while q leaves through the
  top. The same holds for left + right. Each undecided point gets one copy per
  constraint it may violate; left-side copies (bottom, left) and right-side
  copies (top, right) form a bipartite graph whose edges are these conflicts and
  the pairs of copies of the same point. The weight that can be outside at one
  pose is at most the maximum-weight independent set of this graph,
  computed as total weight − maximum flow (König). Soundness needs only that
  the computed flow is feasible.
- **Pair conflicts (small sets only).** Two points cannot both be outside if a
  positive combination λ·c_p + μ·c_q of their violated constraints is ≤ 0 on the
  box; the multiplier is *chosen* in floating point but *checked* exactly. Then
  an exact branch-and-bound bounds the maximum-weight independent set.
- **Wedge coordinates.** Near t = 0, when a box edge in s (or s_y) is an integer
  grid position e, the box is re-parametrised as s = e + t·u (u ∈ [0, U]) plus
  s = S₀ + U·t (S₀ in the original interval). For every t the union of the two
  pieces contains the original interval; the constraints stay polynomial in t
  (degree 5). The new boxes are not subsets of the original box, so **every
  point is re-classified from scratch** in them.

**θ = 0 is checked separately.** When a side of the square lies exactly on a
line of certificate points at θ = 0, those points are *in* at θ = 0 (the square
is closed) but *out* for every small θ > 0. A box whose t-interval contains 0
then keeps them undecided however finely it is split. So:

- In a box with t₀ = 0, *out* means "out for every pose with t > 0": with
  q(r) = d₀ + r·Q(r), it holds if d₀ ≥ 0 and the lower bound of Q is > 0
  (or if the closed-interval test already holds). *In* is still tested on the
  closed interval.
- The poses with θ = 0 (axis-parallel closed squares [s, s+1] × [s_y, s_y+1]
  anywhere in the container) are checked by exact enumeration. Membership in x
  changes only at s = X − 1 and just after s = X, so the capture is constant on
  each cell formed by these break values and the open intervals between them, in
  both axes. Every pair of candidates (break values and midpoints) is evaluated
  with a 2-D prefix sum in 64-bit integers (the total weight is checked to be
  < 2⁶²); the minimising pose is recounted with Python integers.

Floating-point numbers are used only to choose the split direction, the split
limit, when to switch to wedge coordinates, and candidate multipliers. They never
decide acceptance: every acceptance test is an exact integer inequality.

The region is run as independent sub-boxes: 4 × 4 in (s, s_y), × 4 intervals in
t with cuts at 0, 1/1024, 1/64, 1/8, 7/16. Sub-boxes with every s (or s_y) above
L/2 are skipped by symmetry (64 sub-boxes for *n* = 21 and 61, 36 for *n* = 45).
Their union is exactly the root box minus that symmetric part.
The result counts only if the θ = 0 enumeration gives a minimum ≥ 1 **and**
**every** sub-box finishes with no violation and no box at the split limit.

## Reproduce

Python ≥ 3.8 and NumPy (tested with CPython 3.10.5 and NumPy on Windows).
NumPy is used only for the θ = 0 enumeration, with 64-bit integers.

```sh
python fetch_certificate.py                # downloads to certificates/ and checks SHA256
python interval_cover_tests.py             # unit tests and small certificates
python interval_cover_par.py certificates/n21-capture-one.txt --workers 8 --ns 4 --out run_n21.json
python interval_cover_par.py certificates/n45-cover.txt --workers 8 --ns 4 --out run_n45.json
python interval_cover_par.py certificates/n61-cover.txt --workers 8 --ns 4 --out run_n61.json
```

Each run must print the θ = 0 line with `>= T` and end with the line
`Conclusion: capture >= T for every pose in the whole region (exact)`.
The JSON must have `axis.ok: true` and list every sub-box with `done: true`,
`viol: 0`, `fails: 0`.

Our runs (2026-11-01, AMD Ryzen 7 5800X, 15 workers, one certificate after another):

| *n* | sub-boxes | boxes | accepted (bipartite / pair stage) | CPU s | wall s | θ = 0 minimum capture (at s, s_y) |
|---|---|---|---|---|---|---|
| 21 | 64 | 103,472 | 51,264 (39,601 / 1,228) | 3,330 | 411 | 1.0011020573 (0, 0) |
| 45 | 36 | 126,204 | 62,918 (47,149 / 734) | 5,516 | 691 | 1.0032959479 (1.003, 2.06225) |
| 61 | 64 | 76,865 | 38,330 (31,269 / 123) | 4,857 | 518 | 1.0116913222 (0, 3.002) |

The records are `results/interval_par_n21.json`, `results/interval_par_n45.json`
and `results/interval_par_n61.json`.

Do not use `python -O`. (The checks that matter for soundness raise errors
explicitly, but there is no reason to remove the remaining assertions either.)

## Tests

`interval_cover_tests.py` (all counts of mismatches must be 0, and the last line must be `All as expected`):

- classification of points as in/out against random poses in random boxes,
  including wedge-coordinate boxes (3,000 boxes × 20 poses each);
- the opposite-side identity, compared with direct evaluation in exact rationals;
- boxes accepted by the bipartite stage, sampled at random poses;
- the nested bound against exact values of random polynomials at rational points;
- exact coverage of wedge re-parametrisation (every sampled pose of the original
  box lies in one of the new boxes);
- the t > 0 *out* test: whenever it holds, the polynomial is > 0 at exact
  rational r ∈ (0, 1] (1,498 polynomials, many with d₀ = 0);
- the θ = 0 enumeration against direct counting at random axis-parallel poses,
  many of them exactly at break values (60 small certificates × 200 poses);
- four tiny certificates with known answers (two correct, two wrong).

`interval_cover_probe.py <certificate>` looks for false acceptances on a real
certificate: it picks a pose P, computes its capture W(P) directly from the
rotation (without the checker's polynomials), and asks the checker to certify
capture ≥ W(P) + (one weight unit) on a small box with P at a corner,
together with the θ = 0 enumeration. A correct checker must never accept. ("0 accepted" means that no counterexample was found; it is evidence, not a proof that the checker is correct.)
Probes mix random poses, axis-parallel poses with sides on the 1/20 grid, and
near-axis-parallel poses. With the current code (2026-11-01):

- *n* = 21: 60 probes, 58 violations, 2 at the split limit, **0 accepted**;
- *n* = 45: 60 probes, 60 violations, **0 accepted**.

An earlier version (without the separate θ = 0 treatment) gave 0 accepted in
366 probes on *n* = 21.

**Mutation test.** We reduced the weight of one D4 orbit of 8 points
(e.g. (0.63, 1.0), weight 0.0509 each) by 0.0015 in the *n* = 21 certificate.
The axis-parallel square with lower-left corner (0, 0) then captures
0.998102 < 1. `make_mutant.py` makes this file.
The full parallel run on it (current code) reports a θ = 0 minimum of 0.998102
at s = s_y = 0, and 29 violation boxes in 5 of the 64 sub-boxes (all near the
corner, s, s_y < 1.024); it does not conclude. On the unmodified file the same
run concludes.
The record is `results/interval_par_n21_mutant.json`.

## Limitations

- Written with an AI coding assistant and checked only by its own tests;
  no external review.
- Checks the capture step only; the scaling argument is not checked here.
- Three certificates only (*n* = 21, 45, 61), all by Hosono.
- Pure Python apart from the NumPy enumeration; each full run needs about
  7–12 minutes on 15 workers.
- Relies on CPython's arbitrary-precision integers and `fractions.Fraction`.

## License and attribution

Cite as: Hirotaka Nakano, *point-cover-check*, Zenodo, DOI 10.5281/zenodo.23853803 (concept DOI, all versions); v0.1.0: 10.5281/zenodo.23853804.

Code: MIT License, Copyright (c) 2026 Hirotaka Nakano.
The certificates are by Hosono (wand125), MIT License; the *n* = 21 one is
derived from Evan Daniel's certificate (MIT License). They are downloaded,
not included.
