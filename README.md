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

Since v0.2.0 the same checker, unchanged, has also been run on 14 point
certificates by other authors (*n* = 11, 12, 21, 26, 29, 39; all pass) and on
7 negative controls (5 known-bad: all rejected; of 2 published controls whose
status was unknown, one is now a confirmed counterexample and one is
undecided). See [Other certificates](#other-certificates-v020).

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

## Other certificates (v0.2.0)

Run on 2026-11-05 with the checker of v0.1.0 (the acceptance logic of
`interval_cover.py` and `interval_cover_par.py` is unchanged). The tools are in
[`other_certs/`](other_certs/); the records are `other_certs/summary.json`
(one entry per certificate) and `other_certs/results/*.json` (per sub-box).

### Sources

All downloads are pinned to a commit and checked against SHA256
(listed in `other_certs/fetch_and_prep.py`); every file matched.
The files are **not redistributed here**.

- **Evan Daniel** (`evand/square-packing` at `167d842`, `s12/certificates/`):
  seven *n* = 12 certificates, *n* = 11 (680 points) and *n* = 21
  (`s21/s21_lower_4.9950.txt`, 4,604 points).
- **squarepacker** (Ryu Sungjoon, `s12-lower-bound` at `7a96bec`):
  *n* = 12 at 3.969118 and 3.9715 (the registry's verified lower bound
  7943/2000), and the five published controls in `controls/`.
- **wand125** (`square-packing-bounds` at `678aa1d`, `certificates/cert_nXX_LYYY.json`):
  point certificates for *n* = 26, 29, 39. The originals are stated with an
  angle grid and a square of side B = 0.9977; the claim checked here (side 1)
  is weaker than that.

All are D4-invariant, without duplicate points, with non-negative weights.

### Conversion (values unchanged)

`prep_cert.py` multiplies the coordinate denominator D and the coordinates by
an integer k (1 ≤ k ≤ 6): the checker needs D even, and needs its root box
(s up to the least power of two ≥ (L − 1)/2, plus the diagonal) to fit in the
container. `json2txt.py` turns wand125's JSON into the common text format.
`run_batch.py` checks every time that the converted file equals the original
as a set of rational points and weights, with the same container side.

### Results

2 workers, 8 minutes per sub-box. Wall time is for 2 workers;
CPU is the sum over sub-boxes. "θ = 0 minimum" is the exact value from the
axis-parallel enumeration.

| certificate | *n* | points | side | total mass | θ = 0 minimum | sub-boxes | boxes | wall / CPU | result |
|---|---|---|---|---|---|---|---|---|---|
| Daniel `s12_56points_3.8` | 12 | 56 | 19/5 | 11.2 | 1 (zero margin) | 64 | 700 | 1.8 s / 2 s | pass |
| Daniel `s12_uniform_7of81_3.888` | 12 | 81 | 35/9 | 11.571 | 1 (zero margin) | 36 | 1,409 | 2.7 s / 4 s | pass |
| Daniel `s12_lower_3.931795_sparse` | 12 | 224 | 3920/997 | 11.983 | 1.0000022 | 36 | 4,789 | 13 s / 24 s | pass |
| Daniel `s12_lower_3.931795` | 12 | 788 | 3920/997 | 11.933 | 1.0000023 | 36 | 8,334 | 64 s / 125 s | pass |
| Daniel `s12_lower_3.9676` | 12 | 764 | 980/247 | 11.996 | 1.0000033 | 36 | 16,003 | 112 s / 220 s | pass |
| Daniel `s11_lower_3.8143` | 11 | 680 | 3040/797 | 10.815 | 1.0000042 | 36 | 10,997 | 86 s / 163 s | pass |
| Daniel `s12_lower_3.9686` | 12 | 1,736 | 15680/3951 | 11.974 | 1.0000056 | 36 | 28,908 | 330 s / 643 s | pass |
| squarepacker `s12_lower_3.969118` | 12 | 1,736 | 31360/7901 | 11.974 | 1.0000056 | 36 | 31,254 | 338 s / 659 s | pass |
| squarepacker `s12_lower_3.9715` | 12 | 1,736 | 7943/2000 | 11.997 | 1.0000050 | 36 | 32,627 | 350 s / 684 s | pass |
| squarepacker control `56points_critical_760_199` | 12 | 56 | 760/199 | 11.2 | 1 (zero margin) | 64 | 911 | 1.9 s / 3 s | pass (as expected) |
| Daniel `s21_lower_4.9950` | 21 | 4,604 | 5000/1001 | 20.805 | 1.0010177 | 64 | 132,331 | 1,636 s / 3,123 s | pass |
| wand125 `cert_n29_L557` | 29 | 748 | 557/100 | 28.993 | 1.0000030 | 64 | 12,990 | 51 s / 97 s | pass |
| wand125 `cert_n26_L545` | 26 | 1,376 | 109/20 | 25.856 | 1.0000101 | 36 | 17,634 | 126 s / 236 s | pass |
| wand125 `cert_n39_L650` | 39 | 2,724 | 13/2 | 38.519 | 1.0000102 | 36 | 47,017 | 291 s / 548 s | pass |

"pass" means the same as above: the θ = 0 enumeration gives ≥ 1, and every
sub-box finished with no violation and no box at the split limit. The three
certificates with zero margin at θ = 0 also close.

### Negative controls

| certificate | how it is broken | known violation | wall | result | decided by |
|---|---|---|---|---|---|
| squarepacker `control1_lowered_orbit` | heaviest orbit lowered by 100/10⁷ | 0.999985 at θ = 0 | 1,051 s | rejected | θ = 0 enumeration: 199997/200000 at s = s_y = 0 (the box search ran out of memory, see below) |
| squarepacker `control2_scaled_further` | same integers with D = 3999600 | 0.6737611 at θ = 0 | 334 s | rejected | θ = 0 enumeration: 6737611/10⁷ at s = s_y = 0 (box search ran out of memory; see the note below) |
| ours `mut_3.931795_sparse_scale` | all weights × 0.999 | 0.9990022 at θ = 0 | 7 s | rejected | θ = 0 enumeration and violation boxes (20 sub-boxes) |
| ours `mut_3.9686_dropmax0` | heaviest orbit set to 0 | 0.8191874 at θ = 0 | 285 s | rejected | θ = 0 enumeration and violation boxes (8 sub-boxes) |
| ours `mut_4.9950_rot` | one orbit of `s21_lower_4.9950` ((3978, 4522), mass 0.0175) set to 0 | 0.9996176 **only at θ ≈ 9.78°**; the θ = 0 minimum stays 1.0010177 | 616 s | rejected | **box search alone** (5 sub-boxes, 11 violation boxes) |
| squarepacker control `scaled_further_31360_7900` | 3.969118 widened one more step | unknown (see below) | 318 s | **rejected: counterexample confirmed** | box search (5 sub-boxes) |
| squarepacker control `56points_over_1520_397` | 56 points widened one more step | unknown (see below) | 2.7 s | undecided | 50 boxes at the split limit in 1 sub-box |

For `control2` the parallel run stopped with `MemoryError` before writing its
JSON, so `summary.json` lists it as `error`. The θ = 0 line printed before the
crash (minimum 0.6737611000 at s = s_y = 0, 1831 × 1831 candidates) is the
verdict; it is reproduced by the θ = 0 enumeration alone.

`mut_4.9950_rot` is made by `make_rot_mutant.py`: it computes the exact
capture at 1,500 random poses (θ ∈ [5°, 40°]), removes one orbit captured by a
low pose P, and keeps a mutant with capture at P < 1 and θ = 0 minimum ≥ 1. It
cannot be rejected by the θ = 0 enumeration, so it tests whether the box search
finds a violation at θ > 0. (The same construction found no suitable orbit for
`s12_lower_3.931795_sparse`, whose θ = 0 margin is only 2.2·10⁻⁶.)
`interval_cover_probe.py` on wand125 *n* = 29 (30 probes, 15 s each) gave
30 violations, 0 accepted.

**`scaled_further_31360_7900` is not a point cover.** At a corner of a
violation box returned by the box search we recounted the capture from the
original file with the rotation formula (`check_pose.py`, not the checker's
polynomials): θ = 13.1022° (t = 126267424768/2⁴⁰), left vertex on the wall
x = 0 (s = 0, s_y = 10498379116/(2²⁰·7900) ≈ 1.26735), centre
(0.600328, 1.867674). This closed unit square lies in the container and
captures 4997147/5000000 = 0.9994294 < 1; 50 other poses in the same box give
the same value. The README of that repository (v1.1) correctly warns that a
grid failure does not by itself show that a real unit square captures < 1, and
reports 1.0048 for actual poses near the grid point (θ ≈ 3.05°). That remark
stands; the counterexample found here is elsewhere (θ ≈ 13°, against a wall).

```sh
python other_certs/check_pose.py <work>/orig/scaled_further_31360_7900.txt 0 10498379116 126267424768 --box 20 40 7900
```

**`56points_over_1520_397` is undecided.** Its FORMAT notes that an
axis-parallel square captures only 4 of 5 points; the exact θ = 0 minimum is
5/5 (by the enumeration and by a brute-force rational count), so the 0.8 reads
as a test with a shrunk square. The box search stopped in the sub-box
s ∈ [0, 0.43], s_y ∈ [0.43, 0.86], θ ≥ 14.25°, near s = 0, s_y ≈ 0.736524,
θ ≈ 36.87° (t ≈ 1/3), with lower bound 4/5, upper bound 6/5 and 2 points
undecided. There the point (0.94458, 0.99496) lies exactly on the bottom side
and (0.99496, 1.87657) exactly on the right side. 1,600 random poses around it
(radius 10⁻³ to 10⁻⁹) capture 5 or 6 points, never 4. Adjacent sides (bottom
and right) are not handled by the opposite-side stage, and the pair stage did
not close it with its limit raised to 60. Whether the true minimum is exactly 1
or 4/5 at some pose is not known.

### Limitations of this run

- `control1` and `control2` (coordinate denominator 4·10⁶) are rejected by the
  θ = 0 enumeration, but the box search near the counterexample neither accepts
  nor rejects boxes and ran out of memory with 2 workers (`control1`: 2
  sub-boxes reached 76 % of their volume in 8 minutes, about 660,000 boxes
  each). The verdict comes from the exact enumeration; memory use in deep boxes
  has not been fixed or measured. The correct `s12_lower_3.9715` with the same
  denominator passes in 350 s, so a large denominator by itself is not the
  obstacle.
- The wedge constants (U = 8192, re-split size 1/64) and the split limit are
  in coordinate units and were tuned for D = 1000. They affect only efficiency,
  not soundness; whether they are efficient for other denominators has not
  been checked.
- As above, only the capture step is checked: not the scaling argument, and
  not the total mass < *n* (all are below *n*, see the table). None of these
  bounds is new.

### Reproduce (AMD Ryzen 7 5800X, Python 3.10, 2 workers)

```sh
python other_certs/fetch_and_prep.py work      # download, SHA256, convert, build the mutants (~4 min)
python other_certs/run_batch.py work --workers 2 --minutes 8   # all 21 (~1 h 40 min in total)
python other_certs/run_batch.py work --only s12_56points_3.8,mut_4.9950_rot
python other_certs/check_pose.py work/mut_4.9950_rot.txt 2012089/512 4088695/1024 44837/524288 --box 0 0 2002
```

Repeating step 1 in a fresh folder gave byte-identical converted files for all 21.

## Limitations

- Written with an AI coding assistant and checked only by its own tests;
  no external review.
- Checks the capture step only; the scaling argument is not checked here.
- Full-region runs of the three Hosono certificates (*n* = 21, 45, 61) and of
  the 14 certificates listed under [Other certificates](#other-certificates-v020).
- Pure Python apart from the NumPy enumeration; each full run needs about
  7–12 minutes on 15 workers.
- Relies on CPython's arbitrary-precision integers and `fractions.Fraction`.

## License and attribution

Cite as: Hirotaka Nakano, *point-cover-check*, Zenodo, DOI 10.5281/zenodo.23853803 (concept DOI, all versions); v0.1.0: 10.5281/zenodo.23853804.

Code: MIT License, Copyright (c) 2026 Hirotaka Nakano.
The certificates are by Hosono (wand125), MIT License; the *n* = 21 one is
derived from Evan Daniel's certificate (MIT License). They are downloaded,
not included.
