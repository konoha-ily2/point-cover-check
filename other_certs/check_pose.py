# Count, in exact rationals and directly from the original certificate file, the weight captured by the closed unit
# square at one pose (the polynomials of interval_cover are not used).
# Usage: python check_pose.py <certificate (original file)> <s> <sy> <t>
#   s, sy: at θ = 0 the positions of the left and bottom sides (in units of the side length, rationals).
#          The centre is (h + s, h + sy), h = (cos θ + sin θ)/2.
#   t = tan(θ/2) (rational). To read the values of an interval_cover box (s in units of 2^-20 coordinate units,
#   t in units of 2^-40), add --box E K cd.
#   Example: python check_pose.py scaled_further_31360_7900.txt 0 10498379116 126267424768 --box 20 40 7900
import sys, argparse, math
from fractions import Fraction as F

ap = argparse.ArgumentParser()
ap.add_argument("path"); ap.add_argument("s"); ap.add_argument("sy"); ap.add_argument("t")
ap.add_argument("--box", nargs=3, type=int, default=None, help="E K cd: read s, sy in units of 2^-E/cd and t in units of 2^-K")
a = ap.parse_args()
tok = open(a.path).read().split(); sn, sd, D, W, m = map(int, tok[:5]); v = list(map(int, tok[5:]))
pts = [(F(v[3 * i], D), F(v[3 * i + 1], D), F(v[3 * i + 2], W)) for i in range(m)]
L = F(sn, sd)
s, sy, t = F(a.s), F(a.sy), F(a.t)
if a.box:
    E, K, cd = a.box; s, sy, t = s / (2 ** E * cd), sy / (2 ** E * cd), t / 2 ** K
c = (1 - t * t) / (1 + t * t); sn_ = 2 * t / (1 + t * t)
h = (c + sn_) / 2; cx, cy = h + s, h + sy
half = F(1, 2)
verts = [(cx + dx * c - dy * sn_, cy + dx * sn_ + dy * c) for dx, dy in [(half, half), (-half, half), (-half, -half), (half, -half)]]
inside = all(0 <= x <= L and 0 <= y <= L for x, y in verts)
cap = sum(w for x, y, w in pts if abs((x - cx) * c + (y - cy) * sn_) <= half and abs(-(x - cx) * sn_ + (y - cy) * c) <= half)
print(f"container {L} ({float(L):.9f}), theta = {2 * math.degrees(math.atan(float(t))):.6f} deg, centre ({float(cx):.9f}, {float(cy):.9f})")
print(f"square inside the container: {inside}")
print(f"capture = {cap} = {float(cap):.10f}  {'< 1 (violation)' if cap < 1 else '>= 1'}")
