# Make a certificate with an error (to test that the checker finds the error).
# The weight of each of the 8 points in the D4 orbit of the point (X, Y) is reduced by delta. D4 invariance is preserved.
# Usage: python make_mutant.py <point file> <output> [--X 630 --Y 1000 --delta 3/2000]
# By default, the capture of the axis-parallel square with lower-left corner (0.0005, 0.0005) drops from 1.00110 to 0.99810.
import argparse
from fractions import Fraction as F

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("path"); ap.add_argument("out")
    ap.add_argument("--X", type=int, default=630); ap.add_argument("--Y", type=int, default=1000)
    ap.add_argument("--delta", default="3/2000", help="weight to subtract (rational)")
    a = ap.parse_args()
    tok = open(a.path).read().split()
    sn, sd, cd, wd, cnt = map(int, tok[:5]); v = list(map(int, tok[5:]))
    pts = [(v[3 * i], v[3 * i + 1], v[3 * i + 2]) for i in range(cnt)]
    L = F(sn, sd) * cd; assert L.denominator == 1; L = int(L)
    d = F(a.delta) * wd; assert d.denominator == 1; d = int(d)
    x, y = a.X, a.Y
    orb = {(x, y), (L - x, y), (x, L - y), (L - x, L - y), (y, x), (L - y, x), (y, L - x), (L - y, L - x)}
    hit = [p for p in pts if (p[0], p[1]) in orb]
    assert len(hit) == len(orb) and all(w >= d for _, _, w in hit), "the orbit points are not all present, or a weight is too small"
    out = [(X, Y, w - d if (X, Y) in orb else w) for X, Y, w in pts]
    cap = sum(w for X, Y, w in out if 1 <= X <= cd and 1 <= Y <= cd)  # axis-parallel square with lower-left corner (0.5, 0.5)/cd
    open(a.out, "w").write(f"{sn} {sd} {cd} {wd} {len(out)}\n" + "\n".join(f"{X} {Y} {w}" for X, Y, w in out) + "\n")
    print(f"reduced the weight of {len(orb)} points by {a.delta} each -> {a.out}. Capture of the axis-parallel square with lower-left corner (0.5, 0.5)/{cd}: {float(F(cap, wd)):.6f}")
