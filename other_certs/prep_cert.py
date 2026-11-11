# Convert a point-weight certificate by another author into the form read by interval_cover.py (values unchanged).
# Usage: python prep_cert.py <input> <output> [--mutant kind value]
#  - Multiply the coordinate denominator D and the coordinates by k (X/D = kX/kD; the side s_num/s_den is unchanged).
#    k is chosen as described in main.
#  - Hosono's format (side_num side_den coord_den weight_den count) and Daniel's (s_num s_den D W m) have the same layout.
#  - --mutant builds a negative control (a certificate broken on purpose):
#      scale p/q   multiply every weight by p/q (weight denominator times q, numerators times p)
#      dropmax k   set the k-th heaviest D4 orbit (0-based) to 0
#      shift dx    move the D4 orbit with the smallest positive x away from the centre by dx coordinate units
import sys, argparse
from fractions import Fraction as F


def read(path):
    tok = open(path).read().split()
    sn, sd, cd, wd, cnt = map(int, tok[:5]); v = list(map(int, tok[5:]))
    assert len(v) == 3 * cnt
    return sn, sd, cd, wd, [(v[3 * i], v[3 * i + 1], v[3 * i + 2]) for i in range(cnt)]


def orbit(x, y, Lc):
    s = set()
    for (a, b) in [(x, y), (y, x)]:
        for (c, d) in [(a, b), (Lc - a, b), (a, Lc - b), (Lc - a, Lc - b)]:
            s.add((c, d))
    return frozenset(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp"); ap.add_argument("out")
    ap.add_argument("--mutant", nargs=2, default=None)
    a = ap.parse_args()
    sn, sd, cd, wd, pts = read(a.inp)
    # Multiply the coordinate denominator by k (same values). k is the least integer such that cd*k is even and the
    # root-box condition of interval_cover.run holds (smax = least power of two >= half, smax + 1.41422*cd <= L).
    L0 = F(sn, sd)
    k = 1
    while True:
        c = cd * k
        half = (L0 * c - c) / 2; smax = 1
        while smax < half: smax *= 2
        if c % 2 == 0 and (L0 * c).denominator == 1 and smax + F(141422, 100000) * c <= L0 * c:
            break
        k += 1
    cd *= k; pts = [(k * x, k * y, w) for x, y, w in pts]
    Lc = F(sn, sd) * cd
    assert Lc.denominator == 1; Lc = int(Lc)
    note = "same"
    if a.mutant:
        kind, val = a.mutant
        if kind == "scale":
            p, q = map(int, val.split("/")); wd *= q; pts = [(x, y, w * p) for x, y, w in pts]
            note = f"weights x {p}/{q}"
        else:
            d = {}
            for x, y, w in pts: d[(x, y)] = w
            orbs = {}
            for (x, y) in d: orbs.setdefault(orbit(x, y, Lc), 0)
            for o in orbs: orbs[o] = sum(d.get(p, 0) for p in o)
            if kind == "dropmax":
                j = int(val); o = sorted(orbs, key=lambda o: (-orbs[o], sorted(o)))[j]
                for p in o: d[p] = 0
                note = f"orbit {sorted(o)[0]} (mass {F(orbs[o], wd)}) set to 0"
            elif kind == "shift":
                dx = int(val)
                # take the weighted orbit farthest from the centre (nearest a corner) and move it dx toward the corner (mapped by D4)
                o = max((o for o in orbs if orbs[o] > 0), key=lambda o: (min(abs(2 * x - Lc) + abs(2 * y - Lc) for x, y in o), sorted(o)))
                nd = {p: w for p, w in d.items() if p not in o}
                for (x, y) in o:
                    nx = x - dx if 2 * x < Lc else x + dx
                    ny = y - dx if 2 * y < Lc else y + dx
                    nx = min(max(nx, 0), Lc); ny = min(max(ny, 0), Lc)
                    nd[(nx, ny)] = nd.get((nx, ny), 0) + d[(x, y)]
                d = nd
                note = f"orbit {sorted(o)[0]} moved {dx}/{cd} toward the corner"
            else:
                sys.exit("unknown mutant")
            pts = [(x, y, w) for (x, y), w in sorted(d.items())]
    with open(a.out, "w") as f:
        f.write(f"{sn} {sd} {cd} {wd} {len(pts)}\n")
        for x, y, w in pts: f.write(f"{x} {y} {w}\n")
    print(a.out, "k", k, "cd", cd, "pts", len(pts), note)


if __name__ == "__main__":
    main()
