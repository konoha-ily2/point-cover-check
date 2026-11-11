# Build a "broken certificate" whose only certain violations are at poses with θ > 0 (a negative control).
# 1. Compute the capture W(P) of random poses (θ in [θ0, θ1]) exactly in rationals and collect the small ones.
# 2. Set one whole D4 orbit of points captured by P to 0, and look for one with W'(P) < 1 while the θ = 0 minimum
#    (axis_check) stays >= 1. Such a certificate passes the θ = 0 enumeration and certainly fails at P (θ > 0),
#    so it tests whether the box search of the checker finds the violation.
# Usage: python make_rot_mutant.py <converted certificate> <output> [--n 4000] [--seed 1] [--th0 5] [--th1 40]
import sys, os, argparse, random, json
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import interval_cover as ic
from interval_cover_probe import capture


def orbit(x, y, Lc):
    s = set()
    for (a, b) in [(x, y), (y, x)]:
        for (c, d) in [(a, b), (Lc - a, b), (a, Lc - b), (Lc - a, Lc - b)]:
            s.add((c, d))
    return frozenset(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inp"); ap.add_argument("out")
    ap.add_argument("--n", type=int, default=4000); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--th0", type=float, default=5); ap.add_argument("--th1", type=float, default=40)
    a = ap.parse_args()
    L, cd, wd, pts = ic.load(a.inp)
    import math
    rs = random.Random(a.seed); half = F(L - cd, 2)
    t0, t1 = math.tan(math.radians(a.th0) / 2), math.tan(math.radians(a.th1) / 2)
    cand = []
    for _ in range(a.n):
        t = F(rs.randrange(int(t0 * 2 ** 20), int(t1 * 2 ** 20)), 2 ** 20)
        c = (1 - t * t) / (1 + t * t); sn = 2 * t / (1 + t * t)
        room = F(L) - cd * (c + sn)  # upper limit of s (stays inside the container)
        s = F(rs.randrange(0, int(min(room, half) * 1024) + 1), 1024); sy = F(rs.randrange(0, int(min(room, half) * 1024) + 1), 1024)
        cand.append((capture(pts, cd, s, sy, t), s, sy, t))
    cand.sort(key=lambda x: x[0])
    d = {}
    for x, y, w in pts: d[(x, y)] = d.get((x, y), 0) + w
    for W, s, sy, t in cand[:50]:
        c = (1 - t * t) / (1 + t * t); sn = 2 * t / (1 + t * t); H = F(cd, 2); hw = H * (c + sn); cx = hw + s; cy = hw + sy
        inside = [(x, y) for (x, y), w in d.items() if w and abs((x - cx) * c + (y - cy) * sn) <= H and abs(-(x - cx) * sn + (y - cy) * c) <= H]
        orbs = sorted({orbit(x, y, L) for x, y in inside}, key=lambda o: -sum(d[p] for p in o))
        for o in orbs:
            lost = sum(d[p] for p in o if p in inside)
            if W - lost >= wd:
                continue
            nd = {p: (0 if p in o else w) for p, w in d.items()}
            npts = [(x, y, w) for (x, y), w in sorted(nd.items())]
            ax = ic.axis_check(L, cd, wd, npts, F(1))
            if not ax["ok"]:
                continue
            Wn = capture(npts, cd, s, sy, t)
            assert Wn < wd
            with open(a.out, "w") as f:
                f.write(f"{L} {cd} {cd} {wd} {len(npts)}\n")
                for x, y, w in npts: f.write(f"{x} {y} {w}\n")
            info = dict(pose=dict(s=str(s), sy=str(sy), t=str(t), theta_deg=2 * math.degrees(math.atan(float(t)))),
                        capture_before=str(F(W, wd)), capture_after=str(F(Wn, wd)), capture_after_f=float(F(Wn, wd)),
                        orbit=str(sorted(o)[0]), orbit_mass=str(F(sum(d[p] for p in o), wd)), axis_min_after=str(ax["min"]),
                        axis_min_after_f=float(ax["min"]))
            json.dump(info, open(a.out + ".json", "w"), indent=1)
            print(json.dumps(info, ensure_ascii=False)); return
    print("not found")


if __name__ == "__main__":
    main()
