# Tests for interval_cover.py: (1) classifications agree with direct computation, (2) conclusions on small certificates match the known answers.
# Usage: python interval_cover_tests.py
import math, random
from fractions import Fraction as F
import interval_cover as ic


def pose_of(s, sy, t, cd):
    th = 2 * math.atan(t)
    h = cd / 2 * (math.cos(th) + math.sin(th))
    return h + s, h + sy, th


def inside_float(X, Y, cx, cy, th, cd):
    dx, dy = X - cx, Y - cy
    u = dx * math.cos(th) + dy * math.sin(th); v = -dx * math.sin(th) + dy * math.cos(th)
    return abs(u) <= cd / 2 + 1e-9 and abs(v) <= cd / 2 + 1e-9, max(abs(u), abs(v)) - cd / 2


def box_pose(box, rs):
    """A random pose (s, sy, t) in the box (s = S0 + S1*t, sy = Y0 + Y1*t)."""
    s0, s1, u0, u1, y0, y1, v0, v1, m0, m1 = box
    t = (m0 + rs.random() * (m1 - m0)) / (1 << ic.K)
    S0 = s0 + rs.random() * (s1 - s0); S1 = u0 + rs.random() * (u1 - u0)
    Y0 = y0 + rs.random() * (y1 - y0); Y1 = v0 + rs.random() * (v1 - v0)
    return (S0 + S1 * t) / (1 << ic.E), (Y0 + Y1 * t) / (1 << ic.E), t


def test_classify(trials=3000, seed=1, wedge=False):
    """Whether points classified "inside" / "outside" for a box really are so at random poses in the box (0 contradictions expected).
    With wedge=True, test on wedge-coordinate boxes (t near 0, s = e +- t*u around an integer e, etc.)."""
    rs = random.Random(seed); cd = 1000; H = (cd // 2) << ic.E; bad = 0; cnt = [0, 0, 0]
    for _ in range(trials):
        if not wedge:
            ws = rs.choice([1, 8, 64, 512]) << (ic.E - 3); wt = rs.choice([1, 16, 256]) << (ic.K - 14)
            s0 = rs.randrange(0, 2000 << ic.E); y0 = rs.randrange(0, 2000 << ic.E); m0 = rs.randrange(0, 6 << (ic.K - 4))
            box = (s0, s0 + ws, 0, 0, y0, y0 + ws, 0, 0, m0, m0 + wt)
        else:
            e = rs.randrange(1, 2000) << ic.E; f = rs.randrange(1, 2000) << ic.E
            wt = rs.choice([1, 16, 256]) << (ic.K - 34); m0 = rs.choice([0, wt])
            uw = rs.choice([1, 64, 4096]) << ic.E; u0 = rs.randrange(-ic.U, ic.U - uw + 1); v0 = rs.randrange(-ic.U, ic.U - uw + 1)
            if rs.random() < 0.5:
                box = (e, e, u0, u0 + uw, f, f, v0, v0 + uw, m0, m0 + wt)
            else:
                w = 1 << (ic.E - 6)
                box = (e, e + w, ic.U, ic.U, f - w, f, -ic.U, -ic.U, m0, m0 + wt)
        M = ic.basis(box[8], box[9] - box[8])
        s, sy, t = box_pose(box, rs)
        cx, cy, _ = pose_of(s, sy, t, cd)
        X = int(cx + rs.uniform(-800, 800)); Y = int(cy + rs.uniform(-800, 800))
        if wedge and rs.random() < 0.7:  # more points whose sides lie on grid lines
            X = rs.choice([e >> ic.E, (e >> ic.E) + cd, X]); Y = rs.choice([f >> ic.E, (f >> ic.E) + cd, Y])
        c = ic.classify(X, Y, box, M, H); cnt[c + 1] += 1
        for _ in range(20):
            ins, marg = inside_float(X, Y, *pose_of(*box_pose(box, rs), cd), cd)
            if (c == 1 and not ins) or (c == -1 and ins and marg < -1e-6):
                bad += 1
    print(f"classification check{' (wedge coordinates)' if wedge else ''}: outside {cnt[0]}, undecided {cnt[1]}, inside {cnt[2]}, contradictions with random poses {bad}")
    return bad == 0


def test_bip_identity(trials=2000, seed=2):
    """Whether the opposite-sides sum formula (top of stage_bip) agrees exactly with adding the two constraints directly at the same pose."""
    rs = random.Random(seed); H = 500 << ic.E; bad = 0
    ev = lambda e, t: sum(c * t ** k for k, c in enumerate(e))
    for _ in range(trials):
        s, sy = rs.randrange(0, 3000 << ic.E), rs.randrange(0, 3000 << ic.E); t = F(rs.randrange(0, 1 << 20), 1 << 21)
        Xp, Yp, Xq, Yq = (rs.randrange(0, 5001) for _ in range(4))
        ap, bp = (Xp << ic.E) - s, (Yp << ic.E) - sy; aq, bq = (Xq << ic.E) - s, (Yq << ic.E) - sy
        dx, dy = (Xq - Xp) << ic.E, (Yq - Yp) << ic.E
        lhs = ev(ic.polys(ap, bp, H)[3], t) + ev(ic.polys(aq, bq, H)[2], t)
        rhs = ev((dy - 2 * H, -2 * dx, -4 * H, -2 * dx, -dy - 2 * H), t)
        lhs2 = ev(ic.polys(ap, bp, H)[1], t) + ev(ic.polys(aq, bq, H)[0], t)
        rhs2 = ev((dx - 2 * H, 2 * dy, -4 * H, 2 * dy, -dx - 2 * H), t)
        bad += (lhs != rhs) + (lhs2 != rhs2)
    print(f"opposite-sides sum formula: {trials} pairs, mismatches {bad}")
    return bad == 0


def test_bip_sound(trials=400, seed=3):
    """In boxes accepted by stage_bip, whether the capture (floating point) at random poses in the box is really >= T."""
    rs = random.Random(seed); cd = 1000; H = (cd // 2) << ic.E; bad = 0; nacc = 0
    for _ in range(trials):
        ws = rs.choice([1, 8, 64]) << (ic.E - 3); wt = rs.choice([1, 16, 256]) << (ic.K - 16)
        s0 = rs.randrange(0, 2000 << ic.E); y0 = rs.randrange(0, 2000 << ic.E); m0 = rs.choice([0, rs.randrange(0, 1 << (ic.K - 6))])
        box = (s0, s0 + ws, 0, 0, y0, y0 + ws, 0, 0, m0, m0 + wt); M = ic.basis(m0, wt)
        cx, cy, _ = pose_of(s0 / (1 << ic.E), y0 / (1 << ic.E), m0 / (1 << ic.K), cd)
        # put points close to the sides (on grid lines)
        pts = []
        for _ in range(rs.randrange(4, 40)):
            side = rs.randrange(4); u = int(rs.uniform(-520, 520)); o = int(rs.choice([-500, 500]) + rs.uniform(-3, 3))
            X, Y = (int(cx) + u, int(cy) + o) if side < 2 else (int(cx) + o, int(cy) + u)
            pts.append((X, Y, rs.randrange(1, 10)))
        und = [i for i in range(len(pts)) if ic.classify(*pts[i][:2], box, M, H) == 0]
        sure = sum(pts[i][2] for i in range(len(pts)) if ic.classify(*pts[i][:2], box, M, H) == 1)
        tot = sure + sum(pts[i][2] for i in und)
        if tot == 0:
            continue
        Tw = rs.randrange(max(1, tot // 3), tot + 1)
        if not ic.stage_bip(und, pts, box, M, H, tot - Tw):
            continue
        nacc += 1
        for _ in range(30):
            s = (s0 + rs.random() * ws) / (1 << ic.E); sy = (y0 + rs.random() * ws) / (1 << ic.E); t = (m0 + rs.random() * wt) / (1 << ic.K)
            pose = pose_of(s, sy, t, cd)
            cap = sum(w for X, Y, w in pts if inside_float(X, Y, *pose, cd)[0])
            if cap < Tw:
                bad += 1
    print(f"opposite-sides stage check: accepted {nacc} boxes, cases with capture < T at a random pose {bad}")
    return bad == 0 and nacc > 0


def test_bound_exact(trials=3000, seed=5):
    """Whether the nested upper and lower bounds bracket the exact value of q(r) at rational r in [0,1]."""
    rs = random.Random(seed); bad = 0
    for _ in range(trials):
        n = rs.randrange(1, ic.D + 2)
        e = [rs.randrange(-10 ** 6, 10 ** 6) * rs.choice([0, 1, 1]) for _ in range(n)]
        m0 = rs.randrange(0, 7 << (ic.K - 4)); dm = rs.choice([1, 1 << 20, 1 << 36]); M = ic.basis(m0, dm)
        up, lo = ic.bound(e, M, True), ic.bound(e, M, False)
        for r in [F(0), F(1)] + [F(rs.randrange(0, 1 << 30), 1 << 30) for _ in range(8)]:
            t = F(m0, 1 << ic.K) + F(dm, 1 << ic.K) * r
            v = sum(c * t ** k for k, c in enumerate(e)) * (1 << (ic.D * ic.K))
            bad += not (lo <= v <= up)
    print(f"nested upper/lower bounds: {trials} polynomials, cases where the exact value fell outside {bad}")
    return bad == 0


def test_pos_open(trials=4000, seed=7):
    """Whenever pos_open (the "outside" test excluding t = 0, for boxes whose lower end of t is 0) is true, whether q(r) > 0 at rational r in (0, 1] (exact).
    Many polynomials with d0 = 0 (a point exactly on a side at t = 0) are generated."""
    rs = random.Random(seed); bad = 0; npos = 0
    for _ in range(trials):
        n = rs.randrange(2, ic.D + 2)
        e = [rs.randrange(-10 ** 6, 10 ** 6) * rs.choice([0, 1, 1]) for _ in range(n)]
        if rs.random() < 0.6:
            e[0] = 0
        dm = rs.choice([1, 1 << 20, 1 << 36]); M = ic.basis(0, dm)
        if not ic.pos_open(e, M):
            continue
        npos += 1
        for r in [F(1), F(1, 1 << 40)] + [F(rs.randrange(1, 1 << 30), 1 << 30) for _ in range(8)]:
            t = F(dm, 1 << ic.K) * r
            bad += not sum(c * t ** k for k, c in enumerate(e)) > 0
    print(f"outside test excluding t = 0: among {npos} polynomials where it was true, cases with q(r) <= 0 for some r in (0,1]: {bad}")
    return bad == 0 and npos > 0


def test_axis_check(trials=60, seed=8):
    """Whether the minimum of axis_check is <= the capture at random axis-parallel poses (rationals, many exactly at break values) and is actually attained at its pose."""
    rs = random.Random(seed); bad = 0
    for _ in range(trials):
        L = rs.choice([2000, 2500, 3000]); cd = 1000
        pts = [(rs.randrange(0, L + 1), rs.randrange(0, L + 1), rs.randrange(1, 10)) for _ in range(rs.randrange(1, 30))]
        if rs.random() < 0.5:  # snap onto grid lines
            pts = [(x - x % 250, y - y % 250, w) for x, y, w in pts]
        ax = ic.axis_check(L, cd, 1, pts, F(1))
        cap = lambda s, sy: sum(w for x, y, w in pts if s <= x <= s + cd and sy <= y <= sy + cd)
        bad += cap(ax["s"] * cd, ax["sy"] * cd) != ax["min"]
        for _ in range(200):
            def pick():
                if rs.random() < 0.5:
                    v = rs.choice([x for x, _, _ in pts] + [y for _, y, _ in pts]) - rs.choice([0, cd])
                    return min(max(F(v), F(0)), F(L - cd))
                return F(rs.randrange(0, (L - cd) * 64 + 1), 64)
            bad += cap(pick(), pick()) < ax["min"]
    print(f"theta = 0 enumeration: {trials} certificates x 200 poses, captures below the minimum or wrong minimiser {bad}")
    return bad == 0


def test_apex_cover(trials=3000, seed=6):
    """Whether the union of the boxes re-split into wedge coordinates contains every pose (s, sy, t) of the original box (exact in rationals, by solving the re-split formulas)."""
    rs = random.Random(seed); bad = 0; n = 0; one = 1 << ic.E; w = 1 << (ic.E - 6)
    smax_e = 2048 << ic.E
    for _ in range(trials):
        def end():
            e = rs.randrange(1, 1900) * one
            return rs.choice([(e, e + w), (e - w, e), (e + 5, e + w - 3)])
        (s0, s1), (y0, y1) = end(), end()
        m1 = rs.choice([1, 1 << 10, 1 << 22])
        box = (s0, s1, 0, 0, y0, y1, 0, 0, 0, m1)
        kids = ic.apex_pieces(box, smax_e)
        if kids is None:
            continue
        n += 1
        for _ in range(10):
            s = s0 + F(rs.randrange(0, 1 << 20), 1 << 20) * (s1 - s0) if rs.random() < 0.8 else rs.choice([s0, s1])
            y = y0 + F(rs.randrange(0, 1 << 20), 1 << 20) * (y1 - y0) if rs.random() < 0.8 else rs.choice([y0, y1])
            t = rs.choice([F(0), F(m1, 1 << ic.K), F(rs.randrange(0, m1 + 1), 1 << ic.K)])
            def has(lo, hi, w0, w1, x):
                # whether there are S0 in [lo, hi] and S1 in [w0, w1] with x = S0 + S1*t (t >= 0, S0 + S1*t is monotone)
                return lo + w0 * t <= x <= hi + w1 * t
            bad += not any(has(k[0], k[1], k[2], k[3], s) and has(k[4], k[5], k[6], k[7], y) for k in kids)
    print(f"wedge-coordinate coverage: {n} re-split boxes x 10 poses, poses that fall in no child {bad}")
    return bad == 0 and n > 0


def cert(L, pts, wd=1):
    return L, 1000, wd, pts


def d4(L, base):
    out = {}
    for x, y, w in base:
        for (a, b) in [(x, y), (L - x, y), (x, L - y), (L - x, L - y), (y, x), (L - y, x), (y, L - x), (L - y, L - x)]:
            out[(a, b)] = w
    return [(a, b, w) for (a, b), w in sorted(out.items())]


def case(name, L, pts, T, expect):
    L, cd, wd, pts = cert(L, pts)
    dup, sym, mass = ic.check_data(L, cd, wd, pts)
    r = ic.run(L, cd, wd, pts, F(T), 5, log=lambda *a: None)
    ax = ic.axis_check(L, cd, wd, pts, F(T))
    ok_all = r["done"] and not r["viol"] and not r["fails"] and ax["ok"]
    got = "holds for all poses" if ok_all else ("violation found" if r["viol"] or not ax["ok"] else "undecided")
    print(f"{name}: points {len(pts)} mass {float(mass)} D4 {sym} -> {got} (boxes {r['nbox']}, {r['sec']:.1f}s) expected {expect}")
    return got == expect


if __name__ == "__main__":
    ok = test_classify()
    ok &= test_classify(seed=4, wedge=True)
    ok &= test_bip_identity()
    ok &= test_bip_sound()
    ok &= test_bound_exact()
    ok &= test_apex_cover()
    ok &= test_pos_open()
    ok &= test_axis_check()
    # (A) 4 points on a side-2 container at (0.9 +- 0.1): every unit square has at least one point in its interior -> capture >= 1, mass 4 < 5 so s(5) >= 2 (weaker than the known fact, but correct)
    ok &= case("A side 2, 4 points (correct)", 2000, d4(2000, [(900, 900, 1)]), 1, "holds for all poses")
    # (B) 4 points at (0.5, 0.5) etc. on a side-2 container: a square tilted by 45 degrees at the centre captures none of the points -> a violation must be found
    ok &= case("B side 2, 4 corner points (wrong)", 2000, d4(2000, [(500, 500, 1)]), 1, "violation found")
    # (C) 1 point at the centre of a side-2 container: correct (every unit square contains the centre). Zero slack; the 10/31 bounds did not decide it,
    #     but the nested bounds (11/01) close the corner at t = 0
    ok &= case("C side 2, 1 centre point (correct, zero slack)", 2000, [(1000, 1000, 1)], 1, "holds for all poses")
    # (D) 1 point shifted by 0.001 from the centre of a side-2 container: the square [0,1]^2 in the lower-left corner misses it -> a violation must be found (not D4-invariant, but the violation shows up in the reduced range alone)
    ok &= case("D side 2, 1 point shifted 0.001 from the centre (wrong)", 2000, [(1001, 1001, 1)], 1, "violation found")
    print("All as expected" if ok else "Some results differ from the expected ones")
