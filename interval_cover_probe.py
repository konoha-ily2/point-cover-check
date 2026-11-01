# Check that interval_cover.py does not accept a box it must not accept (counterexample search).
# Pick one pose P and compute its capture W(P) directly from the rotation formulas in exact rationals (without interval_cover's polynomials).
# Run `run` on a small box with P at a corner and threshold T = (W(P) + 1)/wd. Since the capture at P is < T, a correct checker
# must end with "violation", "split limit" or "time limit", and must not accept (finished, 0 violations, 0 limit boxes). If it accepts, there is a hole.
# Usage: python interval_cover_probe.py <point file> [--n 40] [--seed 1] [--minutes 1] [--ds 10] [--dt 30]
import argparse, random, time
from fractions import Fraction as F
import interval_cover as ic


def capture(pts, cd, s, sy, t):
    """Weight (exact integer in weight units) captured by the closed square at pose (s, sy, t) (coordinate units, t = tan(theta/2))."""
    c = (1 - t * t) / (1 + t * t); sn = 2 * t / (1 + t * t); H = F(cd, 2)
    hw = H * (c + sn); cx = hw + s; cy = hw + sy
    W = 0
    for X, Y, w in pts:
        dx, dy = X - cx, Y - cy
        if abs(dx * c + dy * sn) <= H and abs(-dx * sn + dy * c) <= H:
            W += w
    return W


def pick(rs, half, kind):
    """How to pick a pose: 0 = random, 1 = theta = 0 with the sides on 1/20 lines, 2 = small theta with the sides near the lines."""
    q = 1 << 10  # step of s: 2^-10 (coordinate units)
    if kind == 0:
        return F(rs.randrange(0, int(half) * q), q), F(rs.randrange(0, int(half) * q), q), F(rs.randrange(0, 7 << 26), 1 << 30)
    s = F(rs.randrange(0, int(half) // 50 + 1) * 50); sy = F(rs.randrange(0, int(half) // 50 + 1) * 50)
    if kind == 1:
        return s, sy, F(0)
    return s + F(rs.randrange(-4, 5), q), sy + F(rs.randrange(-4, 5), q), F(rs.randrange(0, 64), 1 << 30)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("path"); ap.add_argument("--n", type=int, default=40); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--minutes", type=float, default=1)
    ap.add_argument("--ds", type=int, default=10, help="width 2^-ds of the box in s, sy (coordinate units)")
    ap.add_argument("--dt", type=int, default=30, help="width 2^-dt of the box in t")
    a = ap.parse_args()
    L, cd, wd, pts = ic.load(a.path)
    half = F(L, 2) - F(cd, 2)
    rs = random.Random(a.seed); out = {"accepted (hole)": 0, "violation": 0, "limit": 0, "time limit": 0}; caps = []
    ds, dt = F(1, 1 << a.ds), F(1, 1 << a.dt)
    # Boxes whose lower end of t is 0 are classified excluding t = 0, so poses with theta = 0 are judged together with the axis-parallel enumeration (minimum over the whole container)
    ax_min = ic.axis_check(L, cd, wd, pts, F(1))["min"] * wd
    for i in range(a.n):
        s, sy, t = pick(rs, half, i % 3)
        s, sy = min(max(s, F(0)), half - ds), min(max(sy, F(0)), half - ds)
        W = capture(pts, cd, s, sy, t); caps.append(F(W, wd))
        # put P at one of the corners of the box (s, sy at the lower or upper end, t at the lower or upper end; t >= 0)
        bs = (s, s + ds) if rs.random() < 0.5 or s < ds else (s - ds, s)
        by = (sy, sy + ds) if rs.random() < 0.5 or sy < ds else (sy - ds, sy)
        bt = (t, t + dt) if rs.random() < 0.5 or t < dt else (t - dt, t)
        t0 = time.time()
        r = ic.run(L, cd, wd, pts, F(W + 1, wd), a.minutes, max_fail=1, log=lambda *x: None, sub=[*bs, *by, *bt], max2=14)
        res = ("violation" if r["viol"] else "limit" if r["fails"] else "time limit" if not r["done"] else
               "violation" if r["need_axis"] and ax_min < W + 1 else "accepted (hole)")
        out[res] += 1
        print(f"[{i}] kind {i % 3} s={float(s):.4f} sy={float(sy):.4f} t={float(t):.3g} capture {float(F(W, wd)):.6f} -> {res}"
              f" (boxes {r['nbox']}, {time.time() - t0:.1f}s)", flush=True)
    print("Result:", out, f"minimum capture {float(min(caps)):.6f}")
    print("No hole (0 accepted)" if out["accepted (hole)"] == 0 else "HOLE: accepted a box that must not be accepted")
