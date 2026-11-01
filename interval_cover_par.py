# Run the whole region of interval_cover.py in parallel, split into cuts of s, sy and t.
# Usage: python interval_cover_par.py <point file> [--workers 8] [--ns 4] [--minutes 600] [--out output JSON]
# If every sub-box is "finished, 0 violations, 0 limit boxes", the capture is >= T on the whole region (the sub-boxes are built so that their union is exactly the root box).
import sys, json, time, argparse
from fractions import Fraction as F
from multiprocessing import Pool
import interval_cover as ic


def job(args):
    path, T, minutes, sub = args
    L, cd, wd, pts = ic.load(path)
    r = ic.run(L, cd, wd, pts, T, minutes, sub=sub, log=lambda *a: None, max2=14)
    return dict(sub=[str(x) for x in sub], done=r["done"], nbox=r["nbox"], acc=r["acc"], acc2=r["acc2"], acc3=r["acc3"],
                frac=float(r["acc_frac"]), viol=len(r["viol"]), fails=len(r["fails"]), sec=r["sec"],
                fail_boxes=[[list(b), str(lb), str(ub), len(u)] for b, lb, ub, u in r["fails"][:5]])


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("path"); ap.add_argument("--T", default="1")
    ap.add_argument("--workers", type=int, default=8); ap.add_argument("--ns", type=int, default=4)
    ap.add_argument("--tcuts", default="0,1/1024,1/64,1/8,7/16", help="cuts of t (from 0 to 7/16)")
    ap.add_argument("--minutes", type=float, default=600, help="limit per sub-box")
    ap.add_argument("--out", default="interval_cover_par.json")
    a = ap.parse_args()
    L, cd, wd, pts = ic.load(a.path)
    dup, sym, mass = ic.check_data(L, cd, wd, pts)
    ic.need(sym, "not D4-invariant")
    half = F(L, 2) - F(cd, 2); smax = 1
    while smax < half: smax *= 2
    tc = [F(x) for x in a.tcuts.split(",")]
    ic.need(tc[0] == 0 and tc[-1] == F(7, 16) and all(x < y for x, y in zip(tc, tc[1:])), "the cuts of t must be an increasing sequence from 0 to 7/16")
    step = F(smax, a.ns)
    subs = []
    for i in range(a.ns):
        for j in range(a.ns):
            if i * step > half or j * step > half:  # sub-boxes not needed by symmetry (also discarded inside run)
                continue
            for k in range(len(tc) - 1):
                subs.append([i * step, (i + 1) * step, j * step, (j + 1) * step, tc[k], tc[k + 1]])
    # Boxes in the first t-cut (lower end 0) are classified excluding t = 0, so poses with theta = 0 are checked separately by axis-parallel enumeration
    ax = ic.axis_check(L, cd, wd, pts, F(a.T))
    axis = dict(ok=ax["ok"], min=str(ax["min"]), s=str(ax["s"]), sy=str(ax["sy"]), ncand=ax["ncand"])
    print(f"poses with theta = 0: minimum capture {float(ax['min']):.10f} (s = {float(ax['s'])}, sy = {float(ax['sy'])}, "
          f"candidates {ax['ncand'][0]}x{ax['ncand'][1]}) {'>= T' if ax['ok'] else '< T (violation)'}", flush=True)
    print(f"{len(subs)} sub-boxes, {a.workers} workers", flush=True)
    t0 = time.time(); res = []
    with Pool(a.workers) as pool:
        for r in pool.imap_unordered(job, [(a.path, F(a.T), a.minutes, s) for s in subs]):
            res.append(r)
            print(f"{time.time() - t0:7.0f}s [{len(res)}/{len(subs)}] {r['sub']} finished {r['done']} boxes {r['nbox']} "
                  f"volume {r['frac']:.4%} violations {r['viol']} limit {r['fails']} {r['sec']:.0f}s", flush=True)
            json.dump(dict(path=a.path, T=a.T, mass=str(mass), subs=len(subs), axis=axis, res=res), open(a.out, "w"), indent=1)
    ok = all(r["done"] and r["viol"] == 0 and r["fails"] == 0 for r in res) and len(res) == len(subs) and ax["ok"]
    print("Conclusion: capture >= T for every pose in the whole region (exact)" if ok else "Conclusion: some sub-boxes are undecided", flush=True)
