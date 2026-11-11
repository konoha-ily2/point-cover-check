# Run interval_cover_par.py on the certificates by other authors (and on the broken ones), one after another.
# Usage: python run_batch.py <certificate folder> [--workers 2] [--minutes 15] [--only name,...]
#   The folder holds the files converted by prep_cert.py (<name>.txt) and the originals (orig/<name>.txt).
#   Results go to other_certs/results/<name>.json, logs to other_certs/logs/<name>.log, the summary to other_certs/summary.json.
import sys, os, json, time, argparse, subprocess, hashlib
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import interval_cover as ic

# (name, expected, source). expected: pass = capture >= 1 at every pose; fail = some pose captures < 1 (a published control or our own mutant)
JOBS = [
    ("s12_56points_3.8", "pass", "evand@167d842"),
    ("s12_uniform_7of81_3.888", "pass", "evand@167d842"),
    ("s12_lower_3.931795_sparse", "pass", "evand@167d842"),
    ("s12_lower_3.931795", "pass", "evand@167d842"),
    ("s12_lower_3.9676", "pass", "evand@167d842"),
    ("s11_lower_3.8143", "pass", "evand@167d842"),
    ("s12_lower_3.9686", "pass", "evand@167d842"),
    ("s12_lower_3.969118", "pass", "squarepacker@7a96bec"),
    ("s12_lower_3.9715", "pass", "squarepacker@7a96bec"),
    ("56points_critical_760_199", "pass", "squarepacker@7a96bec controls (zero margin; should pass)"),
    ("56points_over_1520_397", "?", "squarepacker@7a96bec controls (fails at 0.8 on the angle grid; the exact θ=0 minimum is 1)"),
    ("scaled_further_31360_7900", "?", "squarepacker@7a96bec controls (grid: 0.98678 at θ≈3.05°; actual poses nearby give 1.0048)"),
    ("control1_lowered_orbit", "fail", "squarepacker@7a96bec controls/3.9715 (0.999985 at θ=0)"),
    ("control2_scaled_further", "fail", "squarepacker@7a96bec controls/3.9715 (0.6738 at θ=0)"),
    ("mut_3.931795_sparse_scale", "fail", "ours: weights of s12_lower_3.931795_sparse times 0.999 (0.9990 at θ=0)"),
    ("mut_3.9686_dropmax0", "fail", "ours: heaviest D4 orbit of s12_lower_3.9686 set to 0 (0.8192 at θ=0)"),
    ("mut_4.9950_rot", "fail", "ours: one orbit of s21_lower_4.9950 set to 0 (0.99962 at θ≈9.78°; θ=0 minimum stays 1.0010)"),
    ("s21_lower_4.9950", "pass", "evand@167d842 s21"),
    # second batch (run with --only). orig/ holds the text converted from the JSON (JSON SHA256 in fetch_and_prep.py)
    ("wand_n29_L557", "pass", "wand125/square-packing-bounds@678aa1d certificates/cert_n29_L557.json (angle grid, certificate for side B=0.9977)"),
    ("wand_n26_L545", "pass", "wand125/square-packing-bounds@678aa1d certificates/cert_n26_L545.json (same)"),
    ("wand_n39_L650", "pass", "wand125/square-packing-bounds@678aa1d certificates/cert_n39_L650.json (same)"),
]


def same_cert(a, b):
    """Whether the original and the converted file are the same as rational points and weights (and container side)."""
    def rd(p):
        tok = open(p).read().split(); sn, sd, cd, wd, cnt = map(int, tok[:5]); v = list(map(int, tok[5:]))
        d = {}
        for i in range(cnt):
            k = (F(v[3 * i], cd), F(v[3 * i + 1], cd)); d[k] = d.get(k, 0) + F(v[3 * i + 2], wd)
        return F(sn, sd), {k: w for k, w in d.items() if w}
    return rd(a) == rd(b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dir"); ap.add_argument("--workers", type=int, default=2); ap.add_argument("--minutes", type=float, default=15)
    ap.add_argument("--only", default=None)
    a = ap.parse_args()
    os.makedirs(os.path.join(HERE, "results"), exist_ok=True); os.makedirs(os.path.join(HERE, "logs"), exist_ok=True)
    sp = os.path.join(HERE, "summary.json")
    summ = json.load(open(sp, encoding="utf-8")) if os.path.exists(sp) else {}
    only = set(a.only.split(",")) if a.only else None
    env = dict(os.environ, PYTHONIOENCODING="utf-8")
    for name, expect, src in JOBS:
        if only and name not in only:
            continue
        path = os.path.join(a.dir, name + ".txt"); orig = os.path.join(a.dir, "orig", name + ".txt")
        if not os.path.exists(path):
            print("missing:", path, flush=True); continue
        L, cd, wd, pts = ic.load(path)
        dup, sym, mass = ic.check_data(L, cd, wd, pts)
        rec = dict(name=name, expect=expect, source=src, npts=len(pts), side=str(F(L, cd)), mass=str(mass), mass_f=float(mass),
                   sym=sym, dup=dup, sha256_orig=hashlib.sha256(open(orig, "rb").read()).hexdigest() if os.path.exists(orig) else None,
                   same_as_orig=same_cert(orig, path) if os.path.exists(orig) and not name.startswith("mut_") else None,
                   workers=a.workers, minutes_per_sub=a.minutes)
        out = os.path.join(HERE, "results", name + ".json"); log = os.path.join(HERE, "logs", name + ".log")
        t0 = time.time()
        with open(log, "w", encoding="utf-8") as f:
            p = subprocess.run([sys.executable, os.path.join(ROOT, "interval_cover_par.py"), path, "--workers", str(a.workers),
                                "--minutes", str(a.minutes), "--out", out], stdout=f, stderr=subprocess.STDOUT, env=env, cwd=ROOT)
        rec["wall_sec"] = round(time.time() - t0, 1); rec["exit"] = p.returncode
        if os.path.exists(out):
            r = json.load(open(out, encoding="utf-8")); res = r["res"]
            rec["axis"] = r["axis"]
            rec["subs"] = r["subs"]; rec["subs_done"] = sum(x["done"] for x in res)
            rec["viol_subs"] = sum(x["viol"] > 0 for x in res); rec["fail_subs"] = sum(x["fails"] > 0 for x in res)
            rec["timeout_subs"] = sum((not x["done"]) and x["viol"] == 0 and x["fails"] < 50 for x in res)
            rec["nbox"] = sum(x["nbox"] for x in res); rec["cpu_sec"] = round(sum(x["sec"] for x in res), 1)
            rec["open_subs"] = [x["sub"] for x in res if not (x["done"] and x["viol"] == 0 and x["fails"] == 0)]
            rec["fail_examples"] = [fb for x in res for fb in x["fail_boxes"]][:5]
            ok = len(res) == r["subs"] and r["axis"]["ok"] and not rec["open_subs"]
            rec["verdict"] = "pass" if ok else ("fail" if (not r["axis"]["ok"]) or rec["viol_subs"] else "open")
        else:
            rec["verdict"] = "error"
        summ[name] = rec
        json.dump(summ, open(sp, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
        print(f"{name}: {rec['verdict']} (expected {expect}) {rec['wall_sec']} s", flush=True)


if __name__ == "__main__":
    main()
