# Download the point-cover certificates (n=21, 45, 61) from pinned commits of Hosono's repository and check their SHA256.
# The certificates are not included in this repository (the data are Hosono's, MIT License; n=21 is derived from Evan Daniel's certificate).
# Usage: python fetch_certificate.py [--n 21 45 61] [--dir certificates]
import argparse, hashlib, os, sys, urllib.request

REPO = "wand125/square-packing-bounds"
CERTS = {
    21: ("39d8ecc74d651b54ec977c331c8f2015b442a6c4", "point_n21_L5/certificates/n21-capture-one.txt",
         "9f631fbae420e0376ea636a232b7680e937deebe205dc035cfb9c43a3d74b456", "n21-capture-one.txt"),
    45: ("678aa1d029423177792f51d3e7e42895b749241e", "point_n45_L7/cover.txt",
         "f7d706aa07506c351d9c4b76082cc391ae3f1fd20f6759496bd7c601bb0ca192", "n45-cover.txt"),
    61: ("678aa1d029423177792f51d3e7e42895b749241e", "point_n61_L8/cover.txt",
         "bcb66c7910ef844ce8d39c423d7a5a419bd530a344689331f7665ba6971d0374", "n61-cover.txt"),
}

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, nargs="+", default=sorted(CERTS), choices=sorted(CERTS))
    ap.add_argument("--dir", default="certificates")
    a = ap.parse_args()
    bad = 0
    for n in a.n:
        commit, path, sha, name = CERTS[n]
        out = os.path.join(a.dir, name)
        url = f"https://raw.githubusercontent.com/{REPO}/{commit}/{path}"
        if os.path.exists(out):
            data = open(out, "rb").read()
            print(f"already present: {out}")
        else:
            print(f"fetching: {url}")
            data = urllib.request.urlopen(url, timeout=60).read()
        h = hashlib.sha256(data).hexdigest()
        if h != sha:
            print(f"n={n}: SHA256 mismatch: {h} (expected {sha}). Do not use it"); bad += 1; continue
        os.makedirs(a.dir, exist_ok=True)
        if not os.path.exists(out):
            open(out, "wb").write(data)
        print(f"n={n}: SHA256 matches: {h} -> {out}")
    sys.exit(1 if bad else 0)
