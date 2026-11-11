# Convert wand125's point certificate (cert_nXX_LYYY.json, atoms = [[x, y, w], ...] as rationals) into the
# text format used by Hosono and Daniel (side_num side_den coord_den weight_den count, then X Y w). Values are unchanged.
# Usage: python json2txt.py <input.json> <output.txt>
import sys, json
from fractions import Fraction as F
from math import lcm

d = json.load(open(sys.argv[1]))
L = F(d["L"]); atoms = [(F(x), F(y), F(w)) for x, y, w in d["atoms"]]
cd = lcm(L.denominator, *[v.denominator for x, y, _ in atoms for v in (x, y)])
wd = lcm(*[w.denominator for _, _, w in atoms])
assert sum(w for *_, w in atoms) == F(d["total_mass"])
with open(sys.argv[2], "w") as f:
    f.write(f"{L.numerator} {L.denominator} {cd} {wd} {len(atoms)}\n")
    for x, y, w in atoms:
        f.write(f"{int(x * cd)} {int(y * cd)} {int(w * wd)}\n")
print(sys.argv[2], "cd", cd, "wd", wd, "atoms", len(atoms), "mass", float(F(d["total_mass"])))
