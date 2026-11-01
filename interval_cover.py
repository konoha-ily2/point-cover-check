# Exact check of a point-cover certificate (lower bound for packing unit squares in a square),
# using integer / rational interval arithmetic only, no floating point.
# What is checked: for every closed unit square contained in the container of side L, the captured weight is >= T.
# Usage: python interval_cover.py <point file> [--T 1] [--minutes 25] [--tmin-num 0 --tmin-den 1]
# The point file has the same format as Hosono's FORMAT.md (side_num side_den coord_den weight_den count, then X Y w lines).
# Method and limitations: see README.md.
#
# How a pose is described (assuming D4-symmetric weights, reduced to 1/8: cx <= L/2, cy <= L/2, theta in [0, pi/4])
#   t = tan(theta/2) in [0, 7/16] (contains tan(pi/8) = 0.41421...), cos = (1-t^2)/(1+t^2), sin = 2t/(1+t^2)
#   h(t) = (cos+sin)/2 (horizontal distance from the centre to the leftmost vertex), cx = h + s, cy = h + sy
#   (s, sy >= 0 is exactly the condition that the square fits in the container)
#   At theta = 0, s is the position of the left side and sy of the bottom side. Coordinates are integers in units of
#   1/coord_den, s and sy are dyadic with step 2^-E, t has step 2^-K.
# A point (X, Y) is inside the square iff the following four quartics in t are all <= 0
# (a = X - s, b = Y - sy, H = half side = coord_den/2)
#   c1 = (a-2H) + (2b-4H)t - 4H t^2 + (2b+4H)t^3 + (-a-2H)t^4      (right side)
#   c2 = -a + (4H-2b)t + 0 t^2 + (-2b-4H)t^3 + a t^4                (left side)
#   c3 = (b-2H) - 2a t + 4H t^2 - 2a t^3 + (-b-2H)t^4               (top side)
#   c4 = -b + 2a t - 8H t^2 + 2a t^3 + b t^4                        (bottom side)
# (u*(1+t^2)^2 and v*(1+t^2)^2 expanded, where u, v are the point's coordinates in the square's frame.)
# The signs of the coefficients of a and b are constant for t in [0,1), so over a box the maximum and minimum
# are attained at the ends of (a, b).
#
# A box has the form s = S0 + S1*t, sy = Y0 + Y1*t (S0, S1, Y0, Y1, t are each intervals). Ordinary boxes have S1 = Y1 = 0.
# In a box with t near 0 where s or sy is near an integer (a grid position), the position at which a side crosses a
# grid line is proportional to (s - integer)/t and spreads over 0..infinity inside the box (a wedge). So a small box
# with an integer e at its end is re-split into s = e + t*u (u in [0, U]) and s = e + U*t + r (r in [0, width])
# ("wedge coordinates", `apex` in run).
# In either form a(t) = (X - S0) - S1*t is linear in t, so c_i is a quintic in t and the choice of ends is unchanged.
# For t, write t = t0 + w*r (r in [0,1]) and q(r) = sum d_j r^j in nested form d0 + r(d1 + r(d2 + ...)); from the
# inside out, bound using r*X <= max(0, X) to get the upper (lower) bound. Everything is computed in exact integers.
import sys, time, argparse
from fractions import Fraction as F
from math import comb

E = 20   # step of s, sy: 2^-E (coordinate units)
K = 40   # step of t: 2^-K
D = 5    # maximum degree of the polynomials in t (5 with wedge coordinates)
U = 8192 << E  # upper limit of u in wedge coordinates (coordinate units; the effect of the slope of a side, 2*a*t, has |a| <= about 2000, so u up to about 4000 suffices)


def need(cond, msg):
    """Check a precondition that soundness depends on. Unlike assert, it is not removed by python -O."""
    if not cond:
        raise ValueError(msg)


def load(path):
    tok = open(path).read().split()
    sn, sd, cd, wd, cnt = map(int, tok[:5]); v = list(map(int, tok[5:]))
    need(min(sn, sd, cd, wd) > 0 and cnt >= 0 and len(v) == 3 * cnt, "bad point file format")
    pts = [(v[3 * i], v[3 * i + 1], v[3 * i + 2]) for i in range(cnt)]
    need(all(w >= 0 for _, _, w in pts), "negative weight")
    L = F(sn, sd) * cd
    need(L.denominator == 1, "the container side must be an integer in coordinate units")
    return int(L), cd, wd, pts


def check_data(L, cd, wd, pts):
    """Exact: coordinate range, duplicates, D4 invariance, total mass."""
    d = {}
    for x, y, w in pts:
        need(0 <= x <= L and 0 <= y <= L, "point outside the container")
        d[(x, y)] = d.get((x, y), 0) + w
    dup = len(d) != len(pts)
    sym = all(d.get((L - x, y)) == w and d.get((y, x)) == w for (x, y), w in d.items())
    mass = F(sum(w for _, _, w in pts), wd)
    return dup, sym, mass


def basis(m0, dm):
    """M[k][j] = C(k,j) m0^(k-j) dm^j 2^((D-k)K): coefficients that rewrite t^k*2^(DK) as a polynomial in r."""
    return [[comb(k, j) * m0 ** (k - j) * dm ** j * (1 << ((D - k) * K)) if j <= k else 0 for j in range(D + 1)] for k in range(D + 1)]


def bound(e, M, upper):
    """Upper (lower) bound of q(r) = d0 + r(d1 + r(d2 + ...)) on r in [0,1], using r*X <= max(0, X) from the inside out.
    (Always <= the earlier bound d0 + max(0, d1 + sum max(d_j, 0)). A corner where a point sits on a vertex,
    d0 = d1 = 0 and d2 < 0 < d3, can also be classified as inside.)"""
    n = len(e)
    d = [sum(e[k] * M[k][j] for k in range(j, n)) for j in range(n)]
    acc = d[-1]
    for j in range(n - 2, -1, -1):
        acc = d[j] + (max(0, acc) if upper else min(0, acc))
    return acc


def pos_open(e, M):
    """Whether q(r) > 0 for every r in (0, 1] (exact). Used for boxes whose lower end of t is 0, to classify a point as
    "outside" for all t > 0, excluding t = 0.
    Write q(r) = d0 + r*Q(r); if d0 >= 0 and the lower bound of Q is > 0, then q(r) >= r*Q(r) > 0 for r > 0. The closed-interval
    test (lower bound of `bound` > 0) is included.
    At a pose with theta = 0 where a side lies on a grid line, the constraint of a point on that line is exactly 0 at t = 0
    and positive for t > 0. With the closed interval, a box containing t = 0 would stay undecided however finely it was
    split (this is what blocked n=45). Poses with t = 0 are checked separately by axis_check."""
    n = len(e)
    d = [sum(e[k] * M[k][j] for k in range(j, n)) for j in range(n)]
    if n == 1:
        return d[0] > 0
    acc = d[-1]
    for j in range(n - 2, 0, -1):
        acc = d[j] + min(0, acc)
    return (d[0] > 0 and d[0] + min(0, acc) > 0) or (d[0] >= 0 and acc > 0)


def polys(a, b, H):
    """Coefficients in t of the four constraints (all integers, scaled by 2^E)."""
    return (
        (a - 2 * H, 2 * b - 4 * H, -4 * H, 2 * b + 4 * H, -a - 2 * H),
        (-a, 4 * H - 2 * b, 0, -2 * b - 4 * H, a),
        (b - 2 * H, -2 * a, 4 * H, -2 * a, -b - 2 * H),
        (-b, 2 * a, -8 * H, 2 * a, b),
    )


def polys_t(a0, a1, b0, b1, H):
    """Coefficients in t (up to degree 5) of the four constraints when a = a0 - a1*t, b = b0 - b1*t."""
    p = polys(a0, b0, H)
    if a1 == 0 and b1 == 0:
        return p
    pa = polys(1, 0, 0); pb = polys(0, 1, 0)  # with H = 0 only the coefficients of a and b remain
    out = []
    for i in range(4):
        e = list(p[i]) + [0]
        for k in range(5):
            e[k + 1] -= a1 * pa[i][k] + b1 * pb[i][k]
        out.append(tuple(e))
    return out


# Ends that give the maximum: c1 (a large, b large), c2 (a small, b small), c3 (a small, b large), c4 (a large, b small). The minimum is the opposite.
HI_A = (True, False, False, True)
HI_B = (True, False, True, False)


def ends(X, Y, box):
    """Large and small ends of a(t) and b(t), as (constant, coefficient of t) pairs. Since t >= 0, the smaller S0 and S1 give the larger a."""
    s0, s1, u0, u1, y0, y1, v0, v1 = box[:8]
    XE, YE = X << E, Y << E
    return ((XE - s1, u1), (XE - s0, u0)), ((YE - y1, v1), (YE - y0, v0))  # (small, large)


def cons_at(X, Y, box, H, i, upper):
    (alo, ahi), (blo, bhi) = ends(X, Y, box)
    a = ahi if HI_A[i] == upper else alo
    b = bhi if HI_B[i] == upper else blo
    return polys_t(a[0], a[1], b[0], b[1], H)[i]


def classify(X, Y, box, M, H):
    """1 = inside for every pose in the box, -1 = outside for every pose, 0 = undecided.
    In a box whose lower end of t is 0, -1 means "outside for every pose with t > 0" (poses with t = 0 are checked separately by axis_check)."""
    open0 = box[8] == 0
    for i in range(4):
        e = cons_at(X, Y, box, H, i, False)
        if pos_open(e, M) if open0 else bound(e, M, False) > 0:
            return -1
    for i in range(4):
        if bound(cons_at(X, Y, box, H, i, True), M, True) > 0:
            return 0
    return 1


def viol_cons(X, Y, box, M, H):
    """Indices of the constraints that may be violated somewhere in the box (upper bound > 0)."""
    return [i for i in range(4) if bound(cons_at(X, Y, box, H, i, True), M, True) > 0]


def plain(box):
    return box[2] == box[3] == box[6] == box[7] == 0


def terms(X, Y, i, box, H):
    """Coefficients of constraint c_i split as c = e0(t) - rho*ws*ea(t) - sigma*wy*eb(t) (rho, sigma in [0,1]); ordinary boxes only."""
    s0, s1, y0, y1 = box[0], box[1], box[4], box[5]
    p0 = polys(0, 0, H)[i]
    ea = [x - z for x, z in zip(polys(1, 0, H)[i], p0)]
    eb = [x - z for x, z in zip(polys(0, 1, H)[i], p0)]
    e0 = polys((X << E) - s0, (Y << E) - y0, H)[i]
    return e0, [-(s1 - s0) * x for x in ea], [-(y1 - y0) * x for x in eb]


def combo_max(tp, tq, lam, mu, M):
    """Upper bound (exact) of lam*c_p + mu*c_q over the box."""
    r = [[lam * x + mu * z for x, z in zip(tp[k], tq[k])] for k in range(3)]
    return bound(r[0], M, True) + max(0, bound(r[1], M, True)) + max(0, bound(r[2], M, True))


def grad(tp, t):
    e0, es, ey = tp
    return (sum(es[k] * t ** k for k in range(5)), sum(ey[k] * t ** k for k in range(5)),
            sum(k * e0[k] * t ** (k - 1) for k in range(1, 5)))


def conflict(tp, tq, M, t):
    """True if a pair of positive multipliers shows that c_p > 0 and c_q > 0 cannot hold at the same time in the box."""
    if combo_max(tp, tq, 1, 1, M) <= 0:
        return True
    gp, gq = grad(tp, t), grad(tq, t)
    pp = sum(x * x for x in gp); pq = sum(x * y for x, y in zip(gp, gq))
    if pp == 0 or pq >= 0:
        return False
    lam = F(-pq / pp).limit_denominator(1 << 20)  # the multiplier may be chosen in floating point (the check is exact)
    return lam > 0 and combo_max(tp, tq, lam.numerator, lam.denominator, M) <= 0


def mwis_le(w, adj, budget):
    """True if every set of points that can be outside at the same time (a conflict-free set) has weight <= budget.
    w is the list of weights, adj[i] is the bit set of points in conflict with i. Returns False as soon as a heavy set is found."""
    n = len(w)
    order = sorted(range(n), key=lambda i: -w[i])
    def rec(k, cand, cur):
        if cur > budget:
            return False
        rest = sum(w[order[j]] for j in range(k, n) if cand >> order[j] & 1)
        if cur + rest <= budget:
            return True
        while k < n and not (cand >> order[k] & 1):
            k += 1
        if k == n:
            return True
        v = order[k]
        return rec(k + 1, cand & ~adj[v] & ~(1 << v), cur + w[v]) and rec(k + 1, cand & ~(1 << v), cur)
    return rec(0, (1 << n) - 1, 0)


def stage2(nu, pos, box, M, H, budget):
    """Use conflicts between undecided points (two points that cannot both be outside); ordinary boxes only.
    True if the weight that can be outside at the same time is shown to be <= budget (= decided + undecided - T),
    i.e. the capture is >= T for every pose in the box."""
    t = box[8] / (1 << K)
    cons = [[terms(pos[i][0], pos[i][1], c, box, H) for c in viol_cons(pos[i][0], pos[i][1], box, M, H)] for i in nu]
    n = len(nu); adj = [0] * n
    for a_ in range(n):
        for b_ in range(a_ + 1, n):
            if all(conflict(tp, tq, M, t) for tp in cons[a_] for tq in cons[b_]):
                adj[a_] |= 1 << b_; adj[b_] |= 1 << a_
    return mwis_le([pos[i][2] for i in nu], adj, budget)


def maxflow(n, edges, src, snk):
    """Dinic's maximum flow (integer capacities). edges are (u, v, capacity)."""
    g = [[] for _ in range(n)]
    for u, v, c in edges:
        g[u].append([v, c, len(g[v])]); g[v].append([u, 0, len(g[u]) - 1])
    flow = 0
    while True:
        lev = [-1] * n; lev[src] = 0; q = [src]
        for u in q:
            for v, c, _ in g[u]:
                if c > 0 and lev[v] < 0:
                    lev[v] = lev[u] + 1; q.append(v)
        if lev[snk] < 0:
            return flow
        it = [0] * n
        def dfs(u, f):
            if u == snk:
                return f
            while it[u] < len(g[u]):
                e = g[u][it[u]]
                if e[1] > 0 and lev[e[0]] == lev[u] + 1:
                    d = dfs(e[0], min(f, e[1]))
                    if d > 0:
                        e[1] -= d; g[e[0]][e[2]][1] += d
                        return d
                it[u] += 1
            return 0
        while True:
            f = dfs(src, 1 << 200)
            if f == 0:
                break
            flow += f


def mwis_bip(wa, wb, conf):
    """Weight of a maximum-weight independent set of a bipartite graph (left wa, right wb, conf = conflicting (i, j) pairs)
    = total weight - minimum vertex cover = total weight - maximum flow."""
    na, nb = len(wa), len(wb)
    if na == 0 or nb == 0 or not conf:
        return sum(wa) + sum(wb)
    src, snk = na + nb, na + nb + 1
    big = sum(wa) + sum(wb) + 1
    edges = [(src, i, w) for i, w in enumerate(wa)] + [(na + j, snk, w) for j, w in enumerate(wb)]
    edges += [(i, na + j, big) for i, j in conf]
    return sum(wa) + sum(wb) - maxflow(na + nb + 2, edges, src, snk)


def stage_bip(nu, pos, box, M, H, budget):
    """Use conflicts between opposite sides (bottom c4 with top c3, left c2 with right c1).
    At the same pose, adding c4(p) + c3(q) cancels s and sy and leaves a polynomial in t and the difference (dx, dy) of the points only:
      c4(p) + c3(q) = (dy-2H) - 2dx t - 4H t^2 - 2dx t^3 + (-dy-2H) t^4      (d = q - p, scaled by 2^E)
      c2(p) + c1(q) = (dx-2H) + 2dy t - 4H t^2 + 2dy t^3 + (-dx-2H) t^4
    (s and sy may depend on t; it works unchanged in wedge-coordinate boxes.)
    If this is <= 0 over the t-range of the box, p cannot leave through the bottom (left) while q leaves through the top (right).
    For each point, make one copy per constraint that may be violated; left = copies of bottom and left sides, right = copies
    of top and right sides. Edges are "conflict between a bottom copy and a top copy", "conflict between a left and a right copy",
    and "the left and right copies of the same point". If for each outside point we pick one violated constraint, the picked
    copies form an independent set, so the weight outside <= the MWIS of this bipartite graph.
    (Two copies on the same side (e.g. bottom and left) cannot be joined by an edge, so those are over-counted. This is sound as a lower bound.)
    If this is <= budget (= decided + undecided - T), the capture is >= T for every pose in the box."""
    left, right = [], []  # (point index, constraint, X, Y, weight)
    for i in nu:
        X, Y, w = pos[i]
        for c in viol_cons(X, Y, box, M, H):
            (left if c in (1, 3) else right).append((i, c, X << E, Y << E, w))
    if not left or not right:
        return sum(x[4] for x in left + right) <= budget
    cache = {}
    def ok(key, e):
        if key not in cache:
            cache[key] = bound(e, M, True) <= 0
        return cache[key]
    conf = []
    for ia, (pa, ca, xa, ya, _) in enumerate(left):
        for ib, (pb, cb, xb, yb, _) in enumerate(right):
            if pa == pb:
                conf.append((ia, ib)); continue
            dx, dy = xb - xa, yb - ya
            if ca == 3 and cb == 2:    # bottom c4 and top c3
                if ok((0, dx, dy), (dy - 2 * H, -2 * dx, -4 * H, -2 * dx, -dy - 2 * H)):
                    conf.append((ia, ib))
            elif ca == 1 and cb == 0:  # left c2 and right c1
                if ok((1, dx, dy), (dx - 2 * H, 2 * dy, -4 * H, 2 * dy, -dx - 2 * H)):
                    conf.append((ia, ib))
    return mwis_bip([x[4] for x in left], [x[4] for x in right], conf) <= budget


def apex_pieces(box, smax_e):
    """Re-split a small box whose lower end of t is 0 and whose ends of s, sy are integers into wedge-coordinate boxes.
    For an integer end e of the interval [s0, s1] of s (same for sy):
      e = s0:  s = e + t*u (u in [0, U])  and  s = S0 + U*t (S0 in [s0, s1])
      e = s1:  s = e - t*u (u in [0, U])  and  s = S0 - U*t (S0 in [s0, s1])
    For every t in [0, m1], the union of the two contains [s0, s1] (the part that sticks out is poses inside the container;
    checking it does not affect soundness). If there is no integer end, no re-split is done. Returns None if nothing was re-split."""
    s0, s1, u0, u1, y0, y1, v0, v1, m0, m1 = box
    one = 1 << E
    def split(a0, a1):
        if a0 % one == 0:
            return [(a0, a0, 0, U), (a0, a1, U, U)]
        if a1 % one == 0:
            return [(a1, a1, -U, 0), (a0, a1, -U, -U)]
        return None
    ps, py = split(s0, s1), split(y0, y1)
    if ps is None and py is None:
        return None
    ps = ps or [(s0, s1, 0, 0)]; py = py or [(y0, y1, 0, 0)]
    out = []
    for a in ps:
        for b in py:
            # check that the poses stay inside the container (s, sy in [0, smax])
            for lo, hi, w0, w1 in (a, b):
                assert lo + min(w0, 0) * m1 // (1 << K) - 1 >= 0 or min(w0, 0) == 0
                assert hi + max(w1, 0) * m1 // (1 << K) + 1 <= smax_e
            out.append((a[0], a[1], a[2], a[3], b[0], b[1], b[2], b[3], m0, m1))
    return out


def axis_check(L, cd, wd, pts, T):
    """Exact enumeration of the minimum capture over poses with theta = 0 (axis-parallel closed unit squares
    [s, s+cd] x [sy, sy+cd], s, sy in [0, L-cd]).
    Membership in x changes only at s = X - cd (enters from this value) and s = X (leaves just after this value), so the
    capture is constant on each pair of {break value, open interval between neighbouring break values} in the two axes.
    Combining all break values and their midpoints gives the minimum. Coordinates are doubled so they stay integers.
    The whole container is examined (the D4 reduction is not used)."""
    import numpy as np
    from bisect import bisect_left, bisect_right
    d = {}
    for x, y, w in pts:
        d[(x, y)] = d.get((x, y), 0) + w
    tot = sum(d.values())
    need(tot < 1 << 62, "the total weight does not fit in int64")
    xs = sorted({x for x, _ in d}); ys = sorted({y for _, y in d})
    hi = L - cd

    def cands(vals):
        cut = sorted({0, hi} | {v for v in vals if 0 <= v <= hi} | {v - cd for v in vals if 0 <= v - cd <= hi})
        return sorted([2 * v for v in cut] + [cut[i] + cut[i + 1] for i in range(len(cut) - 1)])

    def ranges(vals, cs):
        v2 = [2 * v for v in vals]
        lo = [bisect_left(v2, S) for S in cs]; up = [bisect_right(v2, S + 2 * cd) for S in cs]
        return np.array(lo), np.array(up)

    cx, cy = cands(xs), cands(ys)
    lx, hx = ranges(xs, cx); ly, hy = ranges(ys, cy)
    ix = {v: i for i, v in enumerate(xs)}; iy = {v: i for i, v in enumerate(ys)}
    P = np.zeros((len(xs) + 1, len(ys) + 1), dtype=np.int64)
    for (x, y), w in d.items():
        P[ix[x] + 1, iy[y] + 1] += w
    P = P.cumsum(0).cumsum(1)
    need(int(P[-1, -1]) == tot, "the prefix sums do not add up")
    best = None
    for a in range(0, len(cx), 256):
        A, B = hx[a:a + 256], lx[a:a + 256]
        cap = P[np.ix_(A, hy)] - P[np.ix_(B, hy)] - P[np.ix_(A, ly)] + P[np.ix_(B, ly)]
        k = int(cap.argmin()); v = int(cap.flat[k])
        if best is None or v < best[0]:
            best = (v, cx[a + k // len(cy)], cy[k % len(cy)])
    v, S, Sy = best
    # recount the capture of the minimising pose with Python integers
    chk = sum(w for (x, y), w in d.items() if S <= 2 * x <= S + 2 * cd and Sy <= 2 * y <= Sy + 2 * cd)
    need(chk == v, "the recount of the minimum capture does not match")
    return dict(ok=v >= T * wd, min=F(v, wd), s=F(S, 2 * cd), sy=F(Sy, 2 * cd), ncand=(len(cx), len(cy)))


def run(L, cd, wd, pts, T, minutes, tmin=F(0), max_fail=50, log=print, sub=None, use2=True, usebip=True, max2=60, maxbip=400,
        apex=True):
    H = (cd // 2) << E
    need(cd % 2 == 0, "the coordinate denominator must be even (so that the half side H is an integer)")
    half = F(L, 2) - F(cd, 2)
    smax = 1
    while smax < half: smax *= 2
    assert smax + 1.41422 * cd <= L, "the root box contains poses that stick out of the container (L/2 - half side cannot be covered by a power of 2)"
    shalf = int(half * (1 << E))
    tmax = 7 << (K - 4)
    tlo = int(tmin * (1 << K)); need(F(tlo, 1 << K) == tmin and 0 <= tmin <= F(7, 16), "tmin must be a multiple of 2^-K in [0, 7/16]")
    Tw = T * wd  # threshold in weight units (rational)
    root = (0, smax << E, 0, 0, 0, smax << E, 0, 0, tlo, tmax)
    if sub is not None:  # check only part of the range (s, sy in coordinate units, t as is; multiples of 2^-E and 2^-K)
        q = [sub[0] * (1 << E), sub[1] * (1 << E), sub[2] * (1 << E), sub[3] * (1 << E), sub[4] * (1 << K), sub[5] * (1 << K)]
        need(all(F(v).denominator == 1 for v in q), "--sub needs s, sy to be multiples of 2^-E and t a multiple of 2^-K")
        q = [int(v) for v in q]
        need(q[0] <= q[1] and q[2] <= q[3] and 0 <= q[4] <= q[5], "--sub has an interval in the wrong order, or t < 0")
        root = (q[0], q[1], 0, 0, q[2], q[3], 0, 0, q[4], q[5])
    pos = [p for p in pts if p[2] > 0]
    stack = [(root, 0, list(range(len(pos))), 0, F(1))]  # (box, decided weight, undecided points, depth, volume fraction)
    t_start = time.time(); nbox = 0; acc_vol = F(0); fails = []; viol = []; maxdepth = 0; acc = 0; acc2 = 0; acc3 = 0; napex = 0
    last = t_start
    apex_w = 1 << (E - 6)  # size of boxes that are re-split into wedge coordinates (1/64 in coordinate units)
    while stack:
        if time.time() - t_start > minutes * 60:
            log("time limit reached"); break
        box, sure, und, depth, fv = stack.pop()
        s0, s1, u0, u1, y0, y1, v0, v1, m0, m1 = box
        pl = plain(box)
        if pl and (s0 > shalf or y0 > shalf):  # boxes with only cx > L/2 or cy > L/2 are unnecessary by symmetry
            acc_vol += fv; continue
        nbox += 1; maxdepth = max(maxdepth, depth)
        M = basis(m0, m1 - m0)
        nu = []
        for i in und:
            X, Y, w = pos[i]
            c = classify(X, Y, box, M, H)
            if c == 1: sure += w
            elif c == 0: nu.append(i)
        ub = sure + sum(pos[i][2] for i in nu)
        if sure >= Tw:
            acc += 1; acc_vol += fv; continue
        if ub < Tw:
            viol.append((box, F(ub, wd)))
            if len(viol) >= max_fail: log("reached the limit on violation boxes"); break
            continue
        if usebip and len(nu) <= maxbip and stage_bip(nu, pos, box, M, H, ub - Tw):
            acc += 1; acc3 += 1; acc_vol += fv; continue
        if use2 and pl and len(nu) <= max2 and stage2(nu, pos, box, M, H, ub - Tw):
            acc += 1; acc2 += 1; acc_vol += fv; continue
        # Split: s, sy, t are compared by how far a side moves (coordinate units). For the wedge coordinates u, v, the test
        # is decided by the sign of expressions like t*(v + 2a) (dividing by t gives the same scale as a),
        # so they are compared by their widths in coordinate units as they are.
        dt = (m1 - m0) / (1 << K)
        ws, wy = (s1 - s0) / (1 << E), (y1 - y0) / (1 << E)
        wu, wv = (u1 - u0) / (1 << E), (v1 - v0) / (1 << E)
        wt = dt * (1.4 * cd + max(abs(u0), abs(u1), abs(v0), abs(v1)) / (1 << E))
        if apex and pl and m0 == 0 and max(ws, wy) * (1 << E) <= apex_w and wt * (1 << E) <= apex_w:
            kids = apex_pieces(box, smax << E)
            if kids is not None:
                # Wedge-coordinate boxes are not subsets of the original box (s = e + U*t + r etc. stick out a little),
                # so "surely inside" and "surely outside" decided in the original box are not inherited: every point is re-classified.
                napex += 1
                for k in kids:
                    stack.append((k, 0, list(range(len(pos))), depth + 1, fv / len(kids)))
                continue
        if max(ws, wy, wu, wv) < 2 ** -(E - 2) and wt < 2 ** -(K - 30):
            fails.append((box, F(sure, wd), F(ub, wd), [pos[i] for i in nu]))
            if len(fails) >= max_fail: log("reached the limit on split-limit boxes"); break
            continue
        w = max(ws, wy, wu, wv, wt)
        if w == ws:
            mid = (s0 + s1) // 2; kids = [(s0, mid) + box[2:], (mid, s1) + box[2:]]
        elif w == wy:
            mid = (y0 + y1) // 2; kids = [box[:4] + (y0, mid) + box[6:], box[:4] + (mid, y1) + box[6:]]
        elif w == wu:
            mid = (u0 + u1) // 2; kids = [box[:2] + (u0, mid) + box[4:], box[:2] + (mid, u1) + box[4:]]
        elif w == wv:
            mid = (v0 + v1) // 2; kids = [box[:6] + (v0, mid) + box[8:], box[:6] + (mid, v1) + box[8:]]
        else:
            mid = (m0 + m1) // 2; kids = [box[:8] + (m0, mid), box[:8] + (mid, m1)]
        for k in kids:
            stack.append((k, sure, nu, depth + 1, fv / 2))
        if time.time() - last > 60:
            last = time.time()
            log(f"  {last - t_start:6.0f}s boxes {nbox} accepted {acc} (opposite sides {acc3}, conflicts {acc2}) wedge {napex} accepted volume {float(acc_vol):.4%} violations {len(viol)} limit {len(fails)} depth {maxdepth} remaining {len(stack)}")
    done = not stack
    # Boxes whose lower end of t is 0 are classified excluding t = 0 (classify). If the root's lower end of t is 0, axis_check must also pass.
    return dict(nbox=nbox, acc=acc, acc2=acc2, acc3=acc3, napex=napex, acc_frac=acc_vol, viol=viol, fails=fails, maxdepth=maxdepth,
                done=done, sec=time.time() - t_start, smax=smax, need_axis=root[8] == 0)


def fmt_box(box, cd):
    s0, s1, u0, u1, y0, y1, v0, v1, m0, m1 = box
    f = lambda v: float(F(v, 1 << E) / cd)
    g = lambda m: float(F(m, 1 << K))
    w = f"s∈[{f(s0):.6f},{f(s1):.6f}] sy∈[{f(y0):.6f},{f(y1):.6f}] t∈[{g(m0):.6g},{g(m1):.6g}]"
    if not plain(box):
        w += f" u∈[{f(u0):.4g},{f(u1):.4g}] v∈[{f(v0):.4g},{f(v1):.4g}]"
    return w


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("path"); ap.add_argument("--T", default="1")
    ap.add_argument("--minutes", type=float, default=25)
    ap.add_argument("--tmin-num", type=int, default=0); ap.add_argument("--tmin-den", type=int, default=1)
    ap.add_argument("--max-fail", type=int, default=50)
    ap.add_argument("--no-stage2", action="store_true", help="do not use the conflict stage (intervals only)")
    ap.add_argument("--no-bip", action="store_true", help="do not use the opposite-sides stage")
    ap.add_argument("--no-apex", action="store_true", help="do not re-split into wedge coordinates")
    ap.add_argument("--max2", type=int, default=60, help="upper limit on the number of undecided points for the conflict stage")
    ap.add_argument("--maxbip", type=int, default=400, help="upper limit on the number of undecided points for the opposite-sides stage")
    ap.add_argument("--dump", default=None, help="write violation / limit boxes to a JSON file")
    ap.add_argument("--sub", nargs=6, default=None, help="s0 s1 sy0 sy1 t0 t1 (rationals, e.g. 1499/1 1501/1 ... 0 1/1024)")
    a = ap.parse_args()
    L, cd, wd, pts = load(a.path)
    T = F(a.T)
    dup, sym, mass = check_data(L, cd, wd, pts)
    print(f"L={L}/{cd} points {len(pts)} total mass {float(mass):.10f} ({mass}) duplicates {dup} D4-invariant {sym} threshold T={T}")
    if not sym:
        print("not D4-invariant, so this checker's reduction to 1/8 cannot be used"); sys.exit(1)
    sub = [F(x) for x in a.sub] if a.sub else None
    r = run(L, cd, wd, pts, T, a.minutes, F(a.tmin_num, a.tmin_den), a.max_fail, sub=sub, use2=not a.no_stage2,
            usebip=not a.no_bip, max2=a.max2, maxbip=a.maxbip, apex=not a.no_apex)
    print(f"boxes {r['nbox']} accepted {r['acc']} (of which opposite sides {r['acc3']}, conflicts {r['acc2']}) wedge {r['napex']} "
          f"accepted volume {float(r['acc_frac']):.6%} max depth {r['maxdepth']} {r['sec']:.1f}s finished {r['done']}")
    print(f"violations (capture < T for every pose in the box) {len(r['viol'])}, split-limit boxes {len(r['fails'])}")
    for b, ub in r["viol"][:10]: print("  violation", fmt_box(b, cd), "upper bound", float(ub))
    for b, lb, ub, u in r["fails"][:10]: print("  limit", fmt_box(b, cd), "lower bound", float(lb), "upper bound", float(ub), "undecided", len(u))
    if a.dump:
        import json
        json.dump({"viol": [[list(b), str(ub)] for b, ub in r["viol"]],
                   "fails": [[list(b), str(lb), str(ub), u] for b, lb, ub, u in r["fails"]], "E": E, "K": K},
                  open(a.dump, "w"), indent=0)
    n_needed = mass / T
    ax_ok = True
    if r["need_axis"]:
        ax = axis_check(L, cd, wd, pts, T); ax_ok = ax["ok"]
        print(f"poses with theta = 0 (whole container, candidates {ax['ncand'][0]}x{ax['ncand'][1]}): minimum capture {float(ax['min']):.10f}"
              f" (s = {float(ax['s'])}, sy = {float(ax['sy'])}) {'>= T' if ax_ok else '< T (violation)'}")
    if not ax_ok:
        print("Conclusion: capture < T at a pose with theta = 0 (violation)")
    elif r["done"] and not r["viol"] and not r["fails"] and (sub is not None or a.tmin_num != 0):
        print(f"Conclusion: capture >= {T} for every pose in the given range (--sub / --tmin) (exact). Not a proof for the whole region")
    elif r["done"] and not r["viol"] and not r["fails"]:
        print(f"Conclusion: capture >= {T} for every pose (exact). Total mass/T = {float(n_needed):.8f}, so "
              f"more than {float(n_needed):.8f} unit squares do not fit in a container of side less than {L}/{cd}")
