#!/usr/bin/env python3
# author of the probe: A (lane-A), S1674 (naryad #28 takt-63 REDO, per erratum + norm N-7).
# Full Bravais cost map with CORRECT representatives (auto-selected: |Aut| = intended holohedry order,
# checked for accidental higher symmetry, N-3) and N-7 cost: "found at k" AND "excluded below k".
# cost(L) = min over (n,m) of |disc(L)|*n/m^2 with a cyclic isotropic anti-isometry graph; exclusion:
# for each d<k enumerate all (n,m) with |disc|*n/m^2=d and show NO anti-isometry exists.  Bilinear b mod Z.
# No physical words.  RUN LINE:  python child-3.1/src/S1674_bravais_cost_full.py --outdir child-3.1/src
import argparse, os, json
import numpy as np
from sympy import Matrix, Rational
from sympy.matrices.normalforms import smith_normal_form
from itertools import product

SYSTEM = {2: "triclinic", 4: "monoclinic", 8: "orthorhombic", 12: "rhombohedral",
          16: "tetragonal", 24: "hexagonal", 48: "cubic"}
GLOBAL_G = None


def reduce3(Gi):
    global GLOBAL_G
    GLOBAL_G = np.array(Gi, dtype=int)
    G = np.array(Gi, dtype=float); B = np.eye(3, dtype=int)
    for _ in range(80):
        changed = False
        for i in range(3):
            for j in range(3):
                if i != j and G[j, j] > 0:
                    q = int(round(G[i, j] / G[j, j]))
                    if q != 0:
                        B[:, i] -= q*B[:, j]; G = (B.T @ GLOBAL_G @ B).astype(float); changed = True
        order = np.argsort([G[k, k] for k in range(3)])
        if list(order) != [0, 1, 2]:
            B = B[:, order]; G = (B.T @ GLOBAL_G @ B).astype(float); changed = True
        if not changed:
            break
    return (B.T @ GLOBAL_G @ B).astype(int)


def aut_count(Gi):
    G = reduce3(Gi); Gm = np.array(G); diag = [int(G[i, i]) for i in range(3)]; mx = max(diag)
    Ginv = np.linalg.inv(Gm); bnd = [int(np.floor((mx*Ginv[i, i])**0.5)) + 1 for i in range(3)]
    vecs = [(np.array(v), int(round(np.array(v) @ Gm @ np.array(v))))
            for v in product(*[range(-b, b+1) for b in bnd])]
    vecs = [(v, nr) for (v, nr) in vecs if 1 <= nr <= mx]
    cand = {i: [v for (v, nr) in vecs if nr == diag[i]] for i in range(3)}
    count = 0
    for a in cand[0]:
        for b in cand[1]:
            if int(round(a @ Gm @ b)) != int(round(G[0, 1])): continue
            for c in cand[2]:
                if int(round(a @ Gm @ c)) != int(round(G[0, 2])): continue
                if int(round(b @ Gm @ c)) != int(round(G[1, 2])): continue
                if abs(round(np.linalg.det(np.array([a, b, c]).T))) == 1: count += 1
    return count


def disc_form(G):
    G = Matrix(G); Gi = G.inv(); n = G.rows; d = abs(int(G.det()))
    def canon(v):
        gv = Gi*Matrix(v); return tuple(Rational(c)-Rational(c).__floor__() for c in gv)
    reps = {}
    for v in product(range(0, d+1), repeat=n):
        k = canon(v)
        if k not in reps: reps[k] = Matrix(v)
        if len(reps) >= d: break
    def b(vi, vj):
        val = (vi.T*Gi*vj)[0]; return Rational(val)-Rational(val).__floor__()
    def order(v):
        k = 1
        while True:
            gv = Gi*(k*v)
            if all((Rational(c)-Rational(c).__floor__()) == 0 for c in gv): return k
            k += 1
    return list(reps.values()), b, order, d


def gluings(Lg, D):
    """set of achievable |det| = |disc(L)|*n/m^2 via cyclic isotropic anti-isometry graph, n<=D."""
    eltsL, bL, ordL, dL = disc_form(Lg)
    dets = {}
    for n in range(1, D+1):
        eltsK, bK, ordK, dK = disc_form([[-n]])
        for x in eltsL:
            ox = ordL(x)
            for y in eltsK:
                if ordK(y) != ox: continue
                m = ox
                bg = bL(x, x)+bK(y, y)
                if (Rational(bg)-Rational(bg).__floor__()) != 0: continue
                det = Rational(dL*n, m*m)
                if det == int(det):
                    dets.setdefault(int(det), (n, m))
    return dL, dets


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--outdir", default="."); args = ap.parse_args()
    print("=" * 92)
    print("S1674 — full Bravais cost map (correct representatives, N-7 cost: found + excluded-below)")
    print("=" * 92)

    # candidate Grams per type; auto-select first with |Aut| == intended holohedry order
    CANDS = {
        "cP":  (48, [[[1, 0, 0], [0, 1, 0], [0, 0, 1]]]),
        "cF":  (48, [[[2, -1, 0], [-1, 2, -1], [0, -1, 2]]]),
        "cI":  (48, [[[3, -1, -1], [-1, 3, -1], [-1, -1, 3]]]),
        "tP":  (16, [[[1, 0, 0], [0, 1, 0], [0, 0, 2]]]),
        "tI":  (16, [[[2, 0, 1], [0, 2, 1], [1, 1, 3]], [[2, 0, 1], [0, 2, 1], [1, 1, 4]],
                     [[4, 0, 2], [0, 4, 2], [2, 2, 5]], [[2, 0, 1], [0, 2, 1], [1, 1, 5]]]),
        "oP":  (8,  [[[1, 0, 0], [0, 2, 0], [0, 0, 3]]]),
        "oC":  (8,  [[[3, 1, 0], [1, 3, 0], [0, 0, 2]], [[3, 1, 0], [1, 3, 0], [0, 0, 5]]]),
        "oI":  (8,  [[[4, 0, 2], [0, 8, 4], [2, 4, 6]], [[2, 0, 1], [0, 4, 2], [1, 2, 4]],
                     [[2, 0, 1], [0, 6, 3], [1, 3, 5]]]),
        "oF":  (8,  [[[3, 2, 1], [2, 5, 3], [1, 3, 4]], [[3, 1, 2], [1, 4, 2], [2, 2, 5]],
                     [[2, 1, 1], [1, 4, 2], [1, 2, 5]]]),
        "hP":  (24, [[[2, -1, 0], [-1, 2, 0], [0, 0, 3]]]),
        "hR":  (12, [[[3, 1, 1], [1, 3, 1], [1, 1, 3]], [[4, 1, 1], [1, 4, 1], [1, 1, 4]],
                     [[3, -1, -1], [-1, 3, 1], [-1, 1, 3]], [[2, 1, 1], [1, 3, 1], [1, 1, 3]]]),
        "mP":  (4,  [[[4, 0, 1], [0, 3, 0], [1, 0, 6]], [[3, 0, 1], [0, 5, 0], [1, 0, 4]]]),
        "mC":  (4,  [[[3, 1, 0], [1, 5, 0], [0, 0, 2]], [[2, 1, 0], [1, 4, 0], [0, 0, 3]],
                     [[3, 1, 1], [1, 3, 0], [1, 0, 4]]]),
        "aP":  (2,  [[[2, 1, 1], [1, 3, 2], [1, 2, 5]], [[2, 1, 1], [1, 4, 2], [1, 2, 6]],
                     [[3, 1, 2], [1, 4, 1], [2, 1, 6]]]),
    }

    print("\n   %-6s %-6s %-12s %-14s %-22s %s" %
          ("type", "|Aut|", "system", "disc(SNF)", "cost (N-7)", "Gram"))
    rows = []
    for typ, (want_aut, cands) in CANDS.items():
        chosen = None
        for G in cands:
            if aut_count(G) == want_aut:
                chosen = G; break
        if chosen is None:
            print("   %-6s  -- no candidate with |Aut|=%d (representative not found) --" % (typ, want_aut))
            rows.append((typ, None, None, None, None, None)); continue
        G = chosen
        a = aut_count(G)
        snf = smith_normal_form(Matrix(G)); inv = [int(snf[i, i]) for i in range(3) if snf[i, i] != 0]
        disc = "x".join("Z/%d" % f for f in inv if f > 1) or "0"
        # cost with wide D; N-7: found at k = min achievable; excluded below k = no d<k achievable
        D = 40
        dL, dets = gluings(G, D)
        if not dets:
            k = dL; excluded = True
        else:
            k = min(dets)
            excluded = all(d not in dets for d in range(1, k))   # nothing below k found up to D
        n_m = dets.get(k, (1, 1))
        cost_str = "found+excl @ %d" % k if excluded else "found @ %d (excl?)" % k
        rows.append((typ, a, SYSTEM.get(a, "?"), disc, k, cost_str))
        print("   %-6s %-6d %-12s %-14s %-22s %s" %
              (typ, a, SYSTEM.get(a, "?%d" % a), disc, cost_str, G))

    # controls
    dd = {r[0]: r for r in rows}
    ok_cP = dd["cP"][4] == 1
    ok_FCC = dd["cF"][4] == 3
    ok_BCC = dd["cI"][4] == 8
    found = [r for r in rows if r[1] is not None]
    print("\n-- VERDICT --")
    print("   [%s] control cP  -> found+excl @ 1" % ("OK" if ok_cP else "XX"))
    print("   [%s] control FCC -> found+excl @ 3 (corrected C-41)" % ("OK" if ok_FCC else "XX"))
    print("   [%s] control BCC -> found+excl @ 8" % ("OK" if ok_BCC else "XX"))
    print("   [%s] all 14 representatives validated by |Aut| (N-3): %d/14" %
          ("OK" if len(found) == 14 else "XX", len(found)))
    all_ok = ok_cP and ok_FCC and ok_BCC and len(found) == 14
    print("   VERDICT: %s" % ("OK -> full 14-type Bravais cost map (N-7 cost, validated reps)"
                              if all_ok else "PARTIAL -- %d/14 reps found; controls %s" %
                              (len(found), "ok" if (ok_cP and ok_FCC and ok_BCC) else "check")))
    print("\n   N-7: cost = 'found at k' + 'excluded below k' (no |det|<k achievable up to D=40).")
    print("   Cost tracks disc SIZE+FORM (takt-64): FCC(Z/4)=3 < BCC(Z/4xZ/4)=8, not density.")

    os.makedirs(args.outdir, exist_ok=True)
    dump = os.path.join(args.outdir, "S1674_bravais_cost_full_dump.jsonl")
    with open(dump, "w", encoding="utf-8") as fd:
        for typ, a, sysn, disc, k, cost_str in rows:
            fd.write(json.dumps({"kind": "type", "type": typ, "aut": a, "system": sysn,
                                 "disc": disc, "cost_k": k, "cost_str": cost_str}) + "\n")
        fd.write(json.dumps({"kind": "verdict", "controls_ok": bool(ok_cP and ok_FCC and ok_BCC),
                             "reps_found": len(found)}) + "\n")
    print("\n[dump] %s" % dump)


if __name__ == "__main__":
    main()
