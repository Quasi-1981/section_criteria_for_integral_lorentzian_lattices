#!/usr/bin/env python3
# author of the probe: B (lane-B chair a7296aa8), S1756.  Ex-ante: S1756_GENUS_PINNING_EXANTE.md
# QUESTION: does the triple (signature, bilinear disc form b-bar up to isometry, parity) determine
# the GENUS, for rank 4, sig (3,1), |det| in {1,2,4}?
#   KILL side is decisive: two lattices, same signature, b-bar EXPLICITLY isometric, same parity,
#   but different genus invariants  =>  (PIN) is false  =>  P7.Thm1 false as intended.
#   UNIQUENESS side is NOT closable here (full 2-adic canon parked at S1736; SPLAG (29)-(35) not
#   transcribed in the repo).  Absence of a counterexample is reported as "not killed, NOT proved".
# Genus invariants used (both congruence-invariant, both computed natively):
#   (a) partial 2-adic symbol per S1736: SNF 2-scales/ranks, v2(det), odd part of det mod 8
#   (b) Hasse invariant over Q_2 via Hilbert symbols (exact rational arithmetic)
# No physical words.  Handles: 0.
# RUN LINE:  python child-3.1/src/S1756_genus_pinning.py --outdir child-3.1/src
import argparse, os, json, itertools
from fractions import Fraction
from sympy import Matrix, Rational

# ----------------------------------------------------------------------------- basic lattice data

def signature(G):
    import numpy as np
    ev = np.linalg.eigvalsh(np.array(Matrix(G).tolist(), float))
    return (int((ev > 1e-9).sum()), int((ev < -1e-9).sum()))

def is_even(G):
    return all(G[i][i] % 2 == 0 for i in range(len(G)))

def smith_invariants(G):
    """elementary divisors of the Gram => structure of disc(L) = L^v/L"""
    M = Matrix(G)
    from sympy.matrices.normalforms import smith_normal_form
    S = smith_normal_form(M)
    d = [abs(int(S[i, i])) for i in range(S.rows)]
    return [x for x in d if x != 1], d

def disc_elements(G):
    """disc(L) = Z^n/G Z^n with b(x,y) = x^T G^{-1} y mod Z.  Enumerate via elementary divisors."""
    M = Matrix(G)
    Ginv = M.inv()
    n = M.rows
    from sympy.matrices.normalforms import smith_normal_form
    # explicit U,V with U*M*V = S is not exposed; enumerate the quotient by brute force on a box
    divs, full = smith_invariants(G)
    order = 1
    for d in full:
        order *= d
    # representatives: Z^n / M Z^n -- brute force over a box of size max(full)
    B = max(full) if full else 1
    seen = {}
    reps = []
    for v in itertools.product(range(B), repeat=n):
        key = canon_class(Matrix([[x] for x in v]), M, B)
        if key not in seen:
            seen[key] = v
            reps.append(v)
        if len(reps) == order:
            break
    return reps, Ginv, order

def canon_class(vec, M, B):
    """canonical key of vec in Z^n / M Z^n : reduce by lattice M Z^n using HNF of M"""
    from sympy.matrices.normalforms import hermite_normal_form
    H = hermite_normal_form(M)
    n = M.rows
    w = list(vec)
    # reduce greedily using columns of H (upper triangular-ish)
    for i in range(n - 1, -1, -1):
        piv = int(H[i, i])
        if piv == 0:
            continue
        q = w[i] // piv
        if q:
            for r in range(n):
                w[r] -= q * int(H[r, i])
    return tuple(int(x) for x in w)

def bform(x, y, Ginv):
    v = Matrix([[a] for a in x])
    u = Matrix([[a] for a in y])
    val = (v.T * Ginv * u)[0, 0]
    return Rational(val) % 1

# ----------------------------------------------------------------------------- disc form invariants

def disc_profile(G):
    """group elementary divisors + full multiset of b(x,x) + full b-table (for isometry search)"""
    reps, Ginv, order = disc_elements(G)
    diag = sorted([str(bform(x, x, Ginv)) for x in reps])
    divs, full = smith_invariants(G)
    return {"order": order, "elementary_divisors": tuple(divs),
            "diag_multiset": tuple(diag)}, reps, Ginv

def bbar_isometric(G1, G2):
    """EXHAUSTIVE: is there a group isomorphism disc(G1)->disc(G2) preserving b?  |disc| <= 4 here."""
    p1, reps1, Gi1 = disc_profile(G1)
    p2, reps2, Gi2 = disc_profile(G2)
    if p1["order"] != p2["order"] or p1["elementary_divisors"] != p2["elementary_divisors"]:
        return False, None
    n1, n2 = len(reps1[0]), len(reps2[0])
    M1, M2 = Matrix(G1), Matrix(G2)
    B1 = max(smith_invariants(G1)[1]); B2 = max(smith_invariants(G2)[1])
    idx1 = {canon_class(Matrix([[c] for c in r]), M1, B1): i for i, r in enumerate(reps1)}
    idx2 = {canon_class(Matrix([[c] for c in r]), M2, B2): i for i, r in enumerate(reps2)}

    def add1(a, b):
        return canon_class(Matrix([[a[i] + b[i]] for i in range(n1)]), M1, B1)

    def add2(a, b):
        return canon_class(Matrix([[a[i] + b[i]] for i in range(n2)]), M2, B2)

    # brute force over all bijections reps1 -> reps2 fixing 0, checking additivity + b-preservation
    zero1 = canon_class(Matrix([[0]] * n1), M1, B1)
    zero2 = canon_class(Matrix([[0]] * n2), M2, B2)
    nz1 = [r for r in reps1 if canon_class(Matrix([[c] for c in r]), M1, B1) != zero1]
    nz2 = [r for r in reps2 if canon_class(Matrix([[c] for c in r]), M2, B2) != zero2]
    if len(nz1) != len(nz2):
        return False, None
    for perm in itertools.permutations(nz2):
        phi = {zero1: zero2}
        for a, b in zip(nz1, perm):
            phi[canon_class(Matrix([[c] for c in a]), M1, B1)] = \
                canon_class(Matrix([[c] for c in b]), M2, B2)
        # additivity
        ok = True
        for a in reps1:
            for b in reps1:
                ka = canon_class(Matrix([[c] for c in a]), M1, B1)
                kb = canon_class(Matrix([[c] for c in b]), M1, B1)
                if phi[add1(a, b)] != add2(reps2[idx2[phi[ka]]], reps2[idx2[phi[kb]]]):
                    ok = False
                    break
            if not ok:
                break
        if not ok:
            continue
        # b-preservation
        for a in reps1:
            for b in reps1:
                ka = canon_class(Matrix([[c] for c in a]), M1, B1)
                kb = canon_class(Matrix([[c] for c in b]), M1, B1)
                if bform(a, b, Gi1) != bform(reps2[idx2[phi[ka]]], reps2[idx2[phi[kb]]], Gi2):
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return True, phi
    return False, None

# ----------------------------------------------------------------------------- genus invariants

def two_adic_partial(G):
    """S1736 partial 2-adic symbol: 2-scales/ranks from SNF, v2(det), odd part of det mod 8."""
    divs, full = smith_invariants(G)
    scales = {}
    for d in full:
        e = 0
        t = d
        while t % 2 == 0:
            t //= 2
            e += 1
        scales[e] = scales.get(e, 0) + 1
    det = int(Matrix(G).det())
    a = abs(det)
    v2 = 0
    while a % 2 == 0:
        a //= 2
        v2 += 1
    return {"scales": tuple(sorted(scales.items())), "v2_det": v2,
            "odd_det_mod8": a % 8, "det_sign": 1 if det > 0 else -1}

def hilbert_symbol_2(a, b):
    """(a,b)_2 for nonzero rationals, exact."""
    def v2_unit(x):
        x = Fraction(x)
        v = 0
        num, den = x.numerator, x.denominator
        while num % 2 == 0:
            num //= 2; v += 1
        while den % 2 == 0:
            den //= 2; v -= 1
        return v, Fraction(num, den)
    al, au = v2_unit(a)
    be, bu = v2_unit(b)
    def eps(u):
        u = Fraction(u)
        r = (u.numerator * pow(u.denominator, -1, 8)) % 8
        return ((r - 1) // 2) % 2
    def omega(u):
        u = Fraction(u)
        r = (u.numerator * pow(u.denominator, -1, 8)) % 8
        return ((r * r - 1) // 8) % 2
    e = (eps(au) * eps(bu) + al * omega(bu) + be * omega(au)) % 2
    return 1 if e == 0 else -1

def diagonalize_rational(G):
    """rational congruence diagonalisation (exact)."""
    M = Matrix(G).applyfunc(Rational)
    n = M.rows
    d = []
    M = M[:, :]
    idx = list(range(n))
    A = [[Rational(M[i, j]) for j in range(n)] for i in range(n)]
    size = n
    while size > 0:
        # find nonzero diagonal
        p = None
        for i in range(size):
            if A[i][i] != 0:
                p = i; break
        if p is None:
            found = False
            for i in range(size):
                for j in range(i + 1, size):
                    if A[i][j] != 0:
                        for r in range(size):
                            A[r][i] += A[r][j]
                        for c in range(size):
                            A[i][c] += A[j][c]
                        found = True; break
                if found: break
            if not found:
                d.extend([Rational(0)] * size); break
            continue
        A[0], A[p] = A[p], A[0]
        for r in range(size):
            A[r][0], A[r][p] = A[r][p], A[r][0]
        piv = A[0][0]
        d.append(piv)
        for i in range(1, size):
            f = A[i][0] / piv
            for j in range(size):
                A[i][j] -= f * A[0][j]
            for j in range(size):
                A[j][i] -= f * A[j][0]
        A = [row[1:size] for row in A[1:size]]
        size -= 1
    return d

def hasse_2(G):
    """Hasse invariant at p=2 of the rational quadratic form: prod_{i<j} (d_i,d_j)_2"""
    d = diagonalize_rational(G)
    h = 1
    for i in range(len(d)):
        for j in range(i + 1, len(d)):
            h *= hilbert_symbol_2(d[i], d[j])
    return h

def genus_invariants(G):
    inv = two_adic_partial(G)
    inv["hasse_2"] = hasse_2(G)
    inv["sig"] = signature(G)
    return inv

# ----------------------------------------------------------------------------- family

AMBIENTS = {
    "I_{3,1}=diag(1,1,1,-1)": [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-1]],
    "diag(1,1,1,-2)":         [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-2]],
    "diag(1,1,1,-4)":         [[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,-4]],
    "diag(1,1,2,-2)":         [[1,0,0,0],[0,1,0,0],[0,0,2,0],[0,0,0,-2]],
    "A3(+)<-1>":              [[2,-1,0,0],[-1,2,-1,0],[0,-1,2,0],[0,0,0,-1]],
    "A1^2(+)U":               [[2,0,0,0],[0,2,0,0],[0,0,0,1],[0,0,1,0]],
}

def gl4_copies(G, k=3):
    """random-ish unimodular changes of basis, deterministic (no Date/random)."""
    import copy
    out = []
    Us = [
        [[1,1,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]],
        [[1,0,0,0],[0,1,0,0],[1,0,1,0],[0,1,0,1]],
        [[1,0,1,1],[0,1,0,1],[0,0,1,0],[0,0,0,1]],
    ]
    for U in Us[:k]:
        Um = Matrix(U)
        out.append([[int(x) for x in row] for row in (Um.T * Matrix(G) * Um).tolist()])
    return out


def det4(G):
    """exact integer determinant of a 4x4 integer matrix (fast, no sympy)."""
    a, b, c, d = G[0]
    e, f, g, h = G[1]
    i, j, k, l = G[2]
    m, n, o, p = G[3]
    return (a * (f * (k * p - l * o) - g * (j * p - l * n) + h * (j * o - k * n))
            - b * (e * (k * p - l * o) - g * (i * p - l * m) + h * (i * o - k * m))
            + c * (e * (j * p - l * n) - f * (i * p - l * m) + h * (i * n - j * m))
            - d * (e * (j * o - k * n) - f * (i * o - k * m) + g * (i * n - j * m)))

def enumerate_family(B, want_dets=(1, 2, 4)):
    """integral symmetric 4x4, |g_ij| <= B, sig (3,1), |det| in want_dets."""
    found = []
    rng = range(-B, B + 1)
    # diagonal entries bounded separately to keep the search finite and cheap
    for a in range(0, B + 1):
        for b in range(0, B + 1):
            for c in range(0, B + 1):
                for e in range(-B, 1):
                    for off in itertools.product(rng, repeat=6):
                        g01, g02, g03, g12, g13, g23 = off
                        G = [[a, g01, g02, g03],
                             [g01, b, g12, g13],
                             [g02, g12, c, g23],
                             [g03, g13, g23, e]]
                        det = det4(G)
                        if det == 0 or abs(det) not in want_dets:
                            continue
                        if signature(G) != (3, 1):
                            continue
                        found.append(G)
    return found

# ----------------------------------------------------------------------------- main

def bucket_key(G, with_parity=True):
    prof, _, _ = disc_profile(G)
    k = [prof["order"], prof["elementary_divisors"], prof["diag_multiset"], signature(G)]
    if with_parity:
        k.append("even" if is_even(G) else "odd")
    return tuple(str(x) for x in k)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=".")
    ap.add_argument("-B", type=int, default=2)
    args = ap.parse_args()
    print("=" * 96)
    print("S1756 - does (signature, b-bar, parity) determine the GENUS?  rank 4, sig(3,1), |det| in {1,2,4}")
    print("=" * 96)

    named = dict(AMBIENTS)
    for nm, G in list(AMBIENTS.items()):
        for i, C in enumerate(gl4_copies(G)):
            named["%s ~copy%d" % (nm, i + 1)] = C

    print("\n-- named lattices: invariants --")
    print("   %-28s %-8s %-7s %-22s %-8s %s" % ("name", "|det|", "parity", "disc (elem.div.)", "hasse2", "2-scales"))
    recs = {}
    for nm, G in named.items():
        gi = genus_invariants(G)
        prof, _, _ = disc_profile(G)
        recs[nm] = {"G": G, "gi": gi, "prof": prof, "even": is_even(G)}
        print("   %-28s %-8d %-7s %-22s %-8d %s" %
              (nm, abs(int(Matrix(G).det())), "even" if is_even(G) else "odd",
               str(prof["elementary_divisors"]), gi["hasse_2"], str(gi["scales"])))

    # ---------------- controls
    print("\n-- CONTROL K-POS: GL4(Z) copies must carry IDENTICAL genus invariants --")
    kpos_ok = True
    for nm, G in AMBIENTS.items():
        base = recs[nm]["gi"]
        for i in range(1, 4):
            key = "%s ~copy%d" % (nm, i)
            if key not in recs:
                continue
            same = (recs[key]["gi"] == base) and (recs[key]["prof"] == recs[nm]["prof"]) \
                   and (recs[key]["even"] == recs[nm]["even"])
            if not same:
                kpos_ok = False
                print("   [XX] %-28s copy%d DIFFERS -> invariant is not an invariant" % (nm, i))
    print("   [%s] all copies match their base" % ("OK" if kpos_ok else "XX"))

    print("\n-- CONTROL K-NEG (critical): the KNOWN parity pair on the 4th ambient --")
    G_odd = AMBIENTS["diag(1,1,2,-2)"]
    G_even = AMBIENTS["A1^2(+)U"]
    iso, phi = bbar_isometric(G_odd, G_even)
    k_with = (bucket_key(G_odd, True) == bucket_key(G_even, True))
    k_without = (bucket_key(G_odd, False) == bucket_key(G_even, False))
    print("   b-bar explicitly isometric (exhaustive search): %s" % iso)
    print("   same bucket WITH parity in key   : %s   (expected False)" % k_with)
    print("   same bucket WITHOUT parity in key: %s   (expected True)" % k_without)
    kneg_ok = (iso is True) and (k_with is False) and (k_without is True)
    print("   [%s] K-NEG" % ("OK" if kneg_ok else "XX -> MEASUREMENT VOID"))

    print("\n-- CONTROL K-STRUCT: |disc| == |det|, and Z/4 vs (Z/2)^2 separate at |det|=4 --")
    kstruct_ok = all(recs[nm]["prof"]["order"] == abs(int(Matrix(recs[nm]["G"]).det())) for nm in recs)
    z4 = recs["diag(1,1,1,-4)"]["prof"]["elementary_divisors"]
    z22 = recs["diag(1,1,2,-2)"]["prof"]["elementary_divisors"]
    print("   |disc|==|det| everywhere: %s ;  Z/4 %s  vs  (Z/2)^2 %s  -> distinct: %s"
          % (kstruct_ok, z4, z22, z4 != z22))
    kstruct_ok = kstruct_ok and (z4 != z22)
    print("   [%s] K-STRUCT" % ("OK" if kstruct_ok else "XX"))

    print("\n-- CONTROL K-CAL: Hilbert symbol (.,.)_2 against known values --")
    KNOWN = [(-1, -1, -1), (2, 2, 1), (-1, 2, 1), (2, 5, -1), (5, 5, 1),
             (3, 3, -1), (1, 7, 1), (2, 3, -1), (-1, -7, 1)]
    kcal_ok = True
    for a, b, exp in KNOWN:
        got = hilbert_symbol_2(a, b)
        if got != exp:
            kcal_ok = False
            print("   [XX] (%d,%d)_2 = %+d, expected %+d" % (a, b, got, exp))
    h_neg = hasse_2([[3,0,0,0],[0,3,0,0],[0,0,3,0],[0,0,0,-3]])
    print("   known values: %d/%d ; hasse_2(diag(3,3,3,-3)) = %+d (must be -1, proves it CAN fire)"
          % (sum(1 for a, b, e in KNOWN if hilbert_symbol_2(a, b) == e), len(KNOWN), h_neg))
    kcal_ok = kcal_ok and (h_neg == -1)
    print("   [%s] K-CAL" % ("OK" if kcal_ok else "XX"))

    # ---------------- the hunt
    print("\n-- FAMILY ENUMERATION (B=%d) --" % args.B)
    fam = enumerate_family(args.B)
    print("   integral 4x4, |g_ij| <= %d, sig (3,1), |det| in {1,2,4}: %d lattices" % (args.B, len(fam)))
    for nm, G in named.items():
        fam.append(G)

    buckets = {}
    for G in fam:
        k = bucket_key(G, True)
        buckets.setdefault(k, []).append(G)
    print("   buckets by triple (sig, disc profile, parity): %d" % len(buckets))

    print("\n-- HUNT: inside each bucket, do all lattices share the genus invariants? --")
    counterexamples = []
    for k, Gs in buckets.items():
        base = genus_invariants(Gs[0])
        for G in Gs[1:]:
            gi = genus_invariants(G)
            if gi != base:
                # verify b-bar really isometric (exhaustive) before claiming anything
                iso, _ = bbar_isometric(Gs[0], G)
                if iso and is_even(Gs[0]) == is_even(G):
                    counterexamples.append({"bucket": k, "G1": Gs[0], "G2": G,
                                            "inv1": str(base), "inv2": str(gi)})
    print("   candidate counterexamples (same triple, different genus invariants): %d" % len(counterexamples))
    for c in counterexamples[:5]:
        print("      G1=%s" % c["G1"])
        print("      G2=%s" % c["G2"])
        print("      inv1=%s" % c["inv1"])
        print("      inv2=%s" % c["inv2"])

    # ---------------- K-POWER: can the instrument fire AT ALL inside a bucket?
    print("\n-- CONTROL K-POWER (added after the first run caught blindness) --")
    print("   'no counterexample' is worthless unless the genus invariants CAN differ inside a bucket.")
    hasse_seen = set()
    scale_of_group = {}
    scale_is_function_of_group = True
    for G in fam:
        hasse_seen.add(hasse_2(G))
        prof, _, _ = disc_profile(G)
        gkey = prof["elementary_divisors"]
        sval = two_adic_partial(G)["scales"]
        if gkey in scale_of_group and scale_of_group[gkey] != sval:
            scale_is_function_of_group = False
        scale_of_group[gkey] = sval
    split_buckets = 0
    for k, Gs in buckets.items():
        if len({genus_invariants(G)["hasse_2"] for G in Gs}) > 1:
            split_buckets += 1
    print("   distinct hasse_2 values over the whole family : %s" % sorted(hasse_seen))
    print("   2-adic scale profile is a FUNCTION of the disc group (adds nothing): %s"
          % scale_is_function_of_group)
    print("   buckets in which any invariant takes >1 value  : %d" % split_buckets)
    kpower_ok = (len(hasse_seen) > 1) or (not scale_is_function_of_group) or (split_buckets > 0)
    print("   [%s] K-POWER - instrument %s separate inside a bucket"
          % ("OK" if kpower_ok else "XX", "CAN" if kpower_ok else "CANNOT"))

    # ---------------- verdict
    print("\n-- VERDICT --")
    controls = kpos_ok and kneg_ok and kstruct_ok and kcal_ok
    print("   controls (K-POS, K-NEG, K-STRUCT, K-CAL): %s" % ("ALL OK" if controls else "FAILED"))
    if not controls:
        verdict = "VOID"
        print("   VERDICT: MEASUREMENT VOID - controls failed, no conclusion may be drawn.")
    elif not kpower_ok:
        verdict = "VOID-BLIND"
        print("   VERDICT: MEASUREMENT VOID - THE INSTRUMENT IS BLIND.")
        print("            Every genus invariant available to this probe is already a function of")
        print("            (disc group, parity, signature) on this family: hasse_2 is CONSTANT and the")
        print("            2-adic scale profile is determined by the elementary divisors.")
        print("            => a counterexample could exist and this probe could never see it.")
        print("            '0 counterexamples' here is VACUOUS, not evidence. The hunt DID NOT HAPPEN.")
        print("            To move the question one needs an invariant strictly finer than the disc")
        print("            group: the full 2-adic canonical symbol (oddity fusion / sign walking),")
        print("            parked at S1736, or SPLAG (29)-(35) / Table 15.4 transcribed from source.")
    elif counterexamples:
        verdict = "KILL"
        print("   VERDICT: KILL - (PIN) is FALSE: same triple, different genus.")
        print("            => P7.Thm1 is false as intended; the honest criterion needs finer data.")
    else:
        verdict = "NOT-KILLED"
        print("   VERDICT: NOT KILLED, and NOT PROVED.")
        print("            No counterexample inside the scope (B=%d)." % args.B)
        print("            This is 'none found', NOT 'none exists': the uniqueness half needs the")
        print("            full 2-adic canonical symbol (parked at S1736) or SPLAG (29)-(35)/Table 15.4")
        print("            transcribed from the source.  NOT a green light to print.")

    os.makedirs(args.outdir, exist_ok=True)
    dump = os.path.join(args.outdir, "S1756_genus_pinning_dump.jsonl")
    with open(dump, "w", encoding="utf-8") as fd:
        for nm in recs:
            fd.write(json.dumps({"kind": "lattice", "name": nm, "gram": recs[nm]["G"],
                                 "even": recs[nm]["even"],
                                 "elementary_divisors": list(recs[nm]["prof"]["elementary_divisors"]),
                                 "hasse_2": recs[nm]["gi"]["hasse_2"],
                                 "scales": [list(t) for t in recs[nm]["gi"]["scales"]]}) + "\n")
        fd.write(json.dumps({"kind": "controls", "k_pos": bool(kpos_ok), "k_neg": bool(kneg_ok),
                             "k_struct": bool(kstruct_ok)}) + "\n")
        fd.write(json.dumps({"kind": "verdict", "verdict": verdict, "B": args.B,
                             "family_size": len(fam), "buckets": len(buckets),
                             "counterexamples": len(counterexamples)}) + "\n")
    print("\n[dump] %s" % dump)


if __name__ == "__main__":
    main()
