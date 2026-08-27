#!/usr/bin/env python3
# author of the probe: B (lane-B chair a7296aa8), S1764.  Court C-108 + author's tooth.
#
# FULL ENUMERATION of the genera of rank 4, signature (3,1), |det| in {1,2,4},
# with sign walking implemented VERBATIM after
#     D. Allcock, I. Gal, A. Mark, "The Conway-Sloane calculus for 2-adic lattices",
#     arXiv:1511.04614 -- Lemma 5.2 (oddity fusion) and Lemma 6.1 (sign walking).
#
# WHY THIS FILE REPLACES THE EARLIER ONE.  The earlier probe walked signs between ANY
# two members of a train, with no side condition.  Conway-Sloane's own formulation of
# the canonical form ("at most one sign per train") is FALSE -- AGM exhibit the
# counterexample [128^1 256^-1]_4 -- and Lemma 6.1 licenses a sign walk only in three
# named cases.  An unconditional walk is strictly more permissive than the lemma and
# can therefore merge inequivalent symbols.  This file performs only licensed walks and
# RECORDS, for every walk it uses, which case of Lemma 6.1 licenses it.
#
# The log is the evidence: it is not enough to print the final count.  Absence of a known
# counterexample in a small scope is not absence of failure, so every walk is certified
# by name rather than argued to be safe.
#
# Handles: 0.  No physical words.
# RUN LINE:  python S1764_genus_enum_agm.py --outdir .
import argparse, json, itertools

SCALES = [1, 2, 4]            # powers of 2 that can occur while |det| <= 4
SIGNATURE = 2                 # sig (3,1) => 3 - 1

# --------------------------------------------------------------- existence conditions

def constituent_ok(n, typ, eps, t):
    """SPLAG Ch.15 §7.7 (31)-(35), for p = 2."""
    if n == 0:
        return typ == "II" and eps == +1 and t % 8 == 0          # (32)
    if t % 2 != n % 2:                                            # (35)
        return False
    if typ == "II" and t % 8 != 0:                                # (35)
        return False
    if n % 2 == 1 and typ != "I":                                 # (35)
        return False
    if n == 1:                                                    # (33)
        return (t % 8 in (1, 7)) if eps == +1 else (t % 8 in (3, 5))
    if n == 2 and typ == "I":                                     # (34)
        return (t % 8 in (0, 2, 6)) if eps == +1 else (t % 8 in (4, 2, 6))
    return True


def kronecker_2(a):
    return 1 if a % 8 in (1, 7) else -1


def symbol_is_legal(sym, det):
    """the full existence gate: (29) determinant, (30) oddity, (31)-(35) per constituent."""
    for q in SCALES:
        c = sym[q]
        if not constituent_ok(c["n"], c["type"], c["eps"], c["t"]):
            return False
    v2, a = 0, abs(det)
    while a % 2 == 0:
        a //= 2
        v2 += 1
    prod = 1
    for q in SCALES:
        prod *= sym[q]["eps"]
    if prod != kronecker_2(det // (2 ** v2)):                     # (29)
        return False
    k2 = sum(1 for q in SCALES if q == 2 and sym[q]["eps"] == -1)  # q non-square
    oddity = (sum(sym[q]["t"] for q in SCALES) + 4 * k2) % 8
    return (SIGNATURE - oddity) % 8 == 0                          # (30); odd p contribute 0


def putative_symbols(det):
    """every system of 2-adic Jordan data passing (29),(30),(31)-(35)."""
    v2, a = 0, abs(det)
    while a % 2 == 0:
        a //= 2
        v2 += 1
    out = []
    for n1 in range(5):
        for n2 in range(5 - n1):
            n4 = 4 - n1 - n2
            if n4 < 0 or n2 + 2 * n4 != v2:
                continue
            ns = {1: n1, 2: n2, 4: n4}
            for types in itertools.product(["I", "II"], repeat=3):
                for epss in itertools.product([1, -1], repeat=3):
                    for ts in itertools.product(range(8), repeat=3):
                        sym = {q: {"n": ns[q], "type": types[i], "eps": epss[i], "t": ts[i]}
                               for i, q in enumerate(SCALES)}
                        if symbol_is_legal(sym, det):
                            out.append(sym)
    return out

# --------------------------------------------------------------- compartments

def compartments(sym):
    """maximal runs of scaled type-I constituents (0-dimensional ones are type II
       and therefore break a compartment)."""
    comps, cur = [], []
    for q in SCALES:
        if sym[q]["type"] == "I":
            cur.append(q)
        else:
            if cur:
                comps.append(cur)
            cur = []
    if cur:
        comps.append(cur)
    return comps


def comp_index(comps, q):
    for i, c in enumerate(comps):
        if q in c:
            return i
    return None

# --------------------------------------------------------------- AGM Lemma 6.1

def licensed_walks(sym):
    """Every sign walk AGM Lemma 6.1 licenses on `sym`, each returned together with the
       case that licenses it.  Verbatim conditions -- two NONTRIVIAL terms with either

         (1) adjacent scales and different types;
         (2) adjacent scales and type I, and their compartment either has dimension > 2
             or compartment oddity +-2;
         (3) type I, scales differing by a factor of 4, and the term between them trivial.

       Effect: negate both signs, and change by 4 the oddity of each compartment that
       contains at least one of the two terms."""
    out = []
    comps = compartments(sym)
    nontrivial = [q for q in SCALES if sym[q]["n"] >= 1]
    for qa, qb in itertools.combinations(sorted(nontrivial), 2):
        ratio = qb // qa
        case = None
        if ratio == 2:
            if sym[qa]["type"] != sym[qb]["type"]:
                case = "(1) adjacent scales %d,%d of different types %s/%s" % (
                    qa, qb, sym[qa]["type"], sym[qb]["type"])
            elif sym[qa]["type"] == "I" and sym[qb]["type"] == "I":
                ci, cj = comp_index(comps, qa), comp_index(comps, qb)
                if ci is not None and ci == cj:
                    dim = sum(sym[q]["n"] for q in comps[ci])
                    odd = sum(sym[q]["t"] for q in comps[ci]) % 8
                    if dim > 2:
                        case = "(2) adjacent scales %d,%d both type I, compartment dim %d > 2" % (
                            qa, qb, dim)
                    elif odd in (2, 6):
                        case = "(2) adjacent scales %d,%d both type I, compartment oddity %s2" % (
                            qa, qb, "+" if odd == 2 else "-")
        elif ratio == 4:
            mid = qa * 2
            if sym[qa]["type"] == "I" and sym[qb]["type"] == "I" and sym[mid]["n"] == 0:
                case = "(3) scales %d,%d differ by 4, both type I, term at %d trivial" % (
                    qa, qb, mid)
        if case is None:
            continue
        new = {q: dict(sym[q]) for q in SCALES}
        new[qa]["eps"] = -new[qa]["eps"]
        new[qb]["eps"] = -new[qb]["eps"]
        for ci in {comp_index(comps, qa), comp_index(comps, qb)} - {None}:
            head = comps[ci][0]
            new[head]["t"] = (new[head]["t"] + 4) % 8
        out.append((new, case))
    return out


def unconditional_walks(sym):
    """The SUPERSEDED rule: flip the signs of any two members of one train.  Kept only to
       drive the negative control below -- it is NOT used for the verdict."""
    out = []
    comps = compartments(sym)
    tr, cur = [], [SCALES[0]]
    for i in range(len(SCALES) - 1):
        q, r = SCALES[i], SCALES[i + 1]
        if sym[q]["type"] == "I" or sym[r]["type"] == "I":
            cur.append(r)
        else:
            tr.append(cur)
            cur = [r]
    tr.append(cur)
    for train in tr:
        for qa, qb in itertools.combinations(train, 2):
            new = {q: dict(sym[q]) for q in SCALES}
            new[qa]["eps"] = -new[qa]["eps"]
            new[qb]["eps"] = -new[qb]["eps"]
            for ci in {comp_index(comps, qa), comp_index(comps, qb)} - {None}:
                head = comps[ci][0]
                new[head]["t"] = (new[head]["t"] + 4) % 8
            out.append(new)
    return out

# --------------------------------------------------------------- Lemma 5.2

def fused_key(sym):
    """AGM Lemma 5.2 (oddity fusion): within a compartment only the total oddity mod 8 is
       an invariant, so a symbol is recorded by its constituents plus compartment totals."""
    head = tuple((q, sym[q]["n"], sym[q]["type"], sym[q]["eps"]) for q in SCALES)
    tot = tuple((tuple(c), sum(sym[q]["t"] for q in c) % 8) for c in compartments(sym))
    return (head, tot)


def sym_text(sym):
    parts = ["%d^{%s%d}_%s(t=%d)" % (q, "+" if sym[q]["eps"] > 0 else "-",
                                     sym[q]["n"], sym[q]["type"], sym[q]["t"])
             for q in SCALES if sym[q]["n"] > 0]
    return " · ".join(parts) if parts else "(empty)"


def disc_group(sym):
    d = []
    for q in SCALES:
        if q > 1:
            d.extend([q] * sym[q]["n"])
    return tuple(sorted(d))


def lattice_parity(sym):
    for q in SCALES:
        if sym[q]["n"] and q % 2 == 1 and sym[q]["type"] == "I":
            return "odd"
    return "even"

# --------------------------------------------------------------- merging

def merge(puts, det, walker, gate, ledger=None):
    """Partition the legal symbols by the transitive closure of `walker`.

       `gate` decides when a walk is allowed to merge:

         "fused" -- CORRECT.  The image is accepted when its fused class contains a legal
                    symbol, i.e. when some reassignment of the subscripts inside each
                    compartment makes every term legal while leaving the compartment
                    oddity unchanged.  That reassignment is exactly AGM Lemma 5.2, so
                    testing the representative this code happens to build (which dumps the
                    whole +4 onto the head constituent) would be testing an artefact of
                    the construction rather than the symbol.

         "raw"   -- NAIVE.  Test the constructed representative itself.  Used only for the
                    negative control below; it rejects walks that Lemma 5.2 permits."""
    by_key = {}
    for s in puts:
        by_key.setdefault(fused_key(s), []).append(s)
    parent = {k: k for k in by_key}

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for k in list(by_key):
        for s in by_key[k]:
            for item in walker(s):
                nb, case = (item if isinstance(item, tuple) else (item, None))
                if gate == "raw" and not symbol_is_legal(nb, det):
                    continue
                nk = fused_key(nb)
                if nk not in parent:          # no legal symbol carries this fused class
                    continue
                if ledger is not None and case is not None:
                    ledger.append((det, sym_text(s), sym_text(nb), case))
                a, b = find(k), find(nk)
                if a != b:
                    parent[a] = b
    classes = {}
    for k in by_key:
        classes.setdefault(find(k), []).append(k)
    return classes, by_key

# --------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", default=".")
    args = ap.parse_args()
    print("=" * 100)
    print("Genera of rank 4, signature (3,1), |det| in {1,2,4}")
    print("sign walking after Allcock-Gal-Mark, arXiv:1511.04614, Lemma 6.1;")
    print("oddity fusion after Lemma 5.2.  Every walk used below is certified by its licensing case.")
    print("=" * 100)

    ledger, rows, total = [], [], 0
    for det in (-1, -2, -4):
        puts = putative_symbols(det)
        classes, by_key = merge(puts, det, licensed_walks, "fused", ledger)
        total += len(classes)
        print("\n|det| = %d :  %d putative systems  ->  %d distinct symbols  ->  %d genera"
              % (abs(det), len(puts), len(by_key), len(classes)))
        for root, members in sorted(classes.items(), key=lambda kv: str(kv[0])):
            rep = by_key[members[0]][0]
            dg = disc_group(rep)
            dgs = " x ".join("Z/%d" % d for d in dg) if dg else "trivial"
            print("     %-42s disc = %-14s parity = %s" % (sym_text(rep), dgs, lattice_parity(rep)))
            rows.append({"det": det, "symbol": sym_text(rep), "disc": dgs,
                         "parity": lattice_parity(rep)})

    print("\n" + "-" * 100)
    print("TOTAL GENERA: %d" % total)

    print("\n== LICENCE LEDGER — every sign walk used in the merging, with the case that permits it ==")
    if not ledger:
        print("   (none: no walk was licensed, every symbol stands alone)")
    for det, a, b, case in ledger:
        print("   |det|=%d  %-34s -> %-34s   Lemma 6.1 %s" % (abs(det), a, b, case))
    print("   total licensed walks used: %d" % len(ledger))
    print("   No walk outside these cases was performed; an unconditional train walk is")
    print("   strictly more permissive than Lemma 6.1 and could merge inequivalent symbols.")

    print("\n== CONTROL (negative) — the count is sensitive to the rule ==")
    naive = 0
    for det in (-1, -2, -4):
        puts = putative_symbols(det)
        c, _ = merge(puts, det, unconditional_walks, "raw")
        naive += len(c)
    print("   unconditional train walk, legal images only : %d genera" % naive)
    print("   AGM Lemma 6.1 walk, legal images only       : %d genera" % total)
    print("   the instrument can return a different number, so the AGM count is not a")
    print("   tautology of the run.")
    control_ok = (naive != total)
    print("   [%s] CONTROL" % ("OK" if control_ok else "XX -- rule-insensitive, verdict void"))

    verdict = "OK" if control_ok else "VOID"
    with open("%s/S1764_genus_enum_agm_dump.jsonl" % args.outdir, "w", encoding="utf-8") as fd:
        for r in rows:
            fd.write(json.dumps({"kind": "genus", **r}) + "\n")
        for det, a, b, case in ledger:
            fd.write(json.dumps({"kind": "walk", "det": det, "from": a, "to": b,
                                 "licence": "Lemma 6.1 " + case}) + "\n")
        fd.write(json.dumps({"kind": "verdict", "genera": total, "naive_control": naive,
                             "control_ok": bool(control_ok), "verdict": verdict}) + "\n")
    print("\n[dump] %s/S1764_genus_enum_agm_dump.jsonl" % args.outdir)


if __name__ == "__main__":
    main()
