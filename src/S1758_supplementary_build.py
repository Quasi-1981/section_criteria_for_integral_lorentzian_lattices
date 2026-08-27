#!/usr/bin/env python3
# author: B (lane-B chair a7296aa8), S1758.  Naryad #49 step 1 -- build the two supplementary files
# promised in the printed text (P7 §3.5 and §7.4).
#   (a) the 17 putative 2-adic symbol systems and the merge map 17 -> 7   [promised in §3.5]
#   (b) full discriminant-form Grams of the NON-CYCLIC rows of §7.4       [promised in §7.4]
# Numbers come only from already-accepted machinery: S1764 (enumeration) and S1674/S1678 (cost).
# Handles: 0.  No physical words.
# RUN LINE:  python child-3.1/src/S1758_supplementary_build.py --outdir <export/supplementary>
import argparse, os, sys, itertools
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from sympy import Matrix, Rational
from S1764_genus_enum_agm import (putative_symbols, licensed_walks, merge, fused_key,
                                  sym_text as _sym_text, SCALES, disc_group, lattice_parity)
from S1674_bravais_cost_full import disc_form
from S1677_h7_strong import TYPES
from S1678_cost_formula import cost_formula

NONCYCLIC_ROWS = [
    ("hP", "Z/3 x Z/3", 9),
    ("cI", "Z/4 x Z/4", 16),
    ("oC", "Z/2 x Z/8", 16),
    ("hR", "Z/2 x Z/10", 20),
    ("mC", "Z/2 x Z/14", 28),
    ("oI", "Z/2 x Z/2 x Z/24", 96),
]


def sym_text(sym):
    return _sym_text(sym)


def build_a(path):
    lines = []
    lines.append("# Supplementary S1 — the 17 putative 2-adic symbol systems and the merge map 17 → 7")
    lines.append("")
    lines.append("Promised in §3.5 of the paper (P7.Prop1, step (1)).")
    lines.append("")
    lines.append("For rank 4, signature (3,1) and `|det| ∈ {1,2,4}` only the prime 2 carries a")
    lines.append("non-trivial Jordan structure: for odd `p` the lattice is `p`-adically unimodular,")
    lines.append("`q = 1` is a square, so `k_p = 0` and the `p`-excess vanishes. The oddity formula")
    lines.append("therefore reduces to")
    lines.append("")
    lines.append("```")
    lines.append("    Σ_q t_q + 4·k_2  ≡  signature  =  2   (mod 8)")
    lines.append("```")
    lines.append("")
    lines.append("Below, each constituent is written `q^{±n}_{type}(t=oddity)`. Zero-dimensional")
    lines.append("constituents are omitted from the display but are present in the symbol (type II,")
    lines.append("sign +, oddity 0) and do take part in the compartment/train structure.")
    lines.append("")
    lines.append("**`t` here is the oddity of a Jordan constituent** (SPLAG §7.4 notation) and has")
    lines.append("nothing to do with the timelike vector `t` of the criterion (P7.Def4).")
    lines.append("")
    total_put = 0
    total_gen = 0
    for det in (-1, -2, -4):
        puts = putative_symbols(det)
        classes, by_key = merge(puts, det, licensed_walks, "fused")
        groups = {}
        for root, keys in classes.items():
            groups[root] = [s for k in keys for s in by_key[k]]
        total_put += len(puts)
        total_gen += len(groups)
        lines.append("")
        lines.append("## |det| = %d — %d putative systems → %d genera"
                     % (abs(det), len(puts), len(groups)))
        lines.append("")
        lines.append("| # | putative symbol system | merges into genus |")
        lines.append("|--:|:--|--:|")
        gidx = {}
        for i, (ck, members) in enumerate(sorted(groups.items(), key=lambda kv: str(kv[0])), 1):
            gidx[ck] = i
        n = 0
        for ck, members in sorted(groups.items(), key=lambda kv: str(kv[0])):
            for s in members:
                n += 1
                lines.append("| %d | `%s` | **G%d** |" % (n, sym_text(s), gidx[ck]))
        lines.append("")
        lines.append("| genus | canonical representative | disc group | parity |")
        lines.append("|:--|:--|:--|:--|")
        for ck, members in sorted(groups.items(), key=lambda kv: str(kv[0])):
            rep = members[0]
            dg = disc_group(rep)
            dgs = " x ".join("Z/%d" % d for d in dg) if dg else "trivial"
            lines.append("| **G%d** | `%s` | `%s` | %s |"
                         % (gidx[ck], sym_text(rep), dgs, lattice_parity(rep)))
    lines.append("")
    lines.append("## Total")
    lines.append("")
    lines.append("| | putative systems | genera |")
    lines.append("|:--|--:|--:|")
    lines.append("| **rank 4, sig (3,1), `\\|det\\| ∈ {1,2,4}`** | **%d** | **%d** |" % (total_put, total_gen))
    lines.append("")
    lines.append("The merging is done by the two class-preserving operations, taken in the form")
    lines.append("proved by **Allcock–Gal–Mark**, *The Conway–Sloane calculus for 2-adic lattices*,")
    lines.append("`arXiv:1511.04614`:")
    lines.append("")
    lines.append("- **oddity fusion** (their Lemma 5.2) — the subscripts inside a compartment may be")
    lines.append("  reassigned in any way that keeps every term legal and leaves the compartment's")
    lines.append("  oddity unchanged;")
    lines.append("- **sign walking** (their Lemma 6.1) — the signs of two *nontrivial* terms may be")
    lines.append("  negated together, with 4 added to the oddity of each compartment containing one")
    lines.append("  of them, in exactly three cases: adjacent scales of different types; adjacent")
    lines.append("  scales both of type I whose compartment has dimension > 2 or oddity ±2; or")
    lines.append("  scales differing by a factor of 4 with a trivial term between them.")
    lines.append("")
    lines.append("Conway and Sloane's own statement of the canonical form — \"at most one sign per")
    lines.append("train\" — is **not correct**, as Allcock–Gal–Mark show with the counterexample")
    lines.append("`[128¹ 256⁻¹]₄`; their corrected normalisation cuts trains into *signways*. The")
    lines.append("merging above therefore performs only the walks Lemma 6.1 licenses, and the probe's")
    lines.append("log names, for each walk it uses, which of the three cases permits it.")
    lines.append("")
    lines.append("Without the merging the count is wrong: `|det| = 2` would read as four genera")
    lines.append("instead of one.")
    lines.append("")
    lines.append("Reproduce: `python S1764_genus_enum_agm.py` (probe shipped in `src/`).")
    with open(path, "w", encoding="utf-8") as fd:
        fd.write("\n".join(lines) + "\n")
    return total_put, total_gen


def build_b(path):
    lines = []
    lines.append("# Supplementary S2 — discriminant-form Grams of the non-cyclic rows of §7.4")
    lines.append("")
    lines.append("Promised in §7.4 of the paper. For a cyclic `disc(L)` the single number `b̄(x,x)`")
    lines.append("printed in the table determines the form; for the **non-cyclic** rows it does not,")
    lines.append("so the full form is given here.")
    lines.append("")
    lines.append("**How to read it.** For a lattice with Gram `G` in a basis `e₁,e₂,e₃`, the classes")
    lines.append("of the dual basis vectors `e_i^∨` generate `disc(L) = L^∨/L`, and the discriminant")
    lines.append("form on those generators is exactly")
    lines.append("")
    lines.append("```")
    lines.append("    b̄(e_i^∨, e_j^∨) = (G^{-1})_{ij}   mod ℤ")
    lines.append("```")
    lines.append("")
    lines.append("so the matrix `G^{-1} mod ℤ` below **is** the Gram of the discriminant form; the")
    lines.append("relations among the generators are those imposed by `G`, and the elementary divisors")
    lines.append("of `G` give the group. Everything is exact rational arithmetic.")
    lines.append("")
    for key, group, order in NONCYCLIC_ROWS:
        G = TYPES[key]
        M = Matrix(G)
        Ginv = M.inv()
        elts, b, ordf, d = disc_form(G)
        c, wit = cost_formula(G)
        divs = [abs(int(x)) for x in Matrix(G).T.rref()[0].diagonal()] if False else None
        from sympy.matrices.normalforms import smith_normal_form
        S = smith_normal_form(M)
        ed = [abs(int(S[i, i])) for i in range(S.rows) if abs(int(S[i, i])) != 1]
        diag_vals = sorted({Fraction(int(Rational(b(x, x)).p), int(Rational(b(x, x)).q)) % 1
                            for x in elts})
        isotropic = sum(1 for x in elts if Rational(b(x, x)) % 1 == 0)
        lines.append("## %s — `disc(L) = %s`, |disc| = %d, cost = %d, witness (n,m) = (%d,%d)"
                     % (key, group, order, c, wit[0], wit[1]))
        lines.append("")
        lines.append("```")
        lines.append("  lattice Gram G =")
        for row in G:
            lines.append("        [ " + "  ".join("%4d" % v for v in row) + " ]")
        lines.append("")
        lines.append("  discriminant form  b̄ = G^{-1} mod Z  on the dual-basis generators:")
        for i in range(3):
            cells = []
            for j in range(3):
                v = Rational(Ginv[i, j]) % 1
                cells.append("%6s" % (("%d/%d" % (v.p, v.q)) if v != 0 else "0"))
            lines.append("        [ " + "  ".join(cells) + " ]")
        lines.append("")
        lines.append("  elementary divisors of G (the group)   : %s" % (ed,))
        lines.append("  values b̄(x,x) over all of disc(L)      : %s"
                     % ([str(v) for v in diag_vals],))
        lines.append("  # of b-isotropic elements (b̄(x,x)=0)   : %d of %d" % (isotropic, order))
        lines.append("```")
        lines.append("")
    lines.append("Reproduce: the Grams are those of `S1677_h7_strong.TYPES`; the discriminant form")
    lines.append("and the cost are recomputed by `S1674_bravais_cost_full.disc_form` and")
    lines.append("`S1678_cost_formula.cost_formula` (probes shipped in `src/`).")
    with open(path, "w", encoding="utf-8") as fd:
        fd.write("\n".join(lines) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    pa = os.path.join(args.outdir, "S1_genus_enumeration_17_to_7.md")
    pb = os.path.join(args.outdir, "S2_discriminant_form_grams.md")
    tp, tg = build_a(pa)
    build_b(pb)
    print("S1 written: %s   (%d putative -> %d genera)" % (pa, tp, tg))
    print("S2 written: %s   (%d non-cyclic rows)" % (pb, len(NONCYCLIC_ROWS)))
    assert tp == 17 and tg == 7, "GATE FAILED: expected 17 -> 7, got %d -> %d" % (tp, tg)
    print("GATE: 17 -> 7 reproduced OK")


if __name__ == "__main__":
    main()
