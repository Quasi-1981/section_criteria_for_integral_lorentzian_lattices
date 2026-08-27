# Supplementary S1 — the 17 putative 2-adic symbol systems and the merge map 17 → 7

Promised in §3.5 of the paper (P7.Prop1, step (1)).

For rank 4, signature (3,1) and `|det| ∈ {1,2,4}` only the prime 2 carries a
non-trivial Jordan structure: for odd `p` the lattice is `p`-adically unimodular,
`q = 1` is a square, so `k_p = 0` and the `p`-excess vanishes. The oddity formula
therefore reduces to

```
    Σ_q t_q + 4·k_2  ≡  signature  =  2   (mod 8)
```

Below, each constituent is written `q^{±n}_{type}(t=oddity)`. Zero-dimensional
constituents are omitted from the display but are present in the symbol (type II,
sign +, oddity 0) and do take part in the compartment/train structure.

**`t` here is the oddity of a Jordan constituent** (SPLAG §7.4 notation) and has
nothing to do with the timelike vector `t` of the criterion (P7.Def4).


## |det| = 1 — 1 putative systems → 1 genera

| # | putative symbol system | merges into genus |
|--:|:--|--:|
| 1 | `1^{+4}_I(t=2)` | **G1** |

| genus | canonical representative | disc group | parity |
|:--|:--|:--|:--|
| **G1** | `1^{+4}_I(t=2)` | `trivial` | odd |

## |det| = 2 — 4 putative systems → 1 genera

| # | putative symbol system | merges into genus |
|--:|:--|--:|
| 1 | `1^{+3}_I(t=1) · 2^{+1}_I(t=1)` | **G1** |
| 2 | `1^{+3}_I(t=3) · 2^{+1}_I(t=7)` | **G1** |
| 3 | `1^{-3}_I(t=1) · 2^{-1}_I(t=5)` | **G1** |
| 4 | `1^{-3}_I(t=3) · 2^{-1}_I(t=3)` | **G1** |

| genus | canonical representative | disc group | parity |
|:--|:--|:--|:--|
| **G1** | `1^{+3}_I(t=1) · 2^{+1}_I(t=1)` | `Z/2` | odd |

## |det| = 4 — 12 putative systems → 5 genera

| # | putative symbol system | merges into genus |
|--:|:--|--:|
| 1 | `1^{+2}_I(t=0) · 2^{+2}_I(t=2)` | **G1** |
| 2 | `1^{+2}_I(t=2) · 2^{+2}_I(t=0)` | **G1** |
| 3 | `1^{-2}_I(t=2) · 2^{-2}_I(t=4)` | **G1** |
| 4 | `1^{-2}_I(t=4) · 2^{-2}_I(t=2)` | **G1** |
| 5 | `1^{+2}_I(t=2) · 2^{+2}_II(t=0)` | **G2** |
| 6 | `1^{-2}_I(t=6) · 2^{-2}_II(t=0)` | **G2** |
| 7 | `1^{+2}_II(t=0) · 2^{+2}_I(t=2)` | **G3** |
| 8 | `1^{-2}_II(t=0) · 2^{-2}_I(t=6)` | **G3** |
| 9 | `1^{+3}_I(t=1) · 4^{+1}_I(t=1)` | **G4** |
| 10 | `1^{-3}_I(t=5) · 4^{-1}_I(t=5)` | **G4** |
| 11 | `1^{+3}_I(t=3) · 4^{+1}_I(t=7)` | **G5** |
| 12 | `1^{-3}_I(t=7) · 4^{-1}_I(t=3)` | **G5** |

| genus | canonical representative | disc group | parity |
|:--|:--|:--|:--|
| **G1** | `1^{+2}_I(t=0) · 2^{+2}_I(t=2)` | `Z/2 x Z/2` | odd |
| **G2** | `1^{+2}_I(t=2) · 2^{+2}_II(t=0)` | `Z/2 x Z/2` | odd |
| **G3** | `1^{+2}_II(t=0) · 2^{+2}_I(t=2)` | `Z/2 x Z/2` | even |
| **G4** | `1^{+3}_I(t=1) · 4^{+1}_I(t=1)` | `Z/4` | odd |
| **G5** | `1^{+3}_I(t=3) · 4^{+1}_I(t=7)` | `Z/4` | odd |

## Total

| | putative systems | genera |
|:--|--:|--:|
| **rank 4, sig (3,1), `\|det\| ∈ {1,2,4}`** | **17** | **7** |

The merging is done by the two class-preserving operations, taken in the form
proved by **Allcock–Gal–Mark**, *The Conway–Sloane calculus for 2-adic lattices*,
`arXiv:1511.04614`:

- **oddity fusion** (their Lemma 5.2) — the subscripts inside a compartment may be
  reassigned in any way that keeps every term legal and leaves the compartment's
  oddity unchanged;
- **sign walking** (their Lemma 6.1) — the signs of two *nontrivial* terms may be
  negated together, with 4 added to the oddity of each compartment containing one
  of them, in exactly three cases: adjacent scales of different types; adjacent
  scales both of type I whose compartment has dimension > 2 or oddity ±2; or
  scales differing by a factor of 4 with a trivial term between them.

Conway and Sloane's own statement of the canonical form — "at most one sign per
train" — is **not correct**, as Allcock–Gal–Mark show with the counterexample
`[128¹ 256⁻¹]₄`; their corrected normalisation cuts trains into *signways*. The
merging above therefore performs only the walks Lemma 6.1 licenses, and the probe's
log names, for each walk it uses, which of the three cases permits it.

Without the merging the count is wrong: `|det| = 2` would read as four genera
instead of one.

Reproduce: `python S1764_genus_enum_agm.py` (probe shipped in `src/`).
