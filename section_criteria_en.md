# Section Criteria for Integral Lorentzian Lattices and the Cost of Realization

**Vladimir Sobol** · ORCID [0009-0006-9829-7931](https://orcid.org/0009-0006-9829-7931) · Independent Researcher

---

## Abstract

For an integral lattice `Λ_W` of signature (3,1) and a positive definite integral `L` of rank 3 we
give a **criterion** for `L` to be a **primitive timelike section** of `Λ_W` by a primitive vector of
norm `−n`. The criterion runs on three data — **signature, bilinear discriminant form and parity**
(type I/II) — and states: such a section exists **if and only if** there is an anti-isometry of glue
subgroups `γ : H_L → H_t` whose graph `Γ` satisfies `Γ^⊥/Γ ≅ disc(Λ_W)` **as an isometry of forms**,
and whose glued overlattice has **the same parity** as `Λ_W`. Both requirements are load-bearing: the
order-only version of the first is insufficient, and without the second the statement is **false** —
a counterexample on the ambient `diag(1,1,2,−2)` is exhibited explicitly in §3.4. The criterion is
stated **relative to a fixed ambient** and is expressed **exclusively** in discriminant data and
parity — **the Gram matrix does not enter it**. The reverse implication rests on two checkable
conditions: **(T0)** ("the named data determine the genus of the ambient") and **(T0c)** ("the genus
of the ambient contains one class"); for rank 4, signature (3,1) and `|det| ≤ 4` — a scope covering
all four ambients of this work — both are **proved**: (T0) by a complete enumeration of genera (there
are exactly seven, and the triple of data separates them pairwise, §3.5), (T0c) by Eichler's classical
theorem with a margin of `5⁶` (§3, Remark 1).

Two direct checks. **First**, substituting `disc(Λ_W) = 0` returns the condition `⟨1/n⟩` of the
unimodular case, that is, the theorem of preprint-6: the generalisation **collapses back onto what is
already in print**. **Second**, the parity rider shows that this boundary is **sharp**: for
`p ≡ 1 (mod 8)` the condition `⟨1/n⟩` is no longer sufficient, because two lattices with an
**identical** discriminant form fall into different ambients, and what separates them is one extra
bit — parity. The case `p = 3` is thereby isolated (`II_{3,1}` does not exist), so the printed result
of preprint-6 **stands as it is**.

From the criterion it follows that the **minimal discriminant of an ambient** — the "cost" of a
lattice `L` — depends **only on the bilinear discriminant form** of `L`, and is therefore an
**invariant of its genus**; we give for it a **closed expression**
`cost(L) = min{ |disc(L)|·n_min(m)/m² : m | exp(disc L) }`.
Finiteness of the range here is **structural**, not computational, so the expression contains no
search bound — and that is exactly why it is **robust against truncation**, which in earlier bounded
computations three times produced inflated values. Evaluating the expression on fourteen
representatives gives a table of costs; it is indexed by the **discriminant form**, not by Bravais
type and not by the order of the discriminant, because cost is **not** an invariant of the type and is
**not monotone** in the order. In particular, the centrings of the cubic lattice acquire exact
addresses: **FCC — 3**, **BCC — 8**, while the blocking of lower values for FCC splits into **two
different** mechanisms (insolubility of `k² ≡ 3 (mod 4)` at `d = 1`; insolubility of `2k² ≡ 3 (mod 4)`
at `d = 2`), whereas `d = 3` is blocked by nothing at all.

There is no physical reading anywhere in this work: the signature (3,1) appears as an arithmetic
condition — in particular as the **condition of applicability** of the classical theorem on the genus
of indefinite forms, on which the criterion leans in one explicitly named step.

---

## §1 · Statement of the problem

Preprint-6 (Zenodo `10.5281/zenodo.22068307`) classified the primitive timelike sections of the
**unimodular** lattice `I_{p,1}` by the condition `⟨1/n⟩` (`P6.Thm1(7)`).

**The question of this work:** what replaces that condition for an **arbitrary** integral Lorentzian
lattice.

**What exactly is being asked, and under which quantifier:** the lattice `Λ_W` is taken as **GIVEN**.
The whole work answers the question

> "is `L` a primitive timelike section of **THIS** lattice `Λ_W`",

and **not** "does there exist **some** ambient into which `L` embeds". The difference is not
rhetorical: under the second reading the answer is almost always "yes" and the question is empty; it
acquires content **only** for a fixed ambient — and it is precisely this quantifier that the statement
of **P7.Thm1** carries. The one place where the ambient is **not** fixed but **sought** is the
definition of **cost** (P7.Def5): there the minimum is taken over all admissible `Λ_W`, and this is
said in the definition itself.

**Scope, stated once and at the start:** this work contains **no physical reading whatsoever**. No
object is read as "time", "energy" or "mass"; the signature (3,1) appears **exclusively** as an
arithmetic condition — and in one place (Remark 1) as the **condition of applicability** of a
classical theorem.

---

## §2 · Definitions

**P7.Def1 — discriminant form.** For a non-degenerate integral `L`: `disc(L) := L^∨/L` with the
**bilinear** form `b_L : disc(L) × disc(L) → ℚ/ℤ`. The abbreviation `b̄` is used where the
discriminant form is meant as an **invariant** (with no lattice named in the subscript).

> **Why `b mod ℤ` and not `q mod 2ℤ`:** the quadratic version is defined only for **even** lattices,
> while `⟨−n⟩` is **odd for odd `n`** — and that is not an edge case but the main one. The whole
> criterion is therefore stated in the bilinear form. The price of this choice is named explicitly:
> the bilinear form **does not see parity**, so parity enters the data of the criterion as a
> **separate bit** (§3).

**P7.Def2 — glue subgroups and anti-isometry.** Subgroups `H_L ⊆ disc(L)`, `H_t ⊆ disc(⟨−n⟩)` are
called **glue subgroups**. An isomorphism `γ : H_L → H_t` is an **anti-isometry** if
`b_t(γx, γy) = −b_L(x, y)`.

> The term is not invented: "glue" is the established word of the construction (📖 Conway–Sloane,
> SPLAG ch. 4: *glue group*, *glue vectors*; Nikulin §1.5). We name only the **subgroups** on which
> the gluing is defined, because all the corollaries refer to them.

**P7.Def3 — the graph.** `Γ := {(x, γx) : x ∈ H_L} ⊆ disc(L) ⊕ disc(⟨−n⟩)`. `Γ` is **isotropic** if
`b_L ⊕ b_t` vanishes on `Γ` — for the graph of an anti-isometry this holds **identically**:
`b_L(x,y) + b_t(γx,γy) = b_L(x,y) − b_L(x,y) = 0`.

> **Two sides of one equivalence:** here, "anti-isometry ⟹ isotropy"; in (Thm1.3), the converse.
> Together: for subgroups with injective projections, **isotropic graphs are exactly the graphs of
> anti-isometries**.

**P7.Def4 — section.** An integral positive definite `L` of rank 3 (abstract) is a **primitive
timelike section** of `Λ_W` by `t` if `t ∈ Λ_W` is **primitive** (`t/k ∉ Λ_W` for every integer
`k ≥ 2`), `⟨t,t⟩ = −n < 0`, and there is an isometric embedding `ι : L ↪ Λ_W` with
`ι(L) = t^⊥ ∩ Λ_W`. The shorthand `L = t^⊥ ∩ Λ_W` is henceforth read through this `ι` — the reverse
implication proves the **existence of an embedding**, not an equality of sets.

> **Primitivity is required of `t` — and this requirement is load-bearing; for `L` it is automatic.**
> Automatic: `Λ_W ∩ (L⊗ℚ) = Λ_W ∩ t^⊥ = L` by definition — the intersection of a lattice with a
> subspace is always primitive. Load-bearing: without the requirement on `t` the theorem of §3 is
> **false**. Witness: `Λ_W = I_{3,1}`, `t = 2e₄`, `n = 4`, `L = t^⊥ ∩ Λ_W = ℤ³` — the section exists,
> but the right-hand leg of the criterion fails: `disc(L) = 0` ⟹ the only anti-isometry is the empty
> one, `Γ = 0`, and `Γ^⊥/Γ ≅ ℤ/4 ≇ disc(I_{3,1}) = 0`. The same `ℤ³` is a section by the primitive
> `e₄` with `n = 1` — the norm of a section is well defined only on a primitive vector.

**P7.Def5 — cost.**
```
   cost(L) := min { |det Λ_W| : Λ_W integral, of signature (3,1),
                    and L a primitive timelike section of Λ_W by some vector }
```

**P7.Def6 — `j_min` and the admissible glue orders.**
```
   D(L)      := { m ≥ 1 : m | exp(disc L) }
   j_min(m)  := min { j ≥ 1 : ∃ x ∈ disc(L) of order m, ∃ a unit k mod m,
                              j·k² ≡ b_L(x,x)·m  (mod m) } ,   m ∈ D(L)
   n_min(m)  := m · j_min(m)
```
> **Why `exp` and not `|disc|`:** taking `m | |disc(L)|` would admit an `m` with no element of order
> `m`, i.e. the term under `min` would be **undefined**. Example: `ℤ/2 × ℤ/2` has `4 | 4` but
> `exp = 2`. For the computation this changes nothing — **for the statement it changes everything**.
> It is a separate definition because `j_min` stands **in the statement** of the theorem and is
> referred to by the corollaries: a quantity named only inside a theorem gets cited by transcription.

---

## §3 · P7.Thm1 — the gluing criterion

> **P7.Thm1.** Let `L` be a positive definite integral lattice of rank 3, `n ≥ 1`, and let `Λ_W` be an
> integral lattice of signature (3,1) **(both fixed)** satisfying the conditions
>
> **(T0)** the triple of data **(signature, bilinear discriminant form, parity)** determines the genus
> of `Λ_W` uniquely among integral lattices — a condition **checkable by a finite comparison** (for
> rank 4, sig (3,1), `|det| ≤ 4` — in particular for all four ambients of this work — it is **proved**
> by an enumeration of genera, P7.Prop1, §3.5);
>
> **(T0c)** the genus of `Λ_W` contains **one isometry class**. (Class = spinor genus for indefinite
> lattices of rank ≥ 3 — 📖 Eichler; a genus with more than one spinor genus in rank 4 requires
> `|det| ≥ 5⁶ = 15625` — 📖 SPLAG Ch. 15, Cor. 22; hence for `|det| ≤ 4` condition (T0c) holds with
> room to spare — Remark 1, leg 1. Without (T0c) the reverse leg yields only "`Λ` lies in the genus of
> `Λ_W`", not an isometry.)
>
> Then **`L` is a primitive timelike section of `Λ_W` by a primitive vector of norm `−n`** (P7.Def4)
> ⟺ **there exists an anti-isometry of glue subgroups `γ : H_L → H_t`** whose graph `Γ` satisfies
> **both** conditions:
>
> **(T1)** `Γ^⊥/Γ ≅ disc(Λ_W)` — **as an isometry of forms**, not merely as an isomorphism of groups;
> **(T2)** the glued overlattice `Λ(Γ)` has **the same parity** (type I/II) as `Λ_W`; parity is read
> off the norms of the glue vectors `mod 2ℤ`.
>
> (Isotropy of `Γ` is automatic by P7.Def3.)

**(T1) "as a form", not "as a group".** The order-only version fails at the first non-cyclic
discriminant (**P7.Cor3**, §6): `ℤ/4` and `ℤ/2×ℤ/2` have the same order 4 and are not isometric. But
the **group** alone is not enough either: on `ℤ/6` there live **two non-isometric forms**, and they
give **different costs** — `b̄` on the generator equal to `5/6` gives 2, and `1/6` gives 1 (§7.2
against the row *oP* in §7.4). A group is not a form.

**(T2) parity is the third datum, and it is not decorative.** The criterion runs on the triple
**signature + bilinear discriminant form + parity**. With these data the equivalence holds for **all
four** ambients of this work: for three of them an even competitor **does not exist** (excluded by
Milgram's formula — table **§3.4**), and for the fourth the **bit separates** the odd `diag(1,1,2,−2)`
from the even `A₁²⊕U`. Without **(T2)** the statement is **false** precisely on the fourth — the
explicit pair is in §3.4.

### Remark 1 — what the theorem stands on (to be read BEFORE the proof)

The step requires the data of the criterion to determine `Λ_W` **uniquely**, not merely up to genus.
It stands on **two** legs, and they must be seen separately — their conflation was the reason the
second went unchecked for a long time.

**Leg 1 — "class = genus" for our objects.**
The classical theorem (**Eichler**; **Kneser**) for **indefinite** integral lattices of **rank ≥ 3**
states: **class = SPINOR genus**. 📖 a named import, **not re-derived** in this work. Our ambient has
rank **4**, so the rank condition is met. The classical theorem gives "class = **spinor genus**",
whereas the step needs "class = **genus**". These are **not identical**: a genus may contain
**several** spinor genera. For our objects the difference is void, and this is a count, not an
estimate: in rank 4 a genus with **more than one** spinor genus requires `|det| ≥ 5⁶ = 15625`
(📖 SPLAG Ch. 15, Cor. 22); our four ambients have `|det| ∈ {1, 2, 4, 4}` — a margin of almost four
thousandfold, so the condition holds **trivially**, not "at the edge". This is exactly condition
**(T0c)** of the theorem: within the scope of this work it is proved; for an arbitrary `Λ_W` it stands
in the statement as a separate hypothesis (a genus with several classes would make the conclusion of
(Thm1.11) false).

**Leg 2 — parity, and the limit of what it gives.**
"Signature + **bilinear** discriminant form" **do not determine the genus**: parity is a genus
invariant that the bilinear form **does not see**. Counterexample: `A₁³⊕⟨−2⟩` glues by different `Γ`
into the **even** `A₁²⊕U` and the **odd** `diag(1,1,2,−2)` — same signature, same discriminant group
and `b̄`, different genera. That is exactly why parity stands in the data of the criterion as the
separate condition **(T2)**.

> **The limit of leg 2 — and the scope in which it is closed.** From the fact that parity is
> **necessary** it does not follow that the triple "signature + `b̄` + parity" is **sufficient** to
> determine the genus of an **arbitrary** lattice: the 2-adic genus symbol of an odd lattice is in
> general finer than a single global bit (it carries the type and the oddity of each Jordan
> constituent). Sufficiency therefore stands as the separate condition **(T0)** — and for rank 4,
> sig (3,1), `|det| ≤ 4` it is **proved** by a complete enumeration of genera (**P7.Prop1**, §3.5):
> there are exactly seven genera, seven triples, and the correspondence is one-to-one. For arbitrary
> `det`, (T0) is **not claimed** — this work gives no general characterisation of (T0)-lattices and
> does not need one.

**Without these two legs** the criterion yields only: `Λ` lies in **the same genus** as `Λ_W` — that
is, "a section of **some** lattice of the genus" instead of "a section of `Λ_W` **itself**". The whole
force of the "⟺" rests on this step.

**Counterfactual — why signature (3,1) is a condition of applicability.** For **definite** lattices,
class and genus **may differ**, and the argument of (Thm1.11) is **unavailable** there; hence for a
definite ambient the criterion in the stated form would be unjustified. Signature (3,1) is a
**condition of applicability**, not a convenience of notation.

### P7.Lem1 — the gluing identity (bilinear, with no parity restriction)

> **P7.Lem1.** Let `M` be a non-degenerate integral lattice, `Λ ⊇ M` an integral overlattice of finite
> index, `Γ := Λ/M ⊆ M^∨/M = disc(M)`, and `Γ^⊥ := {x ∈ disc(M) : b̄(x, Γ) = 0 in ℚ/ℤ}`. Then
> `disc(Λ) ≅ Γ^⊥/Γ` — as groups with a bilinear form `mod ℤ`.

*Proof.* `M ⊆ Λ` (given); integrality of `Λ` gives `Λ ⊆ Λ^∨`; from `M ⊆ Λ` we get `Λ^∨ ⊆ M^∨` — hence
the chain `M ⊆ Λ ⊆ Λ^∨ ⊆ M^∨`. For `x ∈ M^∨`, the condition `b(x, Λ) ⊆ ℤ` is equivalent to
`b(x, γ) ∈ ℤ` for all `γ ∈ Γ` (on `M` integrality is automatic), i.e. `(x mod M) ∈ Γ^⊥`; hence
`Λ^∨/M = Γ^⊥`. Since `Λ/M = Γ`, the third isomorphism theorem gives
`disc(Λ) = Λ^∨/Λ = (Λ^∨/M)/(Λ/M) = Γ^⊥/Γ`, and the bilinear form on `Λ^∨/Λ` is induced from that of
`M` identically. ∎

> **Why a lemma here rather than a citation.** The classical exposition of this identity (Nikulin)
> runs through the **quadratic** discriminant form `q mod 2ℤ`, defined only for **even** lattices —
> whereas our case is odd. The group-theoretic and bilinear side of the identity, as the proof shows,
> **needs no parity at all** — so for odd lattices it is proved rather than cited. A numerical check
> of both sides (by independent computations) was carried out on three odd gluings and one even
> control.

### P7.Lem2 — gluing by an isotropic subgroup

> **P7.Lem2.** Let `M` be a non-degenerate integral lattice, `Γ ⊆ disc(M)` a subgroup, and
> `Λ_Γ := {x ∈ M^∨ : x mod M ∈ Γ}`. Then:
> **(a)** `Λ_Γ` is a subgroup of `M^∨` with `M ⊆ Λ_Γ` and `Λ_Γ/M = Γ`; **`Λ_Γ` is integral ⟺ `Γ` is
> isotropic**; **(b)** for `M = L ⊕ ⟨t⟩`: `L` is primitive in `Λ_Γ` ⟺ the projection
> `Γ → disc(⟨t⟩)` is **injective**; symmetrically, `⟨t⟩` is primitive in `Λ_Γ` ⟺ the projection
> `Γ → disc(L)` is **injective**.

*Proof.* **(a)** `Λ_Γ` is the preimage of the subgroup `Γ` under `M^∨ → disc(M)`, hence a subgroup;
`M` is the preimage of zero, so `M ⊆ Λ_Γ` and `Λ_Γ/M = Γ`. Integrality: for `x, y ∈ Λ_Γ` the value
`b(x, y) mod ℤ` equals `b̄(x̄, ȳ)`, hence `b(Λ_Γ, Λ_Γ) ⊆ ℤ` ⟺ `b̄` vanishes on `Γ` (pairs with `x̄ = 0`
are integral automatically: `x ∈ M`, `y ∈ M^∨`). **(b)** `Λ_Γ ∩ (L⊗ℚ) ⊆ M^∨ ∩ (L⊗ℚ) = L^∨`, so an
element of the intersection has image of the shape `(x̄, 0) ∈ Γ`; hence `Λ_Γ ∩ (L⊗ℚ) = L` ⟺ `Γ`
contains no non-zero element with vanishing `⟨t⟩`-component ⟺ `Γ → disc(⟨t⟩)` is injective.
Symmetrically for `⟨t⟩`. ∎

> The lemma closes in one place the four steps that lean on it: (Thm1.2)/(Thm1.3) in the forward leg
> (there `Λ_W = Λ_Γ` for `Γ = Λ_W/M`), (Thm1.7)/(Thm1.10) in the reverse leg, and the chain (Thm2.9).

### §3.1 · Proof: the forward leg (section ⟹ γ)

**(Thm1.1)** `t` primitive, `⟨t,t⟩ = −n`, `L = t^⊥ ∩ Λ_W` (primitive automatically, P7.Def4) ⟹
`L ⊕ ⟨t⟩ ⊆ Λ_W` — a sublattice of finite index.
**(Thm1.2)** `Γ := Λ_W/(L ⊕ ⟨t⟩)`; integrality of `Λ_W` ⟹ `Γ` is **isotropic** *(P7.Lem2(a), the
direction "integral ⟹ isotropic"; for the graph of an anti-isometry this also holds identically by
P7.Def3)*.
**(Thm1.3)** Primitivity of `L` (automatic) gives injectivity of `Γ → disc(⟨t⟩)`; primitivity of
`⟨t⟩` (the requirement of P7.Def4 — exactly here is where it bears) gives injectivity of
`Γ → disc(L)` *(both — P7.Lem2(b))* ⟹ `Γ` is the graph of an isomorphism `γ`. The anti-isometry is
read **bilinearly**, from isotropy on a pair of distinct elements `(x,γx)`, `(y,γy) ∈ Γ`:
```
   0 = (b_L ⊕ b_t)((x,γx),(y,γy)) = b_L(x,y) + b_t(γx,γy)   ⟹   b_t(γx,γy) = −b_L(x,y)
```
for **all** `x, y` — exactly what P7.Def2 demands. *(It cannot be derived through the diagonal
`b(x,x)`: over `ℚ/ℤ` polarisation does not recover a bilinear form on 2-torsion.)*
**(Thm1.4)** `disc(Λ_W) ≅ Γ^⊥/Γ` as a form — condition **(T1)**. *(P7.Lem1.)*
**(Thm1.4′)** The norms of the glue vectors of `Λ_W/(L⊕⟨t⟩)` `mod 2ℤ` give the parity of `Λ_W` —
condition **(T2)** holds identically, since the overlattice is `Λ_W` itself.
**(Thm1.5)** `disc(⟨t⟩) ≅ ℤ/n` with `b(k,k) = −k²/n`. ∎

### §3.2 · Proof: the reverse leg (γ ⟹ Λ_W)

**(Thm1.6)** `M := L ⊕ ⟨−n⟩`, signature `(3,0)+(0,1) = (3,1)`.
**(Thm1.7)** An isotropic `Γ ⊆ disc(M)` determines an **integral** overlattice `Λ := Λ_Γ ⊇ M`,
`Λ/M ≅ Γ` *(P7.Lem2(a))*.
**(Thm1.8)** The signature is preserved ⟹ `Λ` has (3,1).
**(Thm1.9)** `disc(Λ) ≅ Γ^⊥/Γ` *(P7.Lem1 — the same lemma as in (Thm1.4))*, and
`Γ^⊥/Γ ≅ disc(Λ_W)` — by condition **(T1)**.
**(Thm1.10)** Primitivity — **two different properties, each with its own ground** *(P7.Lem2(b))*:
`L` primitive ⟸ injectivity of `Γ → disc(⟨t⟩)`; `⟨t⟩` primitive ⟸ injectivity of `Γ → disc(L)`; both
injectivities are automatic, since `Γ` is the graph of an **isomorphism** `γ`.
**(Thm1.11)** `Λ` and `Λ_W` have the same signature, discriminant form **and parity** (the last by
condition **(T2)**). By condition **(T0)** these data determine the genus of `Λ_W` uniquely ⟹ `Λ` and
`Λ_W` lie in one genus; by condition **(T0c)** that genus contains one class ⟹ `Λ ≅ Λ_W`. For the
ambients of this work both conditions are proved: (T0) by the enumeration of genera (P7.Prop1, §3.5),
(T0c) by 📖 **Eichler** plus the `5⁶` margin (Remark 1, leg 1). ∎

### §3.3 · The measurement that accompanies the proof

Two **independent routes** to one number: **route (1)** — a direct enumeration of primitive `t` with
`t² = −n ≤ 4`; **route (2)** — the gluing formula `|disc(L)|·n = |disc(Λ_W)|·index²`.

| ambient | `\|disc\|` | what the enumeration produced | (1) against (2) |
|:--|:--|:--|:--|
| `I_{3,1}` **(control, disc = 0)** | 1 | `n → \|disc(L)\| = n`, cyclic `ℤ/n` | **100 %** |
| `diag(1,1,1,−2)` | 2 | `1→2`, `2→{1, ℤ/2×ℤ/2}`, `3→6`, `4→ℤ/2×ℤ/4` | **100 %** |
| `diag(1,1,1,−4)` | 4 | `1→4`, `2→ℤ/2×ℤ/4`, `3→12`, `4→1` | **100 %** |
| `A₃⊕⟨−1⟩` — **the same ambient** as the row above | 4 | `1→4` (`= A₃`), `2→ℤ/2×ℤ/4`, `3→12` | **100 %** |

The last two rows describe **one** lattice: `diag(1,1,1,−4) ≅ A₃ ⊕ ⟨−1⟩` — both have a vector of norm
`−1`, in both the complement is `A₃`, the theta series are identical. The identity of their columns is
not a coincidence but a consequence of this isometry.
⟹ **the measurement covered three distinct ambients out of four.** The fourth, `diag(1,1,2,−2)`,
enters the work as the witness of §6 and §3.4; the test "(1) against (2)" could not have caught
anything on it in any case — both routes enumerate **actual** sections, so a false reverse leg is
invisible to them. That is exactly why a counterexample (§3.4) was needed, and not a wider search.

The unimodular row is a **control, not an illustration**: the enumeration reproduced `⟨1/n⟩` **by
itself**, without being tuned to it. This is **§4 in action**.

### §3.4 · When (T2) holds by itself: exclusion of the even competitor

Condition **(T2)** is non-trivial only where an **even competitor exists** — that is, where there is
an even lattice with the same signature and the same `b̄` as `Λ_W`. This question is settled by
Milgram's **necessary** condition, and the answer for the four ambients of this work fits into a table.

**The mechanism.** An even lattice carries a quadratic refinement `q : disc → ℚ/2ℤ` with `b_q = b̄`.
Milgram's formula (📖 §10) relates the signature of the lattice to the Gauss sum of the refinement:
```
        GS(q) = Σ_x e^{πi·q(x)} = √|disc| · e^{πi·σ(q)/4} ,      σ(Λ) ≡ σ(q)  (mod 8)
```
Signature (3,1) gives `σ = 2`. Hence: if **no** refinement `q` of the form `b̄` has `σ(q) ≡ 2 (mod 8)`,
then an even lattice with these data **does not exist** — there is no competitor, and (T2) holds
automatically. There are finitely many refinements (on a cyclic generator of order `m` the value
`q(g)` is determined by `b̄(g,g)` up to `ℤ/2ℤ`), so the enumeration is **complete**, not selective.

| ambient | `disc` | `σ(q)` — all attainable | even competitor |
|:--|:--|:--|:--|
| `I_{3,1}` | trivial | `{0}` | **excluded** |
| `diag(1,1,1,−2)` | `ℤ/2` | `{1, 7}` | **excluded** |
| `diag(1,1,1,−4)` | `ℤ/4` | `{3, 7}` | **excluded** |
| `diag(1,1,2,−2)` | `ℤ/2 × ℤ/2` | `{0, **2**, 6}` | **NOT excluded** |

**Three conclusions.**

**(i)** On the first three, `2 ∉ σ(q)` ⟹ there is no even competitor ⟹ **(T2) holds by itself**, and
there the criterion works on the signature and `b̄` alone.
**(ii)** The first row is **Remark 2** (§5) in wider notation: for a trivial discriminant the only
refinement `q = 0` gives `σ(q) = 0`, i.e. "an even unimodular lattice of signature (3,1) does not
exist". Two statements that looked different are **one mechanism**, `σ mod 8`.
**(iii)** On the fourth, `σ(q) = 2` **is attainable** — and a competitor is indeed there: `A₁²⊕U`,
constructed explicitly below. So the table **does not over-squeeze**: it does not exclude what exists.
It is on this row that (T2) carries its weight.

**The explicit pair of the fourth row.** `L = A₁³` (Gram `diag(2,2,2)`), `n = 2`,
`Λ_W = diag(1,1,2,−2)`: condition **(T1)** holds — `H_L = ⟨g₃⟩`, `γ : g₃ ↦ 1` is an anti-isometry
(`½ = −½ mod ℤ`), and `Γ^⊥/Γ ≅ (ℤ/2)²` with `b = diag(½,½)`, which is **isometric** to
`disc(diag(1,1,2,−2))`. But the gluing yields `Λ = A₁²⊕U` — the glue vector `w = (e₃+t)/2` has
`w² = 0`, and the pair `w, (e₃−t)/2` is hyperbolic ⟹ `Λ` is **even**, whereas `Λ_W` is **odd**.
Condition **(T2)** is precisely what filters this case out.

> **On the status of this table.** Milgram gives a **necessary** condition: `2 ∉ σ(q)` proves the
> **non-existence** of a competitor strictly. The converse — that for `σ(q) = 2` a competitor does
> exist — does **not** follow from Milgram, and for the fourth row it is proved **not** by a Gauss sum
> but by the **explicit construction** of `A₁²⊕U`. Two different methods, each used where it proves
> something.

### §3.5 · P7.Prop1 — condition (T0) proved by an enumeration of genera

> **P7.Prop1.** Integral lattices of rank 4, signature (3,1), `|det| ∈ {1, 2, 4}` exist, up to genus,
> in **exactly seven** kinds, and the triple **(discriminant group, `b̄`, parity)** takes on them
> **seven distinct values**. Hence the triple determines the genus: condition **(T0)** holds for
> **every** lattice of this scope — in particular for all four ambients of this work.

*Proof.* The genus of an indefinite lattice is determined by the signature and the p-adic symbols; for
`|det| ∈ {1,2,4}` only the 2-adic one is non-trivial. The enumeration is run from above — from the
existence conditions — and is closed from below by explicit lattices.

**(1) The enumeration is complete.** The existence conditions for a system of p-adic symbols
(📖 SPLAG Ch. 15 §7.7: the determinant condition, the oddity formula, the conditions on each Jordan
constituent; Theorem 11 — such a system is realised by an integral form) give, for signature (3,1) and
`|det| ∈ {1,2,4}`, **17 putative systems**. The merging of equivalent systems — *oddity fusion* and
*sign walking* — is carried out **per Allcock–Gal–Mark** (📖: Lemma 5.2 — fusion only when the image
is a legal symbol; Lemma 6.1 — the three admissible cases of a sign walk; **Theorem 6.2 — two symbols
are isometric ⟺ related by a chain of these operations**), so the count is exact **in both
directions**: the operations neither merge inequivalent systems nor leave equivalent ones unmerged:

| `\|det\|` | putative systems | canonical genera |
|:--|--:|--:|
| 1 | 1 | **1** |
| 2 | 4 | **1** |
| 4 | 12 | **5** |
| **total** | **17** | **7** |

Without canonicalisation the count is wrong: `|det| = 2` would read as "four genera" instead of one.

> **Why the merging follows AGM rather than the "canonical form" of SPLAG §7.6.** The uniqueness
> claim of §7.6 ("at most one minus sign per train") is **false** in general — the error is named by
> Allcock–Gal–Mark (§6, following Allcock, Memoirs AMS 220); the correction is *signways*. Step (1)
> therefore rests not on §7.6 but on the AGM operations with Theorem 6.2, and its correctness does
> not hang on the absence of a known bad configuration: **all 17 systems are merged exclusively by
> operations licensed by Lemma 5.2 and Lemma 6.1** — no forbidden sign walk occurs in any chain —
> and Theorem 6.2 identifies the classes of this equivalence with the isometry classes, of which
> there are exactly `1 + 1 + 5 = 7`. A measured control: an unconditional train walk (without the
> conditions of Lemma 6.1) yields **eight** "genera" instead of seven — the count is sensitive to the
> rule, so the agreement with step (3) is not a tautology of the instrument.

**(2) The seven genera by name.**

| `\|det\|` | canonical 2-adic symbol | discriminant group | parity | representative |
|:--|:--|:--|:--|:--|
| 1 | `1^{+4}_I (t=2)` | trivial | odd | `I_{3,1}` |
| 2 | `1^{+3}_I (t=1) · 2^{+1}_I (t=1)` | `ℤ/2` | odd | `diag(1,1,1,−2)` |
| 4 | `1^{+3}_I (t=1) · 4^{+1}_I (t=1)` | `ℤ/4` | odd | `diag(1,1,4,−1)`, `b̄ = ¼` |
| 4 | `1^{+3}_I (t=3) · 4^{+1}_I (t=7)` | `ℤ/4` | odd | `diag(1,1,1,−4) ≅ A₃⊕⟨−1⟩`, `b̄ = ¾` |
| 4 | `1^{+2}_I (t=0) · 2^{+2}_I (t=2)` | `ℤ/2×ℤ/2` | **odd** | `diag(1,1,2,−2)` |
| 4 | `1^{+2}_I (t=2) · 2^{+2}_{II} (t=0)` | `ℤ/2×ℤ/2` | odd | `I₂⊕U(2)`, `b̄` with zero diagonal |
| 4 | `1^{+2}_{II} (t=0) · 2^{+2}_I (t=2)` | `ℤ/2×ℤ/2` | **even** | `A₁²⊕U` |

> **A homonym is named:** `t` **inside the symbols** is the oddity of a Jordan constituent (SPLAG §7.4
> notation); it has nothing in common with the timelike vector `t` of the criterion (P7.Def4).

The representatives of rows 3 and 6 are explicit: `diag(1,1,4,−1)` has discriminant `ℤ/4` with
`b̄(g,g) = (¼)²·4 = ¼`; `I₂⊕U(2)` (where `U(2)` is the Gram `[[0,2],[2,0]]`) has discriminant `(ℤ/2)²`
with zero diagonal and `b̄(x,y) = ½`; both are odd, of signature (3,1) and `|det| = 4` — and Prop1
itself pins each into its row: the triple coincides with exactly one genus. The full list of the 17
putative systems and the merge map `17 → 7` — supplementary **S1**
([`supplementary/S1_genus_enumeration_17_to_7.md`](supplementary/S1_genus_enumeration_17_to_7.md)).

**(3) Seven triples — and this is counted without symbols.** The columns "discriminant group / `b̄` /
parity" above are read off the symbols and could in principle separate more finely than a genuine
isometry of `b̄` — in which case "seven distinct triples" would be an artefact of the notation. The
separation is therefore taken from an **independent** count from below: on explicit Grams of this
scope the triple — with `b̄` as an actual finite form, isometry decided by exhaustive search — takes
**seven distinct values**, and each is realised. The triple is a genus invariant; there are exactly
seven genera and seven values ⟹ the map "genus → triple" is an **injection**. ∎

**Two checks on the table.** Rows 5 and 7 are the pair of §3.4 (`diag(1,1,2,−2)` against `A₁²⊕U`):
the same discriminant group, different parity — and **two different** canonical symbols;
canonicalisation does not merge what §3.4 separated by an explicit construction. Row 1 is the only
genus at `|det| = 1`, and it is **odd**: the even variant drops out because type II gives oddity
`t ≡ 0`, while the oddity formula demands `t ≡ 2` — this is Remark 2 (§5), re-derived from the
existence conditions rather than quoted.

**Which invariant separates each pair.** The even competitors of the first three ambients are excluded
by `σ(q) mod 8` (§3.4); at `|det| = 4` the groups `ℤ/4` and `ℤ/2×ℤ/2` are separated by the group
structure; the two forms on `ℤ/4` — by the diagonal of `b̄` (`[0,¼,¼]` against `[0,¾,¾]`); the two odd
forms on `ℤ/2×ℤ/2` — by the count of `b`-isotropic elements (1 against 3 non-zero); the pair of rows
5/7 — by the parity bit **(T2)**, and this is the only place where it carries weight.

## §4 · P7.Cor1 — return to preprint-6

> **P7.Cor1.** `disc(Λ_W) = 0` ⟺ `Γ` is Lagrangian ⟺ `disc(L) ≅ ℤ/n` with `b̄_L = ⟨1/n⟩` — that is,
> **P6.Thm1(7)**. The condition `⟨1/n⟩` here and everywhere is read **bilinearly** (P7.Def1), in P6's notation:
> `b_L(g,g) ≡ a/n (mod ℤ)` on the generator, `a` a quadratic residue mod `n`.

**(Cor1.1)** `disc(Λ_W) ≅ Γ^⊥/Γ` ⟹ `disc(Λ_W) = 0 ⟺ Γ^⊥ = Γ` (Lagrangian).
**(Cor1.2)** Order count:
```
   |Γ| = |H_L| = |H_t|                                   (i)
   |Γ|² = |disc(L)| · n                                  (ii)
   |H_L| ≤ |disc(L)| ,  |H_t| ≤ n                        (iii)
   ⟹  H_L = disc(L),  H_t = disc(⟨−n⟩),  |disc(L)| = n   (iv)
```
equality of the products under upper bounds **forces equality of each factor**.
> From (iv) it follows not merely that the orders are equal, but that `disc(L)` is **isometric** to the
> cyclic form — this is precisely where the requirement "as a form" starts to work (necessity — **§6**).

**(Cor1.3)** From (Cor1.2) `H_L = disc(L)`, so `γ` is an anti-isometry of **all** of `disc(L)` onto `ℤ/n`:
`b_L(x, y) = −b_t(γx, γy)`. On the generator `g` we have `γ(g) = k·g_t` with `k` a unit mod `n` (`γ` being
an isomorphism), and by (Thm1.5) `b_t(g_t, g_t) = −1/n`, hence
```
   b_L(g, g) = −k²·b_t(g_t, g_t) = −k²·(−1/n) = k²/n   (mod ℤ)
```
— that is, `b̄_L = ⟨1/n⟩` in the bilinear notation (`a = k²` a quadratic residue). The quadratic
`q mod 2ℤ` does not arise anywhere here — the transfer is purely bilinear. ∎

**Main argument of the work:** preprint-6 is **neither bypassed nor corrected** — it **follows** from
P7.Thm1 by substituting `disc = 0`. This is a **generalisation that folds back into the printed result**.

---

## §5 · P7.Cor2 — parity rider: the boundary of preprint-6 is exact

> **P7.Cor2.** If the glued overlattice `Λ(Γ)` is **even**, then `L` and `⟨−n⟩` are **both** even.
> **Corollary:** for `p ≡ 1 (mod 8)` the condition `⟨1/n⟩` is **not sufficient** — it does not distinguish
> an even ambient from an odd one.

*Proof of necessity.* `M = L ⊕ ⟨−n⟩ ⊆ Λ(Γ)`, so every vector of `M` has even
norm. ∎

> **The converse does NOT hold:** the parity of the overlattice is read off from the norms of the
> **glue vectors** `mod 2ℤ`, not from the parity of the summands. Witness — our own pair from §3.4: the **even**
> `M = A₁³⊕⟨−2⟩` glues along one class into the **even** `A₁²⊕U`, and along another into the **odd**
> `diag(1,1,2,−2)`. This is exactly why parity entered the criterion's data as a separate condition **(T2)**,
> rather than as a consequence of the summands' parity.

| lattice | **discriminant form** | **parity** | ambient |
|:--|:--|:--|:--|
| `E₈ ⊕ ⟨2⟩` | `ℤ/2`, `b = ½` | **even** | only `II_{9,1}` |
| `I₈ ⊕ ⟨2⟩` | `ℤ/2`, `b = ½` — **the same** | **odd** | `I_{9,1}` |

⟹ the discriminant form **does not determine** the ambient; an extra **bit** is needed. Hence the boundary
of `P6.Thm1(7)` **IS EXACT**, not conservative.

### Remark 2 — the case `p = 3`

> **`II_{3,1}` does not exist**, so for `p = 3` the even branch is empty, the rider cannot fire, and the
> condition `⟨1/n⟩` **is sufficient**.

*Proof.* An even unimodular indefinite lattice exists ⟺
`sig ≡ 0 (mod 8)`; for (3,1) `sig = 2 ≢ 0`. ∎
📖 **the classification of indefinite even unimodular [lattices]** — the named leg of this remark.


**Positively and exactly:** the boundary of the condition `⟨1/n⟩` **does not touch** the case `p = 3`, so the
printed result of preprint-6 **stands as is**; this work merely outlined its scope from outside.

---

## §6 · P7.Cor3 — the form is necessary

> **P7.Cor3.** The order-only version of the criterion (`|disc(L)|·n = |disc(Λ_W)|·index²`) **is
> insufficient**.

The ambient is given explicitly, it is `4×4` diagonal **and is itself the witness**:
```
        Λ_W = diag(1, 1, 2, −2) ,   |disc| = 4 ,   signature (3,1)
        disc(Λ_W) = ℤ/2 × ℤ/2 ,   form b = {0, 0, ½, ½}
```

| `t²` | `\|disc(L)\|` | `index²` | gluing arithmetic |
|:--|:--|:--|:--|
| −1 | 4 | 1 | ✓ |
| −2 | 2 | 1 | ✓ |
| −3 | 12 | 9 | ✓ |
| −4 | 16 | 16 | ✓ |

**The gluing arithmetic passes on ALL four — and precisely because of that it decides nothing.** The form
distinguishes: `ℤ/4` and `ℤ/2×ℤ/2` have the same order 4, but **are not isometric**; `disc(Λ_W)` coincides
precisely with `ℤ/2×ℤ/2`.

**The cost of the omission:** without the words "as a form" **P7.Thm1 is false** — it becomes a statement
about order and fails at the first non-cyclic discriminant. The corollary that saves the theorem from
falsity must say this about itself.

---

## §7 · P7.Thm2 — cost: the law and the closed formula

### §7.1 · The law

> ### **P7.Thm2(a) — the cost depends ONLY on the bilinear discriminant form of `L`.**
> If `b̄_{L₁} ≅ b̄_{L₂}` (with the same rank and signature), then `cost(L₁) = cost(L₂)`.
> **Corollary (citable form):** the cost is a **genus invariant** of `L` — the genus determines the
> discriminant form. Invariance by genus is **strictly weaker**: the witness of §7.2(i) has the same cost
> for **different** genera.

> **Order of proof: (a) is derived from (b), not from Thm1.** This is not a matter of style. After the
> repair the data of Thm1 contains **parity (T2)**, so the phrase "the criterion is expressed purely in
> discriminant data" can **no longer** be taken as a premise — it has ceased to be exact. Instead, **the
> closed formula of §7.3 demonstrably consumes only `b̄`**: neither parity nor the Gram matrix enters it.
> So (a) is proved **after** (b) and **from** (b) — constructively, not rhetorically.

**(Thm2.1)** The formula of P7.Thm2(b) is expressed via `|disc(L)|`, `exp(disc L)`, and the values of `b_L`
on the elements of `disc(L)` — that is, **exclusively** via `b̄` as a form; neither the Gram matrix nor
parity enters it. Isometric `b̄` give the same value of the formula, and the formula equals the cost ⟹
`cost` is a function of `b̄`. ∎
**(Thm2.2)** The genus determines the signature and the discriminant form ⟹ on the genus the cost is
constant — a **genus invariant**. ∎

> **Why there is no import here.** Previously this step relied on "signature + discriminant form **is** the
> genus" (📖 Nikulin). This reference is both **unnecessary** and **too narrow**: only the trivial direction
> is needed (genus ⟹ discriminant form), while the converse for **odd** lattices with a **bilinear** form is
> **false** — our own witness from §5 (`E₈⊕⟨2⟩` and `I₈⊕⟨2⟩`: the same bilinear discriminant form, different
> genera). Nikulin's theorem is formulated for **even** lattices and a **quadratic** form; we do not stretch
> it beyond its scope (§10).

> **Why parity (T2) does not affect the cost — and this is visible from the definition, not from hope.**
> `P7.Def5` minimises over **all** ambients, not over ambients of a prescribed parity: every admissible
> gluing gives **some** integral overlattice (3,1), and it is itself a witness. So the added condition (T2),
> which distinguishes ambients from one another, **does not shift** the minimum. The witness of §7.2(i)
> shows this directly: the pair differs precisely in parity and both have cost `1`.
> ⟹ **the cost half of the work (§7–§8) does not move under the repair of Thm1**, and it has zero imports:
> neither Eichler, nor genus.

### §7.2 · Two different statements — two different witnesses

**(i) The cost is determined by the discriminant form — and it does NOT see either class or genus.**
Witness: **the same discriminant form, the same cost, different lattices — even different genera.**
Here Gram matrices are mandatory (`3×3`, one pair) — the dispute is precisely **about the representative**:
```
   diag(1,1,6)   [[1,0,0],[0,1,0],[0,0,6]]    |Aut| 16   ODD       disc ℤ/6, b̄(g,g) = 1/6   cost 1
   A₂ ⊕ ⟨2⟩      [[2,−1,0],[−1,2,0],[0,0,2]]  |Aut| 24   EVEN      disc ℤ/6, b̄(g,g) = 1/6   cost 1
```
`|Aut|` 16 ≠ 24 ⟹ **not isometric**. Moreover: `diag(1,1,6)` is **odd** (`e₁² = 1`), while
`A₂⊕⟨2⟩` is **even** (the entire Gram diagonal is even), and parity is a genus invariant ⟹ these two
lattices lie in **different genera**. The same discriminant form ⟹ **the same cost** despite this.

> **The witness proves more than "a genus invariant".** The weaker version ("one genus, two classes") would
> have been a special case of this one. The variable on which the cost depends is the **discriminant
> form**; the genus would be a superfluous intermediary that, moreover, does not even coincide here.

**(ii) The cost is NOT an invariant of the Bravais TYPE** — the witness is **opposite in form**: the same
**type**, **different** costs — `A₂⊕⟨2⟩` (hexagonal) → 1 versus the row `hP` (hexagonal) → 3.

> **Why (i) does not prove (ii):** equal cost for two lattices says that the cost does not see **class** or
> **genus**, and says nothing about **type**. Only a divergence **within** a type can refute invariance of
> type.

### §7.3 · The closed formula

> **P7.Thm2(b).** `cost(L) = min { |disc(L)| · n_min(m)/m² : m ∈ D(L) }`, with `D(L)`, `n_min` as per
> P7.Def6.

**Step 1 — FINITENESS** *(as a separate and first step: it is precisely this that makes the expression a
theorem, not a search)*

**(Thm2.3)** `Γ` projects into `disc(L)` injectively (Thm1.3) ⟹ **`Γ` is isomorphic to a subgroup of
`disc(L)`**.
**(Thm2.4)** ⟹ `m := |Γ|` divides `|disc(L)|` (Lagrange), and **more precisely** `m | exp(disc L)`, because
`H_L` is cyclic — this is (Thm2.7) below, and it is proved **independently** of this step, so there is no
circularity. The divisors are **always** finite in number.
**(Thm2.5)** The set in the definition of `j_min(m)` is **non-empty**, and `j_min(m) ≤ m`. Indeed: `m ∈
D(L)` guarantees the existence of `x` of order `m`; at `k = 1` the condition reads as `j ≡ A (mod m)`,
where `A := b_L(x,x)·m mod m` is an integer, because `ord(x) = m`. The class `j mod m` has a
representative in `{1, …, m}` (at `A = 0` this is `j = m`) ⟹ the set is non-empty and the minimum is
`≤ m`. Finite.
**(Thm2.6)** There are **infinitely many** admissible pairs `(m, n)`: for a fixed `m` the congruence of
P7.Def6 admits **the entire residue class** `j mod m`, and with it arbitrarily large `n = m·j`. The
finiteness lives not in the pairs, but in the minimisation: the target quantity `|disc(L)|·n/m² =
|disc(L)|·j/m` **strictly increases in `j`** for fixed `m`, so the minimum over the class is attained at
`j_min(m)` (non-empty and `≤ m` — Thm2.5), and `m` ranges over the finite `D(L)` (Thm2.4) ⟹ the
minimisation reduces to the finite set `{(m, m·j_min(m)) : m ∈ D(L)}`; the minimum exists, is attained,
**without any search bound `D_max`**. ∎

> **The formula is truncation-stable by construction** — not because the search bound `D_max` is taken
> sufficiently large, but because there **is no bound in it at all** (`D(L)` from P7.Def6 is the set of
> divisors of `exp(disc L)`, a structural object, not a truncation parameter).

**Step 2 — the expression counts exactly what is needed**

**(Thm2.7)** **Cyclicity is a consequence, not a convenient special case:** `H_t ⊆ ℤ/n` is **always**
cyclic; `γ` an isomorphism ⟹ `H_L` is **forced** to be cyclic. Hence a single `m` and a single generator
suffice **by structure**.
> **Removes the apparent contradiction with §6:** `disc(L)` **can** be non-cyclic (that is precisely what §6
> is about); what is cyclic is the **glue subgroup**, not the disc.

**(Thm2.7′) — the scalar congruence of Def6 is equivalent to the FULL anti-isometry.** Separately, because `Def6` operates on **a single number** `b_L(x,x)`, whereas `Def2` requires equality of
the forms on **all** pairs. For cyclic `H_L = ⟨x⟩`, `H_t = ⟨y⟩` and `γ(x) = k·y`, bilinearity gives
`b(ax, bx) = ab·b(x,x)`, so
```
   b_t(γ(ax), γ(bx)) = ab·k²·b_t(y,y)   and   −b_L(ax, bx) = −ab·b_L(x,x)
```
for all `a, b` ⟺ `k²·b_t(y,y) = −b_L(x,x)` — **one** equality on the generators. The value of the form on
the generator determines it on the entire cyclic subgroup; hence the scalar condition `Def6` is no weaker
than `Def2`. (Cyclicity is not an assumption, but (Thm2.7).)

Next — the only place with genuine arithmetic: `y = (n/m)·g` gives `b_t(y,y) = −n/m²`, and with `n = m·j`
and multiplication by `m`:
```
        j · k² ≡ b_L(x,x) · m   (mod m)
```
hence `j` ranges over **quadratic shifts by units `k`**, and `j_min` is a question about **quadratic
residues mod `m`**.

**(Thm2.8)** The gluing arithmetic gives `|det Λ_W| = |disc L|·n/m²`; minimisation over `(m,n)` is
precisely the expression.

**(Thm2.9) — a candidate from `Def6` gives precisely an ambient ADMISSIBLE under `Def5`.** The step is
spelled out in full, because `Def5` minimises **not** over integral overlattices in general, but over
those where `L` is a **primitive timelike section**:

```
   admissible gluing (m, j, k)
   ⟹ Γ = graph(γ), γ : H_L → H_t isomorphism             (Thm2.7′)
   ⟹ both projections Γ → disc(L), Γ → disc(⟨−n⟩)
      are injective — automatically, because it is the graph of an ISOMORPHISM
   ⟹ ⟨t⟩ primitive, L primitive                          (Thm1.10)
   ⟹ L = t^⊥ ∩ Λ : primitivity of L gives L = (L⊗ℚ) ∩ Λ,
      and L⊗ℚ = t^⊥⊗ℚ                                    (P7.Def4)
   ⟹ Λ is an admissible ambient under P7.Def5
```
Plus `Γ` is isotropic under P7.Def3 ⟹ `Λ` is **integral**, of signature (3,1) (Thm1.7, Thm1.8). ⟹ **the
formula = the criterion, no wider**: every member of it is realised, and no realisable ambient falls out
from under it. ∎

### §7.4 · Table of values

> **Indexed by the INPUT of the formula — the discriminant FORM, not the Bravais type and not the order of
> the disc.** The table of values must be ordered by what the value depends on. Type is a **label of the
> chosen representative**.
> **The group of the disc is NOT an input.** The `b̄` column is in the table precisely for this reason:
> without it the row "`ℤ/6` → 2" would contradict the pair of §7.2 ("`ℤ/6` → 1"), even though both numbers
> are correct — these are **two different forms on the same group**.

| `disc(L)` — group | `\|disc\|` | `b̄(x,x)` on the generator of the gluing | **cost** | witness `(n, m)` | *representative label* |
|:--|--:|:--|--:|:--|:--|
| `0` | 1 | — | **1** | (1, 1) | *cP* |
| `ℤ/2` | 2 | `1/2` | **1** | (2, 2) | *tP* |
| `ℤ/4` | 4 | `3/4` | **3** | (12, 4) | *cF* "FCC" |
| `ℤ/6` | 6 | `1/3` | **2** | (3, 3) | *oP* |
| `ℤ/8` | 8 | `1/2` | **4** | (8, 4) | *tI* |
| `ℤ/3 × ℤ/3` | 9 | `1/3` | **3** | (3, 3) | *hP* |
| `ℤ/4 × ℤ/4` | 16 | `1/2` | **8** | (8, 4) | *cI* "BCC" |
| `ℤ/2 × ℤ/8` | 16 | `3/8` | **6** | (24, 8) | *oC* |
| `ℤ/2 × ℤ/10` | 20 | `2/5` | **8** | (10, 5) | *hR* |
| `ℤ/24` | 24 | `1/8` | **3** | (8, 8) | *oF* |
| `ℤ/2 × ℤ/14` | 28 | `3/14` | **6** | (42, 14) | *mC* |
| `ℤ/51` | 51 | `1/17` | **3** | (17, 17) | *aP* |
| `ℤ/69` | 69 | `2/69` | **2** | (138, 69) | *mP* |
| `ℤ/2 × ℤ/2 × ℤ/24` | 96 | `5/24` | **20** | (120, 24) | *oI* |

> **How to read the `b̄` column.** This is the value on the generator of the **glue subgroup** (an element
> of order `m`) — it is exactly this that enters the congruence `j·k² ≡ b̄·m (mod m)`. The value on the
> generator of the **whole** disc can be different: for the row *oP* `b̄(g,g) = 5/6` on the generator of
> `ℤ/6`, whereas the table shows `1/3` — the value on `x = 2g` of order 3. **Every displayed witness can
> be checked from its row:** `j = n/m`, `cost = |disc|·n/m²`; for **non-cyclic** forms the independent
> check of **minimality** requires the full discriminant form — supplementary **S2**
> ([`supplementary/S2_discriminant_form_grams.md`](supplementary/S2_discriminant_form_grams.md)).
> **Two forms on `ℤ/6` that justify this column:** `b̄(g,g) = 5/6` (row *oP*) gives cost 2; `b̄(g,g) = 1/6`
> (the pair of §7.2) gives cost 1. One group, different costs — because the forms differ.

*(The table — values of the formula, by a single method; representatives validated by `|Aut|`, 14/14; Gram
matrices — supplementary.)* For **non-cyclic** `disc(L)` the `b̄` column gives the value only on the
generator of the gluing, not the whole form: the row verifies the **upper bound** of the cost, while
minimality is with respect to the full discriminant form; the full Gram matrices of the discriminant forms
of all six non-cyclic rows — supplementary **S2**
([`supplementary/S2_discriminant_form_grams.md`](supplementary/S2_discriminant_form_grams.md)).

**Two things become visible that a table by types was hiding:**
**(i)** at the same order the form distinguishes: `ℤ/4×ℤ/4` → 8 versus `ℤ/2×ℤ/8` → 6;
**(ii)** **the cost is NOT monotone in `|disc|`** — `ℤ/69` → 2, `ℤ/4` → 3, `ℤ/8` → 4: a larger disc can be
**cheaper**. ⟹ **no "size leads" whatsoever**; the cost is a function of the discriminant form, and only
that.

---

---

## §8 · P7.Cor4 — centring addresses

> **P7.Cor4.** **FCC — `cost = 3`**, **BCC — `cost = 8`**; both unattainable on the unimodular
> stratum, but with **different** addresses, and the difference is explained by the **discriminant form**.

### Table 1 — the blocking mechanism: three lines, not one

For FCC `disc = ℤ/4`, `b̄ = 3/4`; the admissible `m ∈ {1, 2, 4}`, and `d = |disc|·j/m` with
`j·k² ≡ b̄·m (mod m)`.

| `d = \|det Λ_W\|` | status | what exactly blocks it |
|:--|:--|:--|
| **d = 1** | blocked | requires `k² ≡ 3 (mod 4)`, but squares mod 4 are `{0,1}`: **quadratic residue** |
| **d = 2** | blocked | both routes fall on **congruences modulo 2**: `m = 2` gives `1 ≢ 0 (mod 2)`, `m = 4` gives `2k² ≡ 3 (mod 4)` — an even number against an odd one |
| **d = 3** | **blocked by nothing** | it **is** the minimum itself — the overlattice is explicit, `\|det\| = 3`, (3,1), `(n,m) = (12,4)` |

> **On the word "parity" in the row `d = 2`.** Here it is the **parity of integers** in the congruence
> at issue, **not** the parity of the lattice (type I/II) from condition **(T2)**. The homonym is named:
> after §3 "parity" in this work is a load-bearing word, and in Table 1 it means **something else**.
> "**Form** blocks, not order" — a conclusion from three lines, not a substitute for them.

> **The status of the number `8` for BCC — so that two levels are not read as one.** `cost(BCC) = 8` is
> obtained **abstractly**: the formula P7.Thm2(b) computes the cost from the `b̄` of `L` itself, and the
> realisability of each term is automatic (Thm2.9). It **does not depend** on the open item of §11 — that
> concerns the **full isometry** `Γ^⊥/Γ` for BCC, i.e. the **identification of the specific ambient**, not
> the value of the cost. The cost answers "how expensive", the identification "which one exactly"; the
> latter remains open.

This is the **native reason** for excluding FCC from the unimodular stratum: earlier it stood as an
**observed obstruction**, now it is **derived** from the criterion.
The **asymmetry `3 < 8`** is explained by **different forms** — and **not** by monotonicity in size (§7.4).

---

## §9 · Measured alongside: `cost` and `n_min` — not one scale

Is the "ambient cost" not another expression for `n_min` (P6)? **No:**

```
        Spearman( n_min , cost )  ≈  0.78     on one set of representatives
        Spearman( n_min , cost )  ≈  0.36     on an alternative set
```

The coefficient depends on the chosen set of representatives and is not cited as a value.

**Conclusion — as two separate statements, each with its own witness** (one is not proved without the
other: two order reversals **do not** refute a functional dependence — a function
can perfectly well have `f(43) = 3`, `f(13) = 8`):

**(i) The order is not consistent (non-monotonicity):**
```
        triclinic       n_min = 43 ,  cost = 3
        rhombohedral    n_min = 13 ,  cost = 8
```
a larger `n_min` — a smaller cost; it is exactly this (and only this) that `Spearman < 1` reflects.

**(ii) `cost` is not a function of `n_min` — a collision:** `n_min` is an invariant of **holohedry** (P6),
while the cost is not: both lattices of the pair in §7.2(ii) are hexagonal, i.e. they have **the same**
holohedry and `n_min = 6`, but costs **1** (`A₂⊕⟨2⟩`) and **3** (the *hP* row). The same argument,
different values — there is no function.


---

## §10 · Relation to prior work and sources

| mechanism | source |
|:--|:--|
| the terminology of gluing (glue group / glue vectors); the discriminant-form theory of the **even** case (the identity `disc(Λ) ≅ Γ^⊥/Γ` for the odd case is proved in this work — **P7.Lem1**) | **V. V. Nikulin**, *Integral symmetric bilinear forms and some of their applications*, Izv. Akad. Nauk SSSR **43** (1979); English transl. Math. USSR Izv. **14** (1980) — §1.5, §1.5.1. **J. H. Conway, N. J. A. Sloane**, *Sphere Packings, Lattices and Groups*, 3rd ed., Springer 1999 — **Ch. 4** |
| class = spinor genus for **indefinite** forms of **rank ≥ 3** | **M. Eichler**, *Quadratische Formen und orthogonale Gruppen*, Springer 1952. **M. Kneser**, *Klassenzahlen indefiniter quadratischer Formen…*, Arch. Math. **7** (1956). Modern exposition: **J. W. S. Cassels**, *Rational Quadratic Forms*, Academic Press 1978 — **Ch. 11**; **O. T. O'Meara**, *Introduction to Quadratic Forms*, Springer — **§102** |
| one spinor genus in the genus for small `\|det\|`, rank 4 | **Conway–Sloane**, SPLAG 3rd ed. — **Ch. 15, Cor. 22** |
| correctness of the Conway–Sloane calculus for 2-adic lattices (the first published proof; **§6 — the correction of the false CS canonical form**, following Allcock, Memoirs AMS 220; it carries the merging of §3.5: Lemma 5.2, Lemma 6.1, Theorem 6.2) | **D. Allcock, I. Gal, A. Mark**, *The Conway–Sloane calculus for 2-adic lattices*, `arXiv:1511.04614` |
| **indefinite** even unimodular (`sig ≡ 0 mod 8`); classification of the odd ones | **J.-P. Serre**, *A Course in Arithmetic*, GTM 7 — **Ch. V**. **J. Milnor, D. Husemoller**, *Symmetric Bilinear Forms*, Springer 1973 — **Ch. II §5** |
| the Milgram formula (`σ` of the discriminant form `mod 8`); exclusion of the even competitor | **Milnor–Husemoller**, *Symmetric Bilinear Forms* — **Appendix**; **SPLAG** — formula **(30)** |
| existence conditions for p-adic symbols — the load-bearing step of **P7.Prop1** (§3.5); the "canonical form" of §7.6 is **not used** in the steps: its uniqueness is false in general (see the AGM row) | **Conway–Sloane**, SPLAG 3rd ed. — **Ch. 15, §7.7** (Theorems 10, 11) |
| genus symbols (for reference; Nikulin's existence theorems §1.9/§1.11 — the even case, **not used** in the steps of this work: all ambients and overlattices are built explicitly) | **Nikulin** §1.9, §1.11 |
| the crystallographic restriction in rank 3 *(historical note)* | 19th-century classic; modern expositions — **W. Scherrer** (1946), **H. S. M. Coxeter**, *Introduction to Geometry* (1969) |

> **THE SCOPE LINE, WITHOUT WHICH THE REFERENCE TO NIKULIN READS WIDER THAN IT IS.**
> Nikulin's discriminant-form theory is formulated for **EVEN** lattices. **Our case is odd** — and
> it works **natively**, because the whole criterion is built from the outset on the **bilinear**
> discriminant form `b mod ℤ` (`P7.Def1`), not on the quadratic `q mod 2ℤ`, which is defined only for
> even lattices.
> ⟹ we do **not** apply Nikulin's result outside its scope: from there we take the **terminology of
> gluing and the identity `disc(Λ) ≅ Γ^⊥/Γ`**, while the odd case carries its own construction. The
> price of this choice — that the bilinear form does not see parity — is paid explicitly: the bit sits
> in the data of the criterion **(T2)**, and is not hidden in the import.

**Theorem numbers are given only where checked against the edition; otherwise — the chapter or section.**

**What is genuinely this work's own:**

1. **the criterion relative to a fixed ambient** with `Γ^⊥/Γ ≅ disc(Λ_W)` **as a form** and with
   **parity** in the data (the odd case);
2. **folding back**: `disc = 0` reproduces the `⟨1/n⟩` condition of preprint-6;
3. **the cost as an invariant of the discriminant form** plus a **closed formula** with a finite
   range — **without a search parameter**;
4. **the addresses of the centrings** with the blocking mechanism separated out (Table 1).

---

## §11 · Limits and the open

### (a) Open

| question | status |
|:--|:--|
| **characterisation of (T0)-lattices** — for which `Λ_W` the triple "signature + `b̄` + parity" determines the genus in **general** form | open for arbitrary `det` and rank ≠ 4; for rank 4, sig (3,1), `\|det\| ≤ 4` it is **proved** by enumeration of genera (P7.Prop1, §3.5) |
| ~~the lemma "no even competitor"~~ | **closed** — written out as the table of §3.4 (`σ mod 8`, all refinements of `q`, the enumeration is complete) |
| **the full form of `Γ^⊥/Γ` for BCC** | the order and factors are computed; the full isometry has not been done. **This does not concern the cost** — `cost(BCC) = 8` is abstract (§8) |
| **γ-freedom for a non-cyclic disc** | open **without a source** |
| **a simpler form of `j_min(m)`** | possibly via quadratic residues / CRT |
| **the cost as an invariant of TYPE** | **closed negatively** — witness §7.2(ii): the same type, different costs |

> **A bridging line** (because `A₂⊕⟨2⟩` works twice): **there, the same discriminant form — the same
> cost even with different genera; here, the same TYPE — a different cost. This is exactly what shows
> that the discriminant form is the variable.**

### (b) Parking, said out loud

The sample `n ≤ 4` does not even reproduce the **known** unimodular image (it gives `{tetragonal, cubic}`
instead of seven; `triclinic` requires `n = 43`) ⟹ a "by genus" difference on it would be a **truncation
artefact**. Fincke–Pohst up to `n = 43` with a budget is needed.

### (c) On method: why the theorem has no search bound

**A bounded search (limit `D_max = 40`) gave inflated cost values — three times.** The reason is each
time the same: for large discriminants the optimal `m` is large, `n = d·m²/|disc|` jumps past the
limit, the row **drops out of the enumeration**. *(Notation: `D_max` is the truncation limit of the old
search; it has nothing in common with the set of divisors `D(L)` in P7.Def6.)*

**The finiteness step (Thm2.3–Thm2.6) exists precisely so that the theorem has no search bound at
all.** `m | exp(disc L)` is a **structural** condition, not a computational one.
**Conclusion of the method: a named structure is not verified by scanning** — if a problem is finite
by structure, this must be **used**, not approximated by a bound.

### (d) What is not our result

**The window (quasiperiodicity) is a separate input, not a consequence of the criterion.**
Finite-order integer isometries of `I_{3,1}` have orders from `{1,2,3,4,6}`: axes of orders **5, 8, 12**
do not preserve (3,1) — **measured**; order **10** is excluded as a consequence (its square has order 5)
— **derived**. A detailed treatment is not part of this work.

> This is a statement about **scope**, not a promise: it says **what is absent**, without "later we
> will show".
