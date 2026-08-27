# Supplementary S2 — discriminant-form Grams of the non-cyclic rows of §7.4

Promised in §7.4 of the paper. For a cyclic `disc(L)` the single number `b̄(x,x)`
printed in the table determines the form; for the **non-cyclic** rows it does not,
so the full form is given here.

**How to read it.** For a lattice with Gram `G` in a basis `e₁,e₂,e₃`, the classes
of the dual basis vectors `e_i^∨` generate `disc(L) = L^∨/L`, and the discriminant
form on those generators is exactly

```
    b̄(e_i^∨, e_j^∨) = (G^{-1})_{ij}   mod ℤ
```

so the matrix `G^{-1} mod ℤ` below **is** the Gram of the discriminant form; the
relations among the generators are those imposed by `G`, and the elementary divisors
of `G` give the group. Everything is exact rational arithmetic.

## hP — `disc(L) = Z/3 x Z/3`, |disc| = 9, cost = 3, witness (n,m) = (3,3)

```
  lattice Gram G =
        [    2    -1     0 ]
        [   -1     2     0 ]
        [    0     0     3 ]

  discriminant form  b̄ = G^{-1} mod Z  on the dual-basis generators:
        [    2/3     1/3       0 ]
        [    1/3     2/3       0 ]
        [      0       0     1/3 ]

  elementary divisors of G (the group)   : [3, 3]
  values b̄(x,x) over all of disc(L)      : ['0', '1/3', '2/3']
  # of b-isotropic elements (b̄(x,x)=0)   : 5 of 9
```

## cI — `disc(L) = Z/4 x Z/4`, |disc| = 16, cost = 8, witness (n,m) = (8,4)

```
  lattice Gram G =
        [    3    -1    -1 ]
        [   -1     3    -1 ]
        [   -1    -1     3 ]

  discriminant form  b̄ = G^{-1} mod Z  on the dual-basis generators:
        [    1/2     1/4     1/4 ]
        [    1/4     1/2     1/4 ]
        [    1/4     1/4     1/2 ]

  elementary divisors of G (the group)   : [4, 4]
  values b̄(x,x) over all of disc(L)      : ['0', '1/2']
  # of b-isotropic elements (b̄(x,x)=0)   : 4 of 16
```

## oC — `disc(L) = Z/2 x Z/8`, |disc| = 16, cost = 6, witness (n,m) = (24,8)

```
  lattice Gram G =
        [    3     1     0 ]
        [    1     3     0 ]
        [    0     0     2 ]

  discriminant form  b̄ = G^{-1} mod Z  on the dual-basis generators:
        [    3/8     7/8       0 ]
        [    7/8     3/8       0 ]
        [      0       0     1/2 ]

  elementary divisors of G (the group)   : [2, 8]
  values b̄(x,x) over all of disc(L)      : ['0', '3/8', '1/2', '7/8']
  # of b-isotropic elements (b̄(x,x)=0)   : 4 of 16
```

## hR — `disc(L) = Z/2 x Z/10`, |disc| = 20, cost = 8, witness (n,m) = (10,5)

```
  lattice Gram G =
        [    3     1     1 ]
        [    1     3     1 ]
        [    1     1     3 ]

  discriminant form  b̄ = G^{-1} mod Z  on the dual-basis generators:
        [    2/5    9/10    9/10 ]
        [   9/10     2/5    9/10 ]
        [   9/10    9/10     2/5 ]

  elementary divisors of G (the group)   : [2, 10]
  values b̄(x,x) over all of disc(L)      : ['0', '2/5', '3/5']
  # of b-isotropic elements (b̄(x,x)=0)   : 4 of 20
```

## mC — `disc(L) = Z/2 x Z/14`, |disc| = 28, cost = 6, witness (n,m) = (42,14)

```
  lattice Gram G =
        [    3     1     0 ]
        [    1     5     0 ]
        [    0     0     2 ]

  discriminant form  b̄ = G^{-1} mod Z  on the dual-basis generators:
        [   5/14   13/14       0 ]
        [  13/14    3/14       0 ]
        [      0       0     1/2 ]

  elementary divisors of G (the group)   : [2, 14]
  values b̄(x,x) over all of disc(L)      : ['0', '3/14', '5/14', '3/7', '1/2', '5/7', '6/7', '13/14']
  # of b-isotropic elements (b̄(x,x)=0)   : 2 of 28
```

## oI — `disc(L) = Z/2 x Z/2 x Z/24`, |disc| = 96, cost = 20, witness (n,m) = (120,24)

```
  lattice Gram G =
        [    4     0     2 ]
        [    0     8     4 ]
        [    2     4     6 ]

  discriminant form  b̄ = G^{-1} mod Z  on the dual-basis generators:
        [    1/3    1/12     5/6 ]
        [   1/12    5/24     5/6 ]
        [    5/6     5/6     1/3 ]

  elementary divisors of G (the group)   : [2, 2, 24]
  values b̄(x,x) over all of disc(L)      : ['0', '5/24', '1/3', '3/8', '1/2', '17/24', '5/6', '7/8']
  # of b-isotropic elements (b̄(x,x)=0)   : 8 of 96
```

Reproduce: the Grams are those of `S1677_h7_strong.TYPES`; the discriminant form
and the cost are recomputed by `S1674_bravais_cost_full.disc_form` and
`S1678_cost_formula.cost_formula` (probes shipped in `src/`).
