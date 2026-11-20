# w2a_validation

Validation suite for the white paper

**Two-Coupling One-Loop Matching and the Running Time Unit Under Spatial Blocking in SU(N) Lattice Yang–Mills**
(Weak-Coupling Series, Exercise W2a)

Roger Cooperman — https://orcid.org/0009-0009-5011-2135
A.E. Petrunina V.V. — https://orcid.org/0009-0002-4595-736X
XISOSIS, Vienna, 20 November 2026

- GitHub: https://github.com/teknoquants/Weak-Coupling-Matching-and-the-Running-Time-Unit-Under-Spatial-Blocking-in-SU-N-Lattice-Yang-Mills
- DOI: https://doi.org/10.5281/zenodo.24395118

## Run

```
pip install -r requirements.txt
python run_all.py          # exit code 0 only if every check passes
```

## Checks (46)

| Script | Paper | Content | Checks |
|---|---|---|---|
| `validate_part1_time_unit.py` | Proposition 1, Corollary 2 | H(gE, gB) = ξH(g) as operators (one SU(2) plaquette, character basis); gap scales by ξ; same ground state; scaled kinetic ratio; lifting convention | 5 |
| `validate_part2_matching.py` | Lemma 3, checkpoints | ξ_out = 1 + κx + O(x²) and g_out⁻² = y − A − κ̄ + O(x) (symbolic); κ̄ = 0 from the common decrement; input-clock test; A for SU(2), SU(3) | 6 |
| `validate_part3_drift.py` | Theorem 4, Corollary 5 | Gamma-function formula for the quadratic product; power law K^(κ/A); exponent 0.776663; sign cases; running checkpoint y_n = b + An + O(log n) | 5 |
| `validate_part4_running_clock.py` | Lemma 6, Corollary 7, Theorem 8, Corollary 9 | (1 − ρ)s = κ; obstruction for a repeated multiplier; C_ξ bound on 4000 random admissible steps; exact one-loop ratio; simplified constant 6A\|s\|; telescoping and tail bound | 7 |
| `validate_part5_threshold.py` | Proposition 10, Proposition 11, Theorem 12 | Loss arithmetic; endpoint clock on one plaquette; exp(−c/g^δ) ≤ D_δ g⁴; Σ4⁻ʲ = 4/3; Σg⁴ ≤ g∗⁴ + g∗²/A₋ | 5 |
| `validate_part6_composition_threshold.py` | Theorems 13, 15, Corollary 16, Proposition 18 | Time-dilation covariance of Feshbach and transfer extractions; ξ(g) = 1/c_ref(g) passes (4.1) exactly; κ = 0 and s = −c₁ (symbolic); β_D, g²_match(β_D) = 5.701, step counts 91 and 46 | 7 |
| `validate_numerical_tables.py` | Table 1 | All 48 entries of the P₂/P₄ tables for SU(2), SU(3); gamma formula for P₄; quartic limit 1.139718 | 11 |

## Licence

MIT (see `LICENSE`). Every script carries the MIT header. Copyright (c) 2026 Roger Cooperman; A.E. Petrunina V.V.
