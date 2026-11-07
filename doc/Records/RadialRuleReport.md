# Is MHL the right radial rule, and where do its points go? (lane B, 2026-11-06/07)

Row: doc/OpenWork.md "Is MHL the right radial rule, and where do its points go?".  Measurement only, no lib code.
Data, scripts and the (nR, L) tables: `~/Code/qchem6-runs/RadialRule/` (file list at the end).

## Method

- **One frozen D per system**, XC term only (LDA, `slater,vwn`).  Each candidate mesh is judged against a fine reference
  at that D: |ΔE_xc| against ε·N_atoms and max|ΔH_ij| (the XC matrix, AO frame) against ε, with **ε = 1e-5**.
- **N(ε) is persistent**: at each L the smallest nR that passes AND every larger nR on the ladder passes (a zero
  crossing is not a pass; several single rungs pass by luck, e.g. LiF `pe_m3_4.0_1.0` 100_53 at 4.8e-6 between
  1.6e-4 and 2.6e-5).  N is the smallest point count over the L's.  N(E) and N(H) are also given separately.
- **TeF₆** (PP and all-electron, 7 atoms, tolerance 7e-5 on E): a numpy Becke mesh plus PySCF `nr_rks` at the frozen D.
  It reproduces the C++ molecular mesh to 1e-14 in E_xc and 2e-15 in H (default row −4.4796e-5 / 1.122e-5 reproduced).
  References: (400, 77 GL).  Two independent references (m3 α2 and m2 α1, 400_77_gl vs 600_59) agree to 6.6e-8 in E
  and 4.7e-8 in H.
- **LiF** (trimmed VALENCE_LOWQ_VA, Li displaced 0.1 bohr, Γ, 2 atoms, tolerance 2e-5 on E): patched `rundeck`
  (`solid_scan.patch`, not committed), restart from the converged state (2 iterations), then `XCMeshAudit` per candidate.
  References: (300, 77 GL) at m2 α2 and at m3 α1.  **They agree to 2.8e-9 in E and 2.0e-9 in H.**
  *(Correction to my message of 11-06: the 5.8e-6 / 4.5e-7 I quoted there was a candidate against the RUN's mesh,
  not reference against reference.)*
- Defaults today: molecule m = 3, α = 2, (80, 29), Becke order 2 (`MoleculeSetup.C`); solid m = 2, α = 2, (40, 29),
  Becke order 3 (`MeshParams` default).  MnO was not run (time).

## 1. Radial rules, N(ε = 1e-5)

### TeF₆ PP (Lebedev, Becke order 2)

| rule | N(E) at | N(H) at | **N(both)** at |
|---|---|---|---|
| m2 α4 | 67,048 (50, 23) | 126,468 (60, 29) | **126,468** (60, 29) |
| m2 α2 | 67,192 (50, 23) | 126,540 (60, 29) | 126,540 (60, 29) |
| m2 α1 | 67,336 (50, 23) | 126,588 (60, 29) | 126,588 (60, 29) |
| m2 α0.5 | 67,456 (50, 23) | 147,752 (70, 29) | 147,752 (70, 29) |
| m3 α4 | 93,848 (70, 23) | 168,628 (80, 29) | 168,628 (80, 29) |
| **m3 α2 (today)** | 93,992 (70, 23) | 168,688 (80, 29) | **168,688** (80, 29) |
| m3 α1 | 80,676 (60, 23) | 168,748 (80, 29) | 168,748 (80, 29) |
| m3 α0.5 | 107,716 (80, 23) | 210,992 (100, 29) | 210,992 (100, 29) |
| Log (4 variants) | never | never | never |

### TeF₆ all-electron

| rule | N(E) at | N(H) at | **N(both)** at |
|---|---|---|---|
| m2 α0.5 | 126,648 (60, 29) | 211,664 (70, 35) | **211,664** (70, 35) |
| m2 α1 | 168,784 (80, 29) | 241,480 (80, 35) | 241,480 (80, 35) |
| m3 α4 | 168,628 (80, 29) | 300,460 (100, 35) | 300,460 (100, 35) |
| **m3 α2 (today)** | 240,832 (80, 35) | 301,212 (100, 35) | **301,212** (100, 35) |
| m2 α2 | 168,712 (80, 29) | 301,308 (100, 35) | 301,308 (100, 35) |
| m3 α1 | 210,932 (100, 29) | 301,404 (100, 35) | 301,404 (100, 35) |
| m3 α0.5 | 253,188 (120, 29) | 301,784 (100, 35) | 361,956 (120, 35) |
| m2 α4 | 252,936 (120, 29) | 451,040 (150, 35) | 451,040 (150, 35) |

The AE H floor at L = 35 is 4.32e-6 for every rule (the angular part).

### LiF (Becke order 3; Lebedev L 47 is built as 53, see §4)

| rule | N(E) at | **N(both)** at |
|---|---|---|
| m2 α4 | 29,480 (60, 35) | **130,416** (100, GL 47); 131,932 (80, L65) |
| m2 α2 | 33,604 (60, 35) | 149,276 (100, GL 47); 150,888 (80, L65) |
| m3 α2 | 41,578 (80, 35) | 174,478 (100, L65) |
| m2 α1 / α0.5, m3 α1 / α4 / α0.5 | 37.5k / 40.8k, 45.4k / 37.8k / 48.9k | no pass on the L ≤ 53 ladder (not run at 65) |
| Log (4 variants) | never | never |

**On LiF, H decides the mesh, and H is limited by the angular grid**: E meets ε at 30k points while H needs 130k.
**Today's solid default (40, 29) gives ΔE −1.1e-4 and max|ΔH| 1.8e-3**: it misses ε = 1e-5 by 5.6× in E and 180× in H.
At ε = 1e-4 the picks are m2 α4 59,500 (60, GL 41) and m2 α2 67,996 (60, GL 41).

**Verdict on m:** m = 2 beats m = 3 everywhere: by 25% on TeF₆ PP, by 14-25% on LiF (m2 α2 / α4 against m3 α2), and on TeF₆ AE only with a
small α (m2 α0.5 is 30% under today's m3 α2; m2 α2 is no better).  α matters at the 10-20% level and its best value
differs by system (PP: 4 ≈ 2 ≈ 1; AE: 0.5; LiF: 4), so no single α wins everywhere.

**Log (`RadialKind::Log`) never meets ε**, at any nR up to 300: TeF₆ AE best −9.1e-2 / 8.4e-3 at (300, 35);
PP −1.2e-2 / 6.8e-4; LiF −8.4e-3 / 3.1e-3 at (150, 53).  Its error falls only ×1/4 per doubling of nR (second order),
so it is not a candidate.

## 2. α per element (measurement patch)

TeF₆: α_Te ∈ {1, 2, 4, 8} × α_F ∈ {0.5, 1, 2}, m 2 and 3.  LiF: α_Li, α_F ∈ {0.5, 1, 2, 4} pairs, m 2 and 3, L 53.

- **TeF₆ PP: nothing beyond a rung**: the best is 126,474 (60, 29), the same as uniform m2.
- **TeF₆ AE: best 241,260 (m2, α_Te 1, α_F 2, (80, 35))**: equal to uniform m2 α1 and worse than uniform m2 α0.5.
  α_Te 8 fails AE.  (α_Te 0.5 was not in the grid.)
- **LiF: no per-element pair passes H at L 53** (the angular floor, §4); in E the best is m2 α_Li 4 / α_F 1 at
  49,882 (40, 53), against 29,480 uniform at L 35.
- **So per-element scaling buys nothing measurable** on these two systems.  The rung step (nR +10 to +20) is larger than
  anything α per element moves.

## 3. Becke order, which also explains the box row

| system, mesh | order 1 | **order 2** | order 3 | order 4 |
|---|---|---|---|---|
| TeF₆ PP (80, 29), max\|ΔH\| | 8.0e-5 | **7.4e-6** | 6.8e-5 | 4.4e-4 |
| TeF₆ PP (150, 35) | 1.8e-5 | **1.3e-6** | 1.7e-5 | 1.5e-4 |
| TeF₆ AE (80, 29) | 3.5e-4 | **1.1e-5** | 8.0e-5 | 4.7e-4 |
| LiF (150, 53) | not run | **1.7e-5** | 3.0e-5 | 1.4e-4 |
| LiF (200, 65) | not run | **3.1e-6** | 5.6e-6 | 3.0e-5 |

Order 2 (Becke's published three iterations) is best on both systems.  Order 3 costs ×9 in H on the molecule at the
same mesh, and ×1.8 on LiF.  **That is the box row's "9× molecule"**: lane A's TeF₆ box at (80, 29) gave 6.1-6.7e-5
at order 3 (`x80.out`, `x80m3.out`) and 7.4e-6 at order 2 (`x80m3o2.out`); the molecule here gives 6.8e-5 / 7.4e-6 at
the same orders.  Solids default to order 3 and molecules to 2, so the box and the molecule were never on the same
partition.

**The price on a solid is the build.**  The periodic Becke weights are an ε-converged image series (`beckeEps` 1e-6), and
a softer cell function reaches much further: LiF (200, 65) built in about 90 s at order 3 and in **65 min at order 2**;
(150, 53) took 33 min.  Order 1 was skipped on LiF for that reason.  So order 2 on solids needs the image series made
cheaper first; it is a cost question, not an accuracy one.

## 4. The angular rung on LiF

At nR 200 (m2 α2, order 3): Lebedev 35 gives 4.8e-4, 53 gives 3.0e-5, 65 gives 5.6e-6.  The C++ Lebedev table rounds
**41 and 47 up to 53** (974 directions) and 59 up to 65 (1454): `ang_41` = `ang_47` = `ang_53` bit for bit, and
`ang_59` = `ang_65`.  So the ladder skips the rung LiF needs.  Gauss-Legendre fills it: **GL 47 (1152 directions)
reaches 5.8e-6, which equals Lebedev 65 with 21% fewer directions and is 5× better than Lebedev 53**.  The GL ladder
is not monotone (GL 41 1.8e-5, GL 47 5.8e-6, GL 53 1.05e-5), so a persistent pick needs more than one L.

## 5. Where the points go, and dropping empty shells

Share of built points with |w ρ| > 1e-12: TeF₆ PP 51% (m3 α2 (80, 29)) and 63% (m2), AE 61% / 68%; LiF 66-85%.
On TeF₆ at m = 3, 36% of the shells lie beyond r = 10 bohr (r_max ≈ 1e6); at m = 2, 30% (r_max ≈ 1e4).  On a solid the
Becke cell has already removed the far tail (35% of nominal points), so the empty shells are innermost.

**Drop rule tested:** remove a whole radial shell of one atom when its charge Σ|w ρ| < τ, judged at the frozen D, split
into the shells **outside** the shell of peak charge and the shells **inside** it.

| system, rung | outer, τ = 1e-5: points kept / change in E, H | inner, τ = 1e-8: change in H |
|---|---|---|
| TeF₆ PP m3 α2 (80, 29) | 103,230 of 168,688 (−39%) / 6.0e-8, 2.8e-7 | 1.8e-7 |
| TeF₆ PP m2 α2 (60, 29) | 83,908 of 126,540 (−34%) / 4.6e-7, 9.7e-7 | 7.8e-8 |
| TeF₆ PP m2 α4 (60, 29) | 73,326 of 126,468 (−42%) / 4.5e-7, 1.0e-6 | 6.3e-8 |
| TeF₆ AE m3 α2 (100, 35) | 186,892 of 301,212 (−38%) / 6.3e-7, 1.6e-6 | 1.8e-7 |
| TeF₆ AE m2 α0.5 (70, 35) | 172,142 of 211,664 (−19%) / 9.6e-8, 2.5e-7 | 1.1e-14 |
| TeF₆ AE m3 α4 (100, 35) | 171,666 of 300,460 (−43%) / 1.4e-7, 3.7e-7 | 3.3e-8 |
| LiF m2 α2 (80, 65) | 133,400 of 150,888 (−12%) / 4.1e-8, 2.5e-7 | **2.5e-5** |
| LiF m2 α4 (80, 65) | 113,052 of 131,932 (−14%) / 9.3e-8, 6.1e-7 | **4.3e-5** |
| LiF m3 α2 (100, 65) | 159,614 of 174,478 (−9%) / 1.7e-8, 9.7e-8 | **8.5e-5** |

(Changes are against the same mesh undropped; every outer-dropped rung still passes against the reference.)

- **Outer drop at τ = ε is safe on all three systems**: E moves ≤ 6.3e-7 and H ≤ 1.9e-6 (ε/5) at τ up to 1e-5.
  τ = 1e-4 breaks it (H 1.1-4.3e-5).  It saves 34-43% on TeF₆ with m3 or m2 α ≥ 2, 19% where α is already small
  (m2 α0.5), and 9-14% on LiF.
- **Inner drop by charge is unsafe on the solid**: on LiF a shell carrying 1e-9 e moves H by 3.6e-6 and one carrying
  1e-8 e by 2.5-8.5e-5 (H/τ ≈ 3000), against H/τ ≈ 10 on TeF₆.  A shell close to a nucleus holds little charge
  because it holds little volume, while v_xc φ_i φ_j there can still be large.  **Charge is the wrong criterion
  for inner shells**; a criterion for them would have to bound the H element (|v_xc| φ_max² w per shell), which I did not test.
- **After the outer drop the rules move closer**: TeF₆ PP m2 α4 73k, m2 α2 84k, m3 α4 93k, m3 α2 103k (the m2 lead
  shrinks from 25% to 19-29%).  AE m3 α4 172k, m2 α0.5 172k, m2 α1 181k, m3 α2 187k, m2 α2 204k (all within 19%).
  Most of m = 2's lead on the molecule is tail points that a drop removes anyway.

**Caveats.**  (1) The drop uses ρ at a converged D; in production it would have to be decided before iteration 1, from
the seed (the seed-sizing row tests exactly that question for the rung).  (2) `RadialPoissonCoulomb` needs a rule that
stays uniform in its index; a dropped shell may conflict with a Poisson solve on the same mesh, which I did not check.
(3) Two systems, one solid, LDA, one D each; MnO (d shell, spin) not run.

## What I would propose (for the user to rule; no code written)

1. **Drop outer shells by charge at τ = ε, keep every inner shell.**  This is the largest and safest lever measured:
   34-43% of the molecular mesh, 9-14% of a solid's, at ≤ ε/5 in H.  It needs a density before iteration 1 (seed-sizing row).
2. **m = 2 on the molecule route too**, with α stated per system class until more chemistry is in (PP 2-4, AE 0.5-1).
   Worth 25-30% before the drop, about 20% after it.
3. **Not worth pursuing**: Log, and α per element.
4. **Two OpenWork rows** (filed with this report): the Lebedev table gap (41, 47, 59 round up; LiF needs the 47-59 rung,
   and GL 47 is 21% cheaper than Lebedev 65 at the same H), and order 2 on solids, gated on making the periodic Becke image
   series cheaper (×40 build).  Today's solid default (40, 29) misses ε = 1e-5 by 180× in H on LiF.  That belongs to the
   seed-sizing row, which will set solid rungs from ε.

## Files (`~/Code/qchem6-runs/RadialRule/`)

- TeF₆ tables: `TeF6_{pp,ae}_{check,uniform,perelement}.tsv` (rule nR L pts pts_q dE maxdH dN sec), `order_{pp,ae}.log`,
  shells `TeF6_{pp,ae}_shells.tsv`, drops `TeF6_{pp,ae}_drop{,2,3}.tsv` (`drop3` = the inner/outer split).
- LiF: `LiF_table.tsv` (all 486 candidates vs the reference), `LiF_drop.tsv`, `LiF_<mesh>_shells.tsv`, raw `solid/*.bin`
  ([E, n, H n×n]) and `solid/*.mesh` (x y z w ρ site), candidate lists `solid/scan*.json`.
- Scripts: `rmesh.py`, `setup_sys.py`, `sweep.py`, `ntable2.py`, `order.py`, `step2.py`, `step2b.py`, `step2c.py`,
  `lif_table.py`, `lif_shells.py`, `drop_sum.py`; the C++ measurement patch `solid_scan.patch` (snapshot
  `snap/laneb-rr-20261106-085824`).
