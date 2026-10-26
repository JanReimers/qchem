# Sharp terms and their grid requirements — which Hamiltonian terms carry a Gaussian/erf feature, and what each needs from a grid

Written 2026-10-25 (user: *"list which potential terms have sharp features (gaussians or erf) that must be considered
when setting up grids"*).  Companion to `GPWGrids.md` (the 2026-07 inventory of every grid and what sizes it): that
record says WHICH grids exist; this one says WHICH TERMS put sharp features on them and WHAT EACH REQUIRES.  It is a
reference, not a queue — the open decision it leads to is the OpenWork §3 row "GPW `BallOnly` raster + the C = 2
density floor ALIAS tight pairs".  Evidence: UTLattice_BS `GPW.TightPrimitiveVlongOracle` /
`GPW.TightPrimitiveCollocationOracle`, the grid audit (`grids.auditDensity`, `grids.auditVlocLong`, commit d3a487e1),
and `HubbardUHistory.md` 2026-10-25.

## 1. The terms

| term | sharp feature | exponent that matters | how it is evaluated today | grid requirement | audited? |
|---|---|---|---|---|---|
| **Overlap, kinetic** | basis pairs | \(p=\alpha_i+\alpha_j\) | analytic lattice sums | none | n/a |
| **Local PP, short** (GTH \(C_1\ldots C_4\,e^{-r^2/2r_{loc}^2}\)) | Gaussian × polynomial at each ion | \(\beta_s=1/2r_{loc}^2\) per species (O 8.15, Mn 1.22); integrand \(p+\beta_s\) | analytic lattice sum (`MakeLocalGaussian`) | none on the production path; the κ-sweep FALLBACK (no closed Gaussian form) is grid-based | no |
| **Local PP, long** (\(-Z\,\mathrm{erf}(r/\sqrt2 r_{loc})/r\)) | erf = the ion's Gaussian-smeared charge | \(\beta=1/2r_{loc}^2\) of the SHARPEST species | field built in G space, then **integrated back on a raster** against each pair | BALL: harmonic rule on \(p\beta/(p+\beta)\); RASTER: the PAIR \(p\) | **yes** (`auditVlocLong`) |
| **Separable (KB) nonlocal PP** | projectors \(r^{l+2(i-1)}e^{-r^2/2r_l^2}\) per channel | \(1/2r_l^2\) per channel (Mn d the sharp one) | analytic overlaps \(\langle\chi|\beta\rangle\) | none; the mesh FALLBACK has its own floor (V2.5) | no |
| **Hartree** | the density = products of basis primitives | \(p_{max}=2\alpha_{max}\) | pairs **collocated** onto a raster → FFT → \(4\pi/G^2\) → **integrated back** | RASTER resolves \(p\) (alias); ball judged loosely (the \(1/G^2\) kernel) | **yes** (`auditDensity`) |
| **XC, uniform route** | \(\rho\) and \(v_{xc}\sim\rho^{1/3}\), peaked at the nuclei | pairs \(p\) plus the field sharpness fieldSharp·\(\alpha_{max}\) = \(\tfrac23\alpha_{max}\) (`routingBeta`) | \(\rho\) sampled pointwise on the density raster (raw collocation, exact at the points); \(v_{xc}\)·pair integrated back on it | same raster as Hartree, plus the field routing | covered by `auditDensity` |
| **XC, Becke route** (β-MnO₂'s) | same \(\rho\)/\(v_{xc}\) features | \(\alpha_{max}\) (radial), angular content | atom-centred mesh (nR 40, L 29), OFF the raster | radial + angular convergence of the mesh | no (2026-10-24: mesh ×3 moved nothing) |
| **Kerker/Pulay mixing, cusp-deficit map** | the mixed \(\tilde\rho\) on the density ball | \(p_{max}\) | G space, on the density grid | the density ball | indirectly |
| **+U projectors** | pseudo-atomic radials (16-exponent pools: SHARP) | the pool's largest exponent | analytic overlaps (GPW); closed-form transform + accumulated Parseval (PW) | none (PW: the enlarged sphere for O, converged to 1e-12) | n/a |
| **Ion–ion (Ewald), \(E_{\alpha Z}\) alignment** | the Gaussian Ewald split | the Ewald splitting parameter | analytic, basis-independent | its own convergence; no basis grid | n/a |
| **Linear-response kernels** (Hartree / \(f_{xc}\) response) | as Hartree / XC | as above | same rasters and meshes as the ground state | as Hartree / XC | via the same rasters |

**Two of these share one raster but need different things from it**: Hartree needs the pair resolved where \(4\pi/G^2\)
weights it (a loose ball, a non-aliasing raster); V_loc long needs a ball sized by the harmonic pair × field rule plus
the same non-aliasing raster.  XC on the uniform route adds a third need: a field as sharp as \(\tfrac23\alpha_{max}\),
which is why the routing floors pairs at that sharpness.

**On the plane-wave basis none of the sampling issue exists**: every matrix element is a sum over the exact
\(G-G'\) difference set; only the XC raster (4× the orbital cutoff) is a real-space grid.

## 2. The two requirements of a raster-integrated term

A term that integrates (or collocates) Gaussian PAIRS on a raster has two independent requirements:

- **BALL** — its G-space sum must reach where the integrand's transform dies.  A pair \(p\) against a field of
  Gaussian exponent \(\beta\) sums \(e^{-G^2(1/4p+1/4\beta)}\): the harmonic exponent \(p_h=p\beta/(p+\beta)\)
  (\(\beta=0\): no field, \(p_h=p\)); the tail at the ball edge \(G_b\) is \(e^{-G_b^2/4p_h}\).
- **RASTER** — the pair is SAMPLED in real space, so its spectrum beyond the raster's Nyquist \(G_N\) aliases back.
  The first image sits at \(2G_N\): relative amplitude \(e^{-(2G_N)^2/4p}=e^{-G_N^2/p}\).  `BallOnly` puts \(G_N\) at the
  ball edge (\(G_N^2=2E_{cut}\), so \(\varepsilon\approx e^{-2E_{cut}/p}\)); `AliasFree` doubles it
  (\(\varepsilon\approx e^{-8E_{cut}/p}\)).  \(G_N\) = the raster's inscribed Nyquist radius \(\min_i \pi N_i/|a_i|\).

The tightest pair (\(p=2\alpha_{max}\)) is routed to the finest level, so the finest raster is the one that must pass.

## 3. What was measured (2026-10-25)

One Mn q7 (pure erf local) in a 10-bohr box, single primitives, against closed forms at any exponent:

| \(E_{cut}/p\) | `BallOnly` d-diagonal \(V_{long}\) error | `BallOnly` collocation (low-G \(\tilde\rho\), rel.) |
|---|---|---|
| ≈ 1–2 (the auto floor + its completion rung) | ±1.5–2.7 Ha | 0.04–0.16 |
| ≈ 2.7 | 1e-2 – 5e-2 Ha | — |
| ≥ 5.5 | ≤ 1e-5 Ha | ≤ 6e-8 |
| `AliasFree` at the auto floor (C = 2) | ≤ 3e-4 Ha | ≤ 4e-5 |

So every SIZING rule follows \(\alpha\) (density floor, the V_loc-long harmonic ball, the completion rung at
\(2E_{fine}\)); what fails is the real-space SAMPLING when `BallOnly` is paired with C = 2.  The error is
non-variational: one Mn d primitive at \(\alpha=76.7\) made the β-MnO₂ SCF "converge" at −6917 Ha (Eloc −8764 Ha).

A second, latent condition: the auto floor looks only at the BASIS.  A local-PP field sharper than the basis
(\(\alpha=2\) alone against Mn's erf, whose content reaches \(|G|\approx12\)) is off by 0.03 Ha even alias-free.  In real
bases \(\alpha_{max}\gg\beta_{PP}\) (β-MnO₂: 36 vs 8.15), so it does not bite there.

## 4. How it got this way

`GPWGrids.md` (2026-07-21) records C = 8 = 4 × 2: "resolve a Gaussian of exponent p ... needs Ecut ≈ 4p, measured by the
negCharge probe" × the density's product exponent — and an ALIAS-FREE raster by construction.  Since then C dropped
8 → 2 (doc/GPWPlan 0.5(f2): the raw collocated XC feed removed the negative-lobe problem C = 8 compensated), justified on
NaF's total-energy plateau (sub-2 mHa from C ≈ 1.5), and `BallOnly` became the production raster ("CP2K's bet, the
fold-back controlled by Ecut").  Each step was right on its own premise; together they dropped the ≈ 4p requirement,
and NaF (no tight d) could not show it.  (Pin: fit quality is the grid convergence of ρ, never ΔE_total.)

## 5. The audit

Every raster-integrated term reports one row as it assembles (`grids.auditDensity` at the grid build,
`grids.auditVlocLong` at the V_loc-long assembly): the tightest pair \(p\), \(\beta_{field}\), the ball's Ecut / needed G /
present G / tail, the raster N / needed \(G_N\) / present \(G_N\) / alias, and OK/UNDER.  UNDER at alias > 1e-6; the ball
is judged against the tolerance it was built to (V_loc long: `vlocEps`), Hartree's ball is reported but not judged.
Every UNDER also prints a `[grid audit] UNDER:` console line.  β-MnO₂ baseline (VA + O p 0.18): `BallOnly` → density
45³ alias 1.8e-2, V_loc long 50³ alias 7.1e-3 (both UNDER); `AliasFree` → 1.1e-7 / 2.5e-9 (OK).

Not audited: the κ-sweep and KB-mesh fallbacks, and the Becke mesh (a different kind of grid — radial/angular
convergence, not alias).
