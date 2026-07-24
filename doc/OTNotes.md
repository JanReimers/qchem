# OT (Orbital Transformation) notes — durable findings for the OT increment

CP2K's OT = a preconditioned direct minimiser (same *role* as GDM in this codebase:
`doc/SCFStrategyPlan.md` increment 5).  These notes consolidate what the 2026-07-23/24 GDM
investigation established, so the OT work does not re-derive them.

## 1. E[ρ] IS variational on the GPW collocation path — OT should work
The headline worry ("the fitted GPW energy is non-variational, so a direct minimiser can't descend
it") is **FALSE**.  GDM converges GPW cleanly:
- NaF/Γ, production grid (Ecut=80, BallOnly, raw-XC): DIIS→GDM ladder → **−24.4304 in 22 iters** (the
  CP2K-matching anchor), every geodesic step a descent, gap healthy throughout.
- The apparent "limit cycle" was an artifact (see §3), not physics.

So OT (preconditioned direct-min) descends the same energy.  There is **no non-variational wall** to
design around — provided the potential is built as the *exact adjoint* of the energy quadrature (§2).

## 2. Use the RAW-adjoint V_xc, never the ball-fit — variationality ≠ fit quality
`PW_XC::CalcMatrix` (src/Hamiltonian/Internal/Imp/PWTerms.C) has two routes:
- **RAW** (`applyRawAdjoint`, taken when `itsRhoIsRaw`): the exact discrete adjoint of the `applyRaw`
  collocation → `V_xc = ∂E_xc/∂D` to **machine precision at any grid**.  Gate: `GPW.RawXCConsistencyFD`.
  **This is the route a direct minimiser needs.**
- **BALL** (`DoFit`/`Overlap`, the fallback for a matrix-free seed): `⟨i|v_xc_fit|j⟩` — exact *given the
  fit coefficients*, but the fit of the **nonlinear** `v_xc = ρ^{1/3}` onto the density's `{G}` ball
  (the `relCutoff ≤ 1` clamp) leaves a residual `⟨i|v_xc_fit − v_xc|j⟩`.  FD-measured non-variationality:
  **rel 1.93 under BallOnly vs 6.6e-7 under AliasFree** (`GPW.XCPotentialConsistencyFD`).

Key distinction for OT: **variationality (V = ∂E/∂D exactly) is a different axis from fit accuracy.**
The raw route is variational *without* a good fit; the ball route only reaches it in the complete-fit
limit.  Hartree is variational on either route (linear in ρ → consistent ball truncation; the G-space
Coulomb metric `4π/G²` is already Dunlap-robust, so **no molecular repulsion-metric Dunlap is needed for
GPW Hartree** — that thread is closed).

## 3. Direct-min must DISABLE the mixer, and bail to a MIXED step — not an unmixed one
Two coupled facts a direct minimiser (GDM/OT) must respect:
- **When it owns the density, the mixer is OFF** — per-step (the driver never calls `Mix`) AND post-step
  (the adaptive `WantsReDamp`/`UpdateRelax` must be skipped, or a LinearMixer re-mixes *after* a geodesic
  step and corrupts it).  In-code: one `lineSearch` bool gates the driver, the mixer bypass, and the
  honest `ρ_mix="----"` display (src/SCFIterator/Imp/SCFIterator.C).
- **When the geodesic can't engage, fall back to a MIXED fixed-point step, NEVER an unmixed diagonalize.**
  An unmixed diagonalize runs away on ill-conditioned systems — this is what masqueraded as
  non-variationality: a GDM gated off by too-tight `FDMax` degraded to unmixed diagonalize.  Fixed via
  `tSCFAccelerator::CanLineSearch()` (readiness = seeded + `[F,D] < FDMax`): the iterator runs direct-min
  only when `WantsLineSearch() && CanLineSearch()`, else a mixed step.  `α=1` (LinearMixer passthrough)
  keeps molecular direct-min byte-identical.

## 4. Open design seam for OT — the mixer off-switch (stale-on-resume)
Today the mixer is *bypassed* (not told), so its state (`KerkerMixer::itsMixedRho`) **freezes** through a
direct-min stretch and is **stale if mixing later resumes**.  Doesn't bite now (once GDM engages with a
wide `FDMax` it converges without falling back), but OT — which spends long runs in direct-min — will hit
it.  Options (OOD call, user 2026-07-23): an explicit `Pause()`/`Resume()` on `tDensityMixer` (Tell), or
swap the active mixer for a null mixer during OT.  Do NOT add an `IsUsed()` query to the mixer — "unused"
is a driver property, not a mixer property.

## 5. Scheduling — hand off to OT on ΔE for solids, [F,D] for molecules
`tSCFAcceleratorLadder` takes `ScheduleSignal { Error, EnergyChange }` (the tail hand-off gates on
`[F,D] < switchat` or `|ΔE/E| < switchat`).  For solids, `[F,D]` is contaminated (charge-transfer slosh,
diffuse-mode ill-conditioning, giant-response spikes) and may never settle; the total energy (a
variational scalar) does.  Measured on NaF: energy-gated hand-off → **19 iters vs 44** for `[F,D]` (OT
starts closer to the minimum).  The signal choice is a loop-face policy, owned once by the ladder.

## 6. Naming + engage threshold
- The DIIS/GDM residual **is** `[F,D]`, not the energy — the params are `FDMax`/`FDMin` (renamed from the
  deceptive `EMax`/`EMin`).  (JSON/facade/CLI still say "EMax" — a separate legacy surface, item-1 graduation.)
- The engage threshold (`GDMParams::FDMax`) must be **wide enough that the minimiser reliably engages once
  handed off**, or it degrades to the §3 fallback every step.  For the ladder tail, effectively "always
  engage once handed off" (e.g. `FDMax = 1.0`).

## 7. Debug instruments (env-gated, off by default)
- `GPW_XCROUTE` — logs the V_xc route (RAW vs BALL) + ρ grid range per iteration (PWTerms.C).
- `GPW_GDMTRACE` — logs the geodesic line search DESCENT vs FALLBACK + best margin (SCFIterator.C).
- `GPW_GRIDCHARGE` — ∫ρ_grid vs Tr(DS) + ρ stats.
- Note: collocated ρ dips slightly negative (~−5e-3) even for a PSD density matrix — normal grid
  band-limiting (Gibbs), unlike molecular pointwise ρ; the `ρ>0` guard kink is too small to stop the
  minimiser, but a strictly-PSD collocation would remove it.

## 8. Δρ metric is mixer/driver-dependent — unify with the OT mixer-seam work
The per-iteration **Δρ column / MinΔρ gate is NOT one consistent quantity** (found 2026-07-24):
- `LinearMixer` (molecular) and the direct-min (GDM/OT) driver return `‖ΔD‖_maxabs / N` — REAL-space,
  relative to electron count.
- `KerkerMixer` / `PulayMixer` (the solid DIIS phase) return `‖Δρ̃‖_∞` — a G-SPACE, ABSOLUTE residual.

So a DIIS→GDM run's Δρ column FLIPS meaning at the hand-off (absolute G-space → relative real-space), which
is why the GDM tail Δρ "looks like noise" relative to the DIIS rows.  The `ΔE` column, by contrast, is
uniformly relative (`ΔE/|E|`) — its label was fixed to `ΔE/E` (2026-07-24); the `Δρ` label was left as-is
pending this unification.  **Punted to the OT project** (user 2026-07-24, "not too concerned if it switches
real/G space"): when the mixer off-switch / null-mixer seam lands (§4), unify the residual to ONE definition
(recommend relative `‖Δρ‖/N`), relabel `Δρ/N`, and re-pin the affected MinΔρ anchors.  `ReDampMix` currently
returns an un-normalised change too (a third variant) — fold it in.

## 9. Metals / smearing (from SCFStrategyPlan)
OT + Fermi smearing is CP2K's metals answer; the occupation-face already carries fractional occ.  When
smearing lands, the schedule signal becomes the free energy `A = E − TS` (Mermin), and OT minimises that.
See `doc/SCFStrategyPlan.md` and `doc/GPWPlan1.md` item 4b.
