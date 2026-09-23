# Self-consistent orbital-resolved (U, J) from our own on-site ERIs — the plan

**Born 2026-09-23**, at the top of `doc/` because it is being executed (`doc/OpenWork.md` §1 step 5,
increment 3 remainders).  Retires to `doc/Records/` the day the U-functional lands.  The north star above
it is `doc/BatteryMaterialsRoadmap.md`; the ruling under it is `doc/Pins.md` pin 23.

★ **What is new here, and worth saying once.**  Orbital-resolved DFT+U, determined self-consistently from
the code's *own* on-site ERIs, in a **Gaussian** basis, on transition-metal oxides, on a 14 GB desktop.
Each piece exists somewhere; the combination does not.  ACBN0 was built on plane waves and had to invent a
PAO-3G projection to get four-index on-site integrals at all — we own them natively (`Cache4`/Rk).  CP2K
has the Gaussians and hand-sets U.  QE has the U determination and is plane-wave.  Nobody has the
orbital-resolved, self-contained, Gaussian combination, and the reason to care is not novelty: it is that
a supercell-free U is the only kind that scales to a **composition sweep**, which is what a voltage curve
is.

---

## 1. Background — where we stand (2026-09-22)

**Built and banked** (increments 1–3, `doc/OpenWork.md` §1 step 5): the `Hubbard_U` term (scalar-generic,
spin-native, real-TRIM corner); the manifold as an **input** list of (site, shell, irrep, U) with no code
path assuming the Hubbard atom is the transition metal (pin 23); `ManifoldSymmetry` slots from the site's
coordination point group; three projector flavours (Löwdin column, contracted **atomic**, **ortho-atomic**
with U=0 spectators) matching CP2K's and QE's conventions; `BasisSet::BareCoulombSource` + `ERI4Block` for
the on-site bare integrals; the `ACBN0` estimator; and `SolidCalculation::ConvergeHubbardU`, the outer loop.

**Measured.**  On MnO and NiO, with the same UPF and projector as the oracle:

| | our ACBN0 | seed | ACBN0 paper | hp.x (linear response) | our occupations vs hp.x |
|---|---|---|---|---|---|
| MnO (d⁵) U(Mn 3d) | 10.77 eV | IonicSAD | 4.67 | 0.96 ⚠ | 4.93/0.45 vs 4.988/0.481 |
| NiO (d⁸) U(Ni 3d), ortho | **14.51–14.76** | IonicSAD | 7.63 | **5.27** | 4.920/3.445 vs 4.980/3.368 |
| NiO (d⁸) U(Ni 3d), full ortho set | **13.80–13.99** | IonicSAD | 7.63 | **5.27** | 4.894/3.204 |
| NiO U(O 2p) | 7.27 | IonicSAD | 3.0 | — | |

★ **The seed column exists because the seed CHANGED under these numbers** (`NiOSpec` went SAD → IonicSAD on
2026-09-23 when the "Ni²⁺ cannot be seeded" claim was retracted), and a table that cannot say which
configuration produced it is the failure mode the VA/VB span files exist to prevent.  **Re-measured on both
seeds, the answer does not move**: ortho 14.60/14.67 (SAD, 63 iterations, E = −109.2322219) vs 14.51/14.76
(IonicSAD, 114 iterations, E = −109.2322528) — 31 µHa apart, U_eff within 0.7 %, occupations to four
decimals; the full ortho set agrees to 0.7 % on Ni 3d and to **0.02 %** on O 2p.  ⇒ the banked numbers are a
property of the converged state, not of the guess that found it.  ⚠ The IonicSAD ortho run needed **114**
iterations where SAD needed 63, and reported "NOT converged" at the old 80-iteration default — **the third
time an iteration cap produced a wrong conclusion this week**; the probe default is now 200.

★★ **The finding that matters: the manifold is right and the functional is wrong.**  Our projected
occupations agree with hp.x to 1–2 % on the same cell, same pseudopotential, same projector — so the
factor 2.6 in the *value* is the U functional alone, not the projector, manifold or PP.

⛔ **TWO DIFFERENT COMPARISONS LIVE IN THIS TABLE AND THEY ANSWER DIFFERENT QUESTIONS.  Never merge them.**
- **ours ÷ published ACBN0** (1.8–2.8 across MnO-d, MnO-O2p, NiO-d, NiO-O2p) measures **PROJECTOR
  COMPLETENESS** — the paper's PAO-3G-projected plane-wave states carry ~60 % of their norm and ours carry
  ~100 %.  It is *not* evidence about screening, and a kernel fitted to close it would be fitting out a
  basis artefact.
- **ours ÷ an INDEPENDENT oracle** (hp.x, cRPA — never ACBN0) is the only thing that can test the
  **SCREENING** hypothesis.  Gate 3 holds that comparison, and an earlier version of gate 3 mixed the two,
  which is the circularity that had to be retracted on 2026-09-23.
**ACBN0's \f$\bar N^2\f$ renormalisation is not a screening model** — it vanishes as the basis approaches
completeness, which is exactly the regime a Gaussian code lives in.

**⚠ Three traps already paid for** (do not re-pay them):
1. **hp.x's MnO number is not usable** — d⁵ high-spin is filled-majority/empty-minority, \f$\chi_0=-0.045\to\chi=-0.043\f$
   (96 % of the bare response survives), so \f$\chi_0^{-1}-\chi^{-1}\f$ nearly cancels.  NiO (d⁸) is the healthy
   oracle: \f$\chi_0=-0.113\to\chi=-0.071\f$.
2. **An hp.x number is conditioned on its starting U.**  The NiO 5.267 eV is \f$U_{\rm LR}(U_{\rm in}=3\,{\rm eV})\f$
   (the decks inherit `U Ni-3d 3.0` from QE's benchmark); the MnO decks carry none.  Read the `HUBBARD` block
   before comparing.  A response is a function of the state it linearises about.
3. **LDA NiO at U=0 loses AFM-II — and it is PHYSICS, not seeding** (strengthened 2026-09-23).  Three runs
   now reach the same non-magnetic fixed point by different paths: SAD seed 79 iterations to
   E = −109.2693031, SAD+MOM *bit-identical* (so MOM holds the collapse rather than preventing it), and
   **IonicSAD 48 iterations to −109.2693702** — 67 µHa apart, same state, collapse signature
   \f$N\uparrow/N\downarrow\f$ = 4.2889/4.2889 in both.  ⇒ **two different seeds, different trajectories, one
   fixed point**: LDA NiO at Γ on this cell genuinely has the non-magnetic solution as its SCF fixed point,
   so gate 1 is asking a physics question and no amount of seeding will answer it.  The
   **ACBN0 outer loop is a non-result on NiO** — sites
   decouple to \f$U_{\rm eff}=-1.05\f$ eV, Hartree runs to 2.04× its floor.  MnO's clean monotone 8-step loop
   was a property of d⁵, not of the loop.  ⇒ **magnetic robustness under a changing U is a gate, not a detail.**

---

## 2. The decision, and why the APPLICATION decides it

Four routes to a U, not three — the fourth is `doc/OpenWork.md` step 5 item 5's parked fallback:

| | route | needs a supercell? | status |
|---|---|---|---|
| (a) | ACBN0 as is | **no** | ⛔ REFUTED — overshoots on both materials and both manifolds; its loop drives U *away* from the literature (MnO 10.8 → 16.8) |
| (b) | **ACBN0 with a screened interaction** | **no** | ★ THE ROUTE |
| (c) | hp.x values as input | **yes** (q-mesh ≡ supercell) | a bridge: legitimate under pin 12 (computed input ≠ hand-set knob), but one QE run per material *and per composition* |
| (d) | our own finite-difference linear response | **yes** | parked on cost; scriptable over the projector we already own |

★★ **The supercell column is the whole argument, and it is a TRANSLATIONAL-SYMMETRY statement** (user
asked, 2026-09-23).  (c) and (d) measure a *response* to \f$\alpha\hat P_I\f$ — a potential shift on ONE
site's manifold — and want the matrix \f$\chi_{IJ}=dn_I/d\alpha_J\f$ including \f$J\neq I\f$, because U
comes from a difference of INVERSES.  Under PBC the perturbation replicates coherently onto every equivalent
site, so a primitive cell returns \f$\sum_J\chi_{IJ}\f$ (the q=0 component) — one number where a matrix was
needed.  An isolated perturbation is genuinely not translationally symmetric, and there are exactly two
cures.  **Supercell**: break it explicitly, to a period long enough that the response to the nearest image
has decayed (the classic 2×2×2).  **DFPT**: never break it — perturb monochromatically with
\f$\alpha\hat Pe^{i\mathbf q\cdot\mathbf R}\f$, which IS Bloch-periodic and lives in the primitive cell,
then Fourier-sum \f$\chi_{IJ}(\mathbf R)=N_q^{-1}\sum_{\mathbf q}\chi_{IJ}(\mathbf q)e^{-i\mathbf q\cdot\mathbf R}\f$.
An \f$N_q\f$ mesh is MATHEMATICALLY EQUIVALENT to an \f$N_q\f$ supercell — the same dichotomy as frozen-phonon
vs DFPT phonons, which is what DFPT was invented for.  ⇒ our own `nq=1` non-result was q=Γ only, i.e. all
sites perturbed IN PHASE (`IntegrationTests/QE/README.md`: *"the q-mesh is the supercell size of the
linear-response method"*).  **Either cure costs a supercell-equivalent per composition.**

★ **(a) and (b) PERTURB NOTHING.**  They evaluate on-site integrals on the ground-state density in whatever
cell you already ran — no perturbation, no images, nothing to break, so the question never arises.  NiO's U
came out of a 4-atom cell in ~10 minutes.  A voltage curve needs U at many compositions; only the ACBN0
family delivers that without a supercell each time.  (a) is refuted ⇒ **the application selects (b)**, and
it selects it independently of the physics argument.

**What (b) actually requires, honestly.**  The bare \f$F^0\approx27\f$ eV is right (it is a bare integral on
the atomic Slater scale) and \f$\bar N^2\approx0.64\f$ measures basis completeness, so the missing physics is
the dielectric screening \f$\varepsilon^{-1}\f$.  ⛔ **\f$\varepsilon\f$ must be COMPUTED, never entered** — a
hand-set \f$\lambda\f$ or \f$\varepsilon\f$ is pin 12 one layer down, and the "~1/4 for a TMO" that motivated
this route was a *sizing* argument (does the gap have the right magnitude to be screening?), not a proposed
value.  Three ways to get it, in increasing order of cost and correctness:
- **Thomas–Fermi** from the valence density: no empirical input, but it is a *metallic* screening model and
  these are insulators — it will over-screen.  Cheap; expect it to be wrong, and measure *how* wrong.
- **\f$\varepsilon_\infty\f$ from our own response** to a field.
- **cRPA**: the principled route, the one the literature uses for screened U, and a solver increment.

★ **CP2K is the reference implementation for the machinery, not for the value.**  `dft_plus_u.F` only
*applies* U (`u_ramping` ramps toward a hand-set `u_minus_j_target` as an SCF aid — not a determination),
so there is nothing to learn there about values.  But `qs_linres_*.F`, `response_solver.F` and especially
`qs_linres_polar_utils.F` ("Polarizability calculation by dfpt … Berry phase operator … periodic Raman")
are **periodic DFPT response in a Gaussian basis**.  That is exactly the "solver increment we have no seam
for", implemented in our own basis type rather than QE's plane waves — worth reading before we conclude
anything about what a computed \f$\varepsilon\f$ costs us.

⚠ **The honest size of (b).**  ACBN0's renormalisation *is* the paper's screening model; replacing it with
a real screened kernel is not tuning ACBN0, it is substituting a different mechanism and keeping the
bookkeeping.  If the screened kernel ends up needing cRPA, ask explicitly whether we are building cRPA with
extra steps — and if so, whether (d) is cheaper after all.  **Step 3 below is designed to answer that
before we commit.**

---

## 3. The target: a Li\f$_x\f$Mn\f$_2\f$O\f$_4\f$ voltage curve without a U per configuration

The user's scheme: compute (U, J) at the three ordered compositions — **λ-MnO₂** (Mn⁴⁺, d³),
**LiMn₂O₄** (nominally Mn³·⁵⁺), **Li₂Mn₂O₄** (Mn³⁺, d⁴) — with no supercells, and use those for every
Li\f$_x\f$Mn\f$_2\f$O\f$_4\f$ configuration the cluster expansion needs.  Feasible exactly because (b) is
supercell-free (§2).

★ **One refinement, and it falls out of what is already built.**  Do not interpolate U on the global
composition \f$x\f$ — **assign U per Mn SITE by its local oxidation state**.  U is a property of the local
electronic configuration, not of a composition average, and pin 23 already makes the term take a per-site
list, so this needs no new capability:
- λ-MnO₂ → \f$U({\rm Mn}^{4+}), J({\rm Mn}^{4+})\f$
- Li₂Mn₂O₄ → \f$U({\rm Mn}^{3+}), J({\rm Mn}^{3+})\f$
- LiMn₂O₄ → **a transferability CHECK, not a third interpolation node**: if it charge-orders into Mn³⁺ +
  Mn⁴⁺, the per-site U's it produces must reproduce the two end members' values.  That is a falsifiable
  prediction, which an interpolation on \f$x\f$ is not.

★★ **AND U(O 2p) AT ALL THREE COMPOSITIONS, AS A FIRST-CLASS TARGET — not a spectator** (user, 2026-09-23).
Three reasons it is on the critical path for the VOLTAGE specifically, not just for the band structure:
(1) pin 23 exists because β-MnO₂'s decisive correction was on **O-p_z, not Mn-d**; (2) we measured
U(O 2p) = **7.27 eV** on NiO — against the paper's 3.0, and nowhere near negligible; (3) delithiation is
formally Mn³⁺→Mn⁴⁺ but a real fraction of the hole lands on **O 2p (ligand hole / oxygen redox)**, so
U(O 2p) moves the computed voltage directly.  The O sublattice also stops being equivalent once Li is
partially removed.  ★ **This costs nothing extra**: the `orthofull` arm already gives every listed manifold
its own estimate from ONE SCF — the NiO run returned Ni 3d, Ni 4s, O 2s and O 2p together.  List O 2p (and
O 2s) at U=0 and read all of them off each of the three runs.

Then each Mn in a CE training supercell takes its U from its own local environment, and **the CE training
runs need no new U calculations at all**.

**Risks, in the order they are likely to bite:**
1. ⛔ **Magnetic robustness — and it is NOT an argument for non-collinear.**  Ruled 2026-09-23 (user):
   **collinear is enough for a room-temperature voltage curve**; non-collinear order is not on this path.
   The risk that remains is narrower and real: the spinel Mn sublattice is the **pyrochlore lattice —
   geometrically frustrated**, so there are many near-degenerate COLLINEAR states, and an SCF can slide
   between them.  Two ways that bites: (a) the imposed order collapses under a U change, as it did on NiO
   (trap 3) — and NiO was an *unfrustrated* rocksalt AFM, so this is strictly harder; (b) worse for our
   purpose, different Li configurations land in DIFFERENT magnetic states, which pollutes exactly the
   energy DIFFERENCES the cluster expansion is fitted to.  ⇒ what matters is that the magnetic state be
   **consistent and reproducible across runs**, not that it be the true ground state.  **Gate 1 below
   exists for this and nothing else.**
2. **Mn³⁺ d⁴ high-spin is Jahn–Teller active** (e_g¹) — it is the whole story of LiMn₂O₄'s structural
   transition, and Li₂Mn₂O₄ is tetragonally distorted because of it.  Orbital-resolved U on e_g is both the
   best test of orbital resolution and the most dangerous place to apply it: Macke's FeS₂ warning is that
   correcting a hybridised e_g wrecked the lattice parameter.  Do U at **fixed geometry** first; U ↔ JT
   distortion is a coupled problem and must not be entered accidentally.
3. **LiMn₂O₄ charge ordering is delicate at the DFT level.**  If it does not charge-order, the
   transferability check is weaker — say so rather than reading the average as a third node.
4. **O 2p is not a spectator by assumption.**  We measured 7.27 eV on NiO (paper 3.0), and pin 23 exists
   because β-MnO₂'s decisive correction was on O-p_z.  Carry O 2p in the manifold list at U=0 and *measure* it.

---

## 4. The work, in order

Gates 1–3 are cheap enough for now-to-Oct-5; the long unattended runs are sized for the **Oct 6–20** window
(user away).  ⛔ Unattended runs go through **`scripts/memsafe -p`** (cgroup + OOM shield), never bare.

**Prerequisites (no physics, do them first).**
- **A Li valence basis, and q1 vs q3 is a REAL TEST, not a formality** (user, 2026-09-23).  There is no
  `LI` block in any `valence_lowq_*.bsd`.  GTH LDA offers Li **q1** (2s¹ only; 1s frozen into the core) and
  **q3** (1s²2s¹ explicit).  The tension is specific to a cathode: Li is nearly fully ionised, so q1's frozen
  core is being asked to describe an ion whose valence electron has LEFT — exactly the regime where a frozen
  core is least justified, because the 1s sees a different potential once 2s is gone.  q3 cannot have that
  problem but costs functions on every Li site in a 14-atom cell.  ⇒ mint BOTH with
  `valgen --nmax 60 --floor`, and settle it on a number the voltage cares about: **the Li intercalation
  energy** (E[LiMn₂O₄] − E[λ-MnO₂] − E[Li]) computed both ways.
  ✅ **BOTH MINTED AND VALIDATED 2026-09-23.**  q1: `--q 1 --shell 0:5:0.03:2`, converged, E = −0.189367 Ha,
  0.11 mHa below the pool floor — committed to `valence_lowq_{va,sph}.bsd` as the VALENCE variant, matching
  those files' own convention (Mn is q7 not q15, Na is q1 not q9).  q3: `--q 3 --shell 0:8:0.05:60`,
  converged, E = −4.23584 Ha, gap 0.107 mHa — validated but NOT committed, because a `.bsd` block is keyed
  by ELEMENT and the two cannot coexist in one file.  ⛔ **Two pieces of plumbing the discriminator needs
  before it can run**, both found while minting: (i) a basis-file variant that can carry q3's `LI` block
  (a `BasisSetData` enum value + its two map entries — cheap, but wire it with a real consumer, not
  speculatively); (ii) `CleanupCandidates.md` row **D-SEED1** — the atomic seed library is keyed by
  (Z, functional) and returns the FIRST match, so "neutral Li" is ambiguous between q1's 1 electron and
  q3's 3, decided silently by file order.  The A/B is not trustworthy until that throws or is keyed on q.  If q1 and q3 agree there, q1 is free
  throughput for the CE training set; if they do not, q3 is mandatory and we have learned why.  Seed density
  too (Li⁺ is a stripped cation, so `HasAtomicSpinPair` correctly calls it non-magnetic).
- **The Mn seed ions are already checked and they are fine** (2026-09-23): Mn³⁺ (d⁴) and Mn⁴⁺ (d³) both
  generate cleanly at `--nmax 60` — converged, moments 4.000 / 3.000, ⟨r⟩ 1.078 / 1.021 bohr, shrinking with
  charge as a cation should.  ⇒ **IonicSAD is available for all three spinel compositions**, which matters
  for gate 1: the ionic seed is the basin chooser, and a neutral-superposition seed would be a much worse
  guess for Mn⁴⁺ (four electrons away in a 7-electron PP) than it was for Ni²⁺ (two).
  ⚠ **A retraction that belongs here** (it produced a work item that no longer exists): "Ni²⁺ cannot be
  seeded, high-spin d⁸ is minority-d³ in a five-fold shell and the atom occupies whole irreps" was WRONG.
  The evidence — a non-aufbau run with ⟨r⟩ = 2.91 bohr for a cation — predated the iteration-cap fix made in
  the same session; at `--nmax 120` Ni²⁺ converges cleanly (charge 8.001, moment 2.000, ⟨r⟩ 1.008, between
  Mn³⁺'s 1.078 and Ni³⁺'s 0.938).  NiO now uses IonicSAD like MnO.  What is TRUE is milder: a partially
  filled MINORITY shell under a filled majority is a long descent (two orbitals 2e-6 Ha apart straddle the
  boundary), while a partially filled MAJORITY over an empty minority is not — Mn³⁺, Mn⁴⁺, Ni³⁺ (d⁷) and
  Co³⁺ (d⁶) all converge by 60.  **Both halves of that were the cap, twice.**
- **The three spinel structures in `materials.json`.**  Primitive cells: λ-MnO₂ 12 atoms (4 Mn, 8 O),
  LiMn₂O₄ 14 (2 Li, 4 Mn, 8 O), Li₂Mn₂O₄ 16 — plus whatever magnetic decoration gate 1 settles on.
  ⚠ Lattice constants are anchors: take them from a named source and say which.

**Gate 1 — does the magnetic state survive a U change?**  (the NiO lesson; blocks everything downstream.
⚠ **Seeding is NOT the lever** — trap 3 now has three runs and two different seeds converging to the same
non-magnetic fixed point, so a better seed will not buy the order.  The levers are the mixer's
magnetisation channel (§4 row N3), kT, and U itself.)
Run λ-MnO₂ and LiMn₂O₄ at fixed U = 0, 2, 4 eV and watch the integrated site moment.  If the order dies as
it did on NiO, no loop on this material means anything and the fix (mixer preconditioning in the
magnetisation channel, §4 row N3; or a different ordering) comes first.  **Cheap: three short SCFs each.**

**Gate 2 — size the runs.**  One converged SCF per composition at Γ, timed and RSS-logged.  Estimate to
beat: ~316 basis functions for LiMn₂O₄ against MnO's 118, so expect 3–8× MnO's ~6 min ⇒ 20–50 min per SCF,
⇒ an 8-step outer loop is an overnight run per composition.  **Measure it; do not plan on the estimate.**

**Gate 3 — BUILD THE ORACLE SET, then test the screening hypothesis.**  ⚠ **This gate replaces an earlier
version that was CIRCULAR** (caught 2026-09-23): it proposed fitting one \f$\varepsilon^{-1}\f$ to move all
four measured manifolds onto "their targets", but two of those targets were the ACBN0 paper's OWN values
(MnO O-2p 2.68, NiO O-2p 3.0).  Fitting to those fits out PAO-3G incompleteness — precisely the thing we
established is NOT screening.  **A target must be independent of ACBN0.**  Against independent oracles:

| manifold | ours | independent oracle | ratio |
|---|---|---|---|
| NiO Ni-3d | 13.89 | `hp.x` **5.27** (linear response) | **2.64** ← the only SOLID point |
| MnO Mn-3d | 10.77 | literature 4–7 eV (cRPA / LR-cDFT) | 1.5–2.7 (a range) |
| NiO O-2p | 7.27 | cRPA \f$\bar U_{2p}\gtrsim4\f$ eV, "independent of the TM" (ACBN0 paper's own citation [97]) | ≲1.8 (a bound) |
| MnO O-2p | 7.36 | same | ≲1.8 (a bound) |

⛔ **Two things fall out, both unwelcome.**  (1) We have **ONE solid oracle point**, not four — a range and
two bounds otherwise — so "does a single factor fit?" is UNDERDETERMINED and cannot refute anything yet.
(2) What structure there is looks like **TWO CLUSTERS, not one factor**: d ≈ 2.6, p ≈ 1.8.

⛔ **And a physics argument runs the WRONG WAY.**  A bulk \f$\varepsilon\f$ screens everything equally.  Make
it q-resolved — \f$U\sim\sum_q|\rho_\phi(q)|^2\varepsilon^{-1}(q)v(q)\f$ — and a MORE localized orbital has
broader \f$\rho_\phi(q)\f$, samples larger q, where \f$\varepsilon^{-1}(q)\to1\f$: so TM 3d should need
**less** screening correction than O 2p.  We observe it needing **more** (2.6 vs 1.8).  That is backwards for
a dielectric picture, and it is not explained by orbital differences — our occupations match the oracle's to
1–2 %.  ⇒ **this is the specific observation that could refute (b)**, and it is available for the price of
assembling the oracle set.

★★ **AND THE POINTS ARE ALREADY BUILT: ABINIT DOES cRPA *AND* LINEAR-RESPONSE J** (user asked 2026-09-23;
verified in `~/Code/abinit/src`).  Two capabilities QE does not give us, each hitting one of gate 3's two
weaknesses:
- **`ucrpa`** (+ `ucrpa_bands`, `ucrpa_window`) — **constrained RPA**, which is *methodologically
  independent of linear response*.  This is the one that matters: the O-2p "oracle" above is a literature
  **bound** (cRPA \f$\gtrsim4\f$ eV) and cRPA is exactly the method that produced it, so `ucrpa` turns a bound
  into a VALUE on our own materials — and gives an independent Mn-3d number on MnO, where hp.x is broken by
  the d⁵ shell (trap 1).  Two of the four soft rows become solid.
- **`lruj`** (`src/98_main/lruj.F90`, "Linear Response U **and J**", citing Cococcioni & de Gironcoli PRB 71
  035105) driven by `macro_uj` / `pawuj_det` ("Determine U (or J) parameter", with/without a compensating
  charge bath).  **hp.x gives U only**, which is precisely why §5's first open question — *does J transfer
  the way U does?* — has had no oracle.  ABINIT answers it.
- ⚠ **BUT IT IS A DIFFERENT-PP ORACLE, AND THE TABLE MUST SAY SO.**  ABINIT's Hubbard is **PAW-only** (every
  Hubbard file lives under `src/65_paw`; the keyword is literally `usepawu`), and our GTH is
  norm-conserving.  The discipline that made the QE comparison worth anything was `gth2upf` — the oracle ran
  OUR pseudopotential and OUR projector, which is why the occupations matched to 1–2 % and the disagreement
  could be pinned on the functional.  We cannot do that here without a PAW dataset.  ⇒ ABINIT numbers are a
  SECOND OPINION, not a matched comparison; cRPA's methodological independence is worth more here than PP
  matching, but **label every oracle row `matched` or `different-PP`** or the next reader will average them.
- Cost note: ABINIT is an MPI build like the rest — `mpirun` always, and it additionally needs
  `--force-mpirun` (`CLAUDE.md`).  cRPA needs bands/windows chosen, which is a real input-convergence
  question of its own; treat the first run as a recipe hunt, not a number.

★ **PySCF is the OTHER install, and it is for the KERNEL, not for the values** (so: second, not first).
Free, `pyscf.pbc` gives periodic Gaussians with k-points and density fitting, and — the part we would
actually use — a four-index ERI engine with **range-separated (`omega`) integrals**, i.e. the exact object
route (b)'s screened kernel needs checking against.  It buys nothing for gate 3.
⚠ **Not `pip install` on this box as it stands** (checked 2026-09-23): Python is **3.14.4** and
`/usr/lib/python3.14/EXTERNALLY-MANAGED` is present, so PEP 668 makes apt-managed Python refuse installs
into system site-packages.  `python3-pip` provides BOTH `pip3` and `pip` (they are the same thing here —
there is no python2), but the name is not the issue; the venv is:
```
sudo apt install python3-venv python3-pip
python3 -m venv ~/Code/pyscf-env && ~/Code/pyscf-env/bin/pip install pyscf
```
Inside a venv `pip` and `pip3` are identical, so the question stops mattering.  ⚠ **Python 3.14 is new
enough that a PySCF manylinux wheel may not exist yet** — if pip falls back to building from source it will
want cmake + libcint, which is where the project's standing "source builds are the DEFAULT" policy takes
over anyway (and is arguably what we want, since the point of PySCF here is to read its integral engine).
`python3-numpy` / `python3-scipy` are in apt if a source build needs them outside the venv.

★★ **FIRST RESULT, 2026-09-23 — THE ONE-FACTOR TEST PASSES ITS FIRST REAL TRIAL, AND IT CORRECTS THE
ARGUMENT ABOVE.**  Run on paper with the newly installed PySCF (2.14.0, `~/Code/pyscf-env`), using its
range-separated (`omega`) four-index integrals on OUR OWN contractions read off the `[+U radial]` banner —
`scripts/gate3_screening_test.py`:

| manifold | bare \f$F^0\f$ | independent target | needed ratio | needed \f$\omega\f$ |
|---|---|---|---|---|
| Ni 3d | 24.29 eV | 5.27 (hp.x, **matched PP**) | 0.217 | **1.039 a.u.** |
| O 2p | 20.97 eV | ≳4 (cRPA **bound**) | 0.191 | **0.975 a.u.** |

**The two manifolds want the same screening length to 6 %.**  ⇒ ⛔ **MY "TWO CLUSTERS, AND THE LOCALIZATION
ARGUMENT RUNS THE WRONG WAY" OBJECTION WAS APPLIED TO THE WRONG QUANTITY** and is withdrawn.  It compared
ratios of our **\f$U_{\rm eff}\f$**, which already carries ACBN0's \f$\bar N^2\f$ renormalisation — i.e. the
projector-completeness effect — so it was measuring completeness and screening mixed together.  Against the
**BARE \f$F^0\f$**, which is what a screened kernel actually modifies, the picture is one factor.  And the
localization physics comes out RIGHT, not backwards: the same \f$\omega\f$ screens the compact Ni 3d less
(ratio 0.217) than the diffuse O 2p (0.191), which is the required direction.
⚠ **Three things keep this from being a verdict.**  (1) \f$1/\omega\approx1\f$ bohr is a SHORT length —
inside the 3d orbital itself and well inside the 3.94 bohr Ni–O bond — so whatever this is, calling it
"screening by the medium" needs an argument; an `erfc` cut-off at 1 bohr reshapes the on-site
self-interaction rather than dressing it.  (2) TWO points, one of which is a BOUND, is not a fit.  (3) PySCF's
bare \f$F^0\f$ = 24.29 eV against our own run's reported bare \f$\bar U\f$ = 27.3 eV — **11 % apart, and
unexplained**; ours is a density-matrix-weighted eq-10 average rather than the plain
\f$(2l+1)^{-2}\sum_{mm'}(mm|m'm')\f$, which probably accounts for it, but two codes 11 % apart on "the same"
number is exactly what a cross-check exists to catch.  **Resolve (3) before quoting any of this.**
⇒ this RAISES the value of the ABINIT runs below: turning the O 2p bound into a cRPA VALUE, and adding
Mn 3d / Mn O-2p rows, is what turns a suggestive two-point coincidence into a test.

**So the gate is: get more independent points BEFORE building a kernel.**  ★ The reframing that makes this
affordable: **(c)/(d) are unaffordable PER COMPOSITION but perfectly affordable ONCE.**  Use them for what
they are good at — a calibration set on a handful of materials where linear response is healthy (`hp.x`
run properly, i.e. with its \f$U_{\rm in}\f$ declared) — then test the screened kernel against it, then run
(b) in production across the composition sweep.  No amount of cleverness about (b) fixes a one-point oracle
set.  **Refutation criterion:** if d and p need systematically different factors AND the localization sign
stays wrong, the residual is not bulk screening and (b) is refuted; (d) becomes the route and this plan
changes shape.

**Gate 4 — SUPERCELL vs q-MESH IN OUR CODE, and the k/irrep parallel axis** (user, 2026-09-23:
*"they may not be equivalent on our code"*).  §2's equivalence is MATHEMATICAL.  Computationally the two
diverge, and on this box they diverge in the same direction for three independent reasons:
- **Asymptotics favour the q-mesh by \f$N_q^2\f$.**  A supercell is one \f$N_qN\f$ problem — \f$O((N_qN)^3)\f$ in
  the dense-algebra part.  A q-mesh is \f$N_q\f$ problems of size \f$N\f$ — \f$O(N_qN^3)\f$.  At \f$N_q=8\f$ that
  is 64× in the cubic term.
- ⛔ **RAM is our binding constraint, and the supercell route is RAM-hostile.**  Grids scale with cell
  volume, matrices as \f$(N_qN)^2\f$; the tracker's own estimate says **the 32-atom MnO supercell is
  UNTESTED**.  A q-mesh runs one small cell at a time and fits trivially.  On 14 GB this may be the
  difference between slow and *impossible* — which is a fact to MEASURE, not to reason about.
- ⛔ **A supercell perturbation breaks the symmetry we lean on.**  Our imposed runs fold the Becke mesh
  11.55× on MnO AFM-II (12 ops); an isolated \f$\alpha\hat P_I\f$ destroys that by construction.  The q-mesh
  route keeps each primitive-cell run's symmetry.
- ★ **And the q-mesh is embarrassingly parallel over q** — which is the k/irrep parallel axis we already
  own a row for.

★★ **k-PARALLELISM IS ON THE CRITICAL PATH FOR EVERY ROUTE, (b) INCLUDED.**  Not just for (c)/(d): the
ACBN0 estimator consumes converged orbitals, we measured **~10 % k-sensitivity** in MnO's U between Γ and
2×2×2, and a trustworthy production U therefore needs a k-mesh SCF per composition.  §2 row **KP** is parked
"*BY AGREEMENT not now … it should stay on the list until we are suitably embarrassed*" — the parking
condition was "no multi-k rows exist to measure the payoff on".  **This plan is the embarrassment.**  What
KP needs is already scoped there: the pre-warm exists, and one k-dependent write remains inside the loop
(`tDynamic_HT_Imp::GetMatrix`'s `mutable CacheMap itsCache` keyed by `Irrep`, plus `itsByL/itsByLSeen`).
⚠ **And there is a PREREQUISITE DEFECT**: the shifted-MP fold lowers Si by 1.02 mHa (§4, found 2026-09-20)
— "a 1 mHa fold error on the fractional-k mesh is a sampling defect that every future k-mesh run inherits;
fix before KP-1".  Multi-k U values inherit it too.

**The measurement this gate asks for** (cheap, and it settles an assumption the tracker has carried
untested): run **one 32-atom MnO 2×2×2 supercell SCF at Γ**, logging wall and peak RSS.  That single number
says whether route (d)-by-supercell is available to us at all, or only in principle.

**Then, and only then: the screened kernel.**  A new integral type on `BareCoulombSource` (the pseudo-wall
pin allows exactly this), with \f$\varepsilon\f$ computed — Thomas–Fermi first as the cheap bound, reading
`qs_linres_polar_utils.F` before committing to anything more.

**The Oct 6–20 long runs** (only what the gates have justified): the three compositions × the outer loop,
plus a NiO re-run as the control once the magnetic-robustness fix exists, plus the k-mesh arms gate 4 says
we need (Γ-only was ~10 % on MnO's U) — and, if gate 4's 32-atom measurement says the supercell route fits,
the calibration LR run that turns gate 3's one solid oracle point into several.  ⛔ `scripts/memsafe -p`,
never bare: these are unattended.

---

## 5a. UNFINISHED — started and not completed (distinct from §5, which is things we do not KNOW)

⚠ **These will rot silently if nobody looks.**  Each says what exists, what is missing, and where the
evidence is.  A fresh session should clear or re-park them before starting new work.

1. **The 11 % bare-\f$F^0\f$ discrepancy — resolve FIRST, it gates every screening number.**  PySCF reports
   24.29 eV for the Ni 3d shell-averaged bare \f$F^0\f$ on our own contraction; our own run's banner reports
   bare \f$\bar U\f$ = 27.3 eV.  The likely cause is that ours is the density-matrix-weighted eq-10 average
   and PySCF's is the plain \f$(2l+1)^{-2}\sum_{mm'}(mm|m'm')\f$ — **likely is not verified**, and two codes
   11 % apart on nominally the same quantity is what a cross-check exists to catch.  Gate 3's ω values are
   computed from the PySCF number, so they inherit it.  Evidence: gate 3's first-result block above;
   `scripts/gate3_screening_test.py`.
2. ⛔ **NOTHING CHECKS THE BASIS AGAINST THE PSEUDOPOTENTIAL** (added 2026-09-23 — this one had been said
   aloud and never written down, which is how it nearly got lost).  A run declares its PP variant through
   `SolidCalcOptions::species` (`{"Li",3}`) while the basis comes from a `.bsd` file whose blocks are keyed
   by ELEMENT only.  So a q3 run can be handed the q1 `LI` block — a basis built and validated for a
   one-electron valence, describing three — **silently**, with no diagnostic anywhere.  It is latent today
   only because every element in `valence_lowq_*` happens to be the variant its runs use.  The shape of the
   fix: the `.bsd` files already carry rich headers, so a machine-readable per-element provenance line
   (invisible to a Gaussian94 reader) lets the factory assert each block's q against the run's declared
   valence, and THROW on a mismatch.  Same family as `CleanupCandidates.md` **D-SEED1** (the seed library
   keyed on (Z, functional)), and it becomes live the moment `valence_semicore.bsd` exists.
3. **Li q3 is validated but uncommitted.**  `--q 3 --shell 0:8:0.05:60`, converged, E = −4.23584 Ha, gap
   0.107 mHa.  It cannot share a file with q1, so it needs `valence_semicore.bsd` + a `BasisSetData` enum
   value and its two map entries.  Left unwired deliberately (no consumer yet, and a dead file is worse than
   a recorded command) — the command is in the `.bsd` header and in §4's prerequisites.  ⇒ item 2 above
   should land WITH it, not after.

---

## 5. Open questions (write the answer here when it is earned)

- Does **J** transfer the way U does?  ACBN0 gives J for free; hp.x does not give us a J to check it
  against.  The atomic limit (\f$J\approx1\f$ eV for 3d) is the only oracle we have — and our bare
  \f$\bar J\approx7.5\f$ eV is NOT Hund's J (it carries eq 13's self-terms), which is why only
  \f$U_{\rm eff}=\bar U-\bar J\f$ is quotable.
- Is the per-site-oxidation-state assignment stable when two Mn sites are crystallographically equivalent
  but electronically inequivalent (charge ordering)?  That is a symmetry-breaking question and the
  imposed-symmetry machinery has an opinion — check it does not average the two.
- If a computed \f$\varepsilon\f$ needs DFPT, we will have BUILT the machinery that makes a DFPT-based (d)
  nearly free — hp.x's own method with \f$\alpha\hat P\f$ in place of the electric field, CP2K's
  `qs_linres_*` as the Gaussian reference.  ⇒ **(b) and (d) share most of their cost**, which makes the
  choice between them a hedge rather than a bet.  `doc/OpenWork.md` step 5 item 5's "never build DFPT for
  this" was costed against QE's plane-wave implementation and assumed DFPT would be built ONLY for U; if (b)
  forces it anyway, reopen that ruling rather than inherit it.
- The spherical atom resolves a partially-filled degenerate shell by picking orbitals, not by occupying the
  shell uniformly — energetically converged but symmetry-broken.  Harmless for a seed (spherically averaged
  anyway); worth a thought for the **+U atomic radial**, which takes "the lowest occupied l orbital" and on a
  broken shell that need not be the spherical average.  Measure before caring: the captured-norm line already
  prints (NiO 0.9999, MnO VA 0.991).
- GGA before any value comparison with the PBE literature (the paper's 7.63/3.0 and the 4–7 eV range are
  both PBE).  Still open, still gating the *value* comparisons, not the *method* work.
