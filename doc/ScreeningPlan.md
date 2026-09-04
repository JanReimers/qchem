# The lattice screener — a DIP seam for the collocation tolerance

**STATUS: THE SEAM IS BUILT (2026-09-04).  §5's prize is NOT -- that is the next increment.**
`src/BasisSet/Molecule/LatticeScreener.C` carries the interface and both concretes; the two collocation
faces take a `const LatticeScreener&`; `GPW_DAWARE_SCREEN` survives only as
`qchem::RunPolicy::DAwareScreen`, which now selects a screener OBJECT instead of a branch in the box walk.
D-aware remains the default and the whole suite is unchanged.  §4 was settled first, as this file
demanded, and the answer is in §4 below -- it is not any of A/B/C.

## 1. What is being decided

Which tolerance sizes a collocation box, and kills a (pair, offset) term:

| today | \f$\varepsilon_{ij}=\varepsilon/|c_{ij}|\f$, with \f$c_{ij}=\mathrm{fold}\cdot\mathrm{Re}[D_{ij}e^{-ikR}]\f$ |
|---|---|
| the alternative | \f$\varepsilon\f$ flat — geometry only, which is what CP2K does |

## 2. Why it needs a seam rather than a flag (user, 2026-08-28)

> *"I much prefer optionality to work through virtual dispatch rather than if statements in the low level
> code.  This makes the big code base much easier to work with in the long run."*

⇒ **DEPENDENCY INVERSION (SOLID::DIP).**  The box walk asks an abstract `LatticeScreener` for a tolerance;
it does not know, and must not branch on, which rule answers.  A `bool` consulted inside the innermost
geometry code is the thing being removed, not the thing being added — and it is worth saying WHY beyond
taste, because this particular flag has a record:

⛔ **THE D-AWARE TOLERANCE IS THE MOST DEFECT-DENSE IDEA IN THE COLLOCATION PATH.**  Its bug list, all
found the hard way and all recorded elsewhere in these docs:
- the **shifted-MP defect** — the integrate-back screened on \f$\mathrm{Re}[D\overline{e^{ikR}}]\f$ where
  it needed \f$|D|\f$, so at quarter-integer \f$k\f$ every odd-offset term was discarded and \f$H\f$ came
  out real.  **4.1 Ha**, and it survived in a SECOND copy until 2026-08-27 because the value cache made
  that copy unreachable;
- **`FoldScreenMax`** — an orbit reduction that used the signed average instead of the max cost a **2.3 Ha**
  variational collapse on the imposed O₂ triplet;
- the **union-vs-per-component box**, an entire design decision (one box must serve a shell pair, so it is
  sized by the tightest \f$\varepsilon\f$ present and each component re-screens) that exists ONLY because
  the tolerance is per-pair;
- the **reduced-vs-full tier split** in the T3 fold gates, re-tiered 2026-08-28: at general \f$k\f$ the
  offset phase varies across an orbit, so the two arms pick different union tolerances and agree only at
  the screening tier.

Every one of those is a consequence of the tolerance depending on \f$D\f$.  A seam does not fix them, but
it makes the alternative reachable, testable and switchable — and it stops the next screening idea being a
fifth branch in the same function.

## 3. The interface — AS BUILT

Two virtuals, one dispatch per **task** (a whole shell-pair block at one offset), no raw pointers, no
per-task allocation.  `src/BasisSet/Molecule/LatticeScreener.C`:

```cpp
//! The screen's verdict for ONE (shell pair, offset) task.  Caller-owned and REUSED across the offsets
//! of a shell pair, so screening a task allocates nothing after the first.
class ScreenPlan
{
public:
    double Tolerance()                const;   //!< the shared box's eps; 0.0 = the task is dead
    bool   Dead()                     const;
    bool   Keeps(size_t a, size_t b)  const;   //!< does component (a,b) of the block ride the box?
    ...                                        //!< Reset/Keep/SetTolerance: the screener's side
private:
    std::vector<char> itsKeep;                 //!< NOT vector<bool> — see below
    size_t itsNJ; double itsEps;
};

//! The tolerance rule for a (shell pair, offset) collocation task.  STATELESS and RUN-LIFETIME.
class LatticeScreener
{
public:
    virtual ~LatticeScreener() = default;
    //! \param cij  the nI x nJ block of screening WEIGHTS (0 = the term is absent).
    //! \param plan OUT: the union tolerance + the survivor mask.
    virtual void   Screen(const ivec3_t& n, double pf, const rmat_t& cij, ScreenPlan& plan) const = 0;
    //! No Screen may answer below this, so the static task list stays a strict superset.
    virtual double Floor() const = 0;
};
```

Two concretes: **`GeometryOnlyScreener`** (flat \f$\varepsilon\f$, weights ignored — CP2K's rule) and
**`DAwareScreener`** (\f$\varepsilon/|c_{ij}|\f$ — today's, and the default).

### Four things the first cut of this section had wrong

- ⛔ **`GeometryOnlyScreener` does NOT "never kill".**  The \f$pf\ge-\ln\varepsilon\f$ prefactor test IS
  the geometry screen — CP2K's fixed `radius_list`.  Shipping "never kills" would have silently walked
  every far-offset task.  Pinned by `M_LatticeScreener.GeometryOnlyStillKillsOnThePrefactor`.
- ★ **BATCHING over a whole task, not per (i,j).**  One virtual call per task instead of \f$2n_In_J\f$ —
  and, more importantly, it moves the **union rule inside the screener**, where it belongs.  That rule was
  open-coded *identically at both call sites* before the seam, which is precisely the duplicated-decision
  shape §2's bug list is made of.
- ★ **`IsGeometryOnly()` WAS DROPPED — it had no consumer.**  Its only candidate is `BoxTasks()`, i.e. §5,
  which is not built yet.  Rather than ship a virtual whose contract we would be guessing at, note that the
  invariance is a *testable property of the object* (`M_LatticeScreener.GeometryOnlyIgnoresTheWeights`)
  rather than a claim it makes about itself.  When §5 arrives the natural shape is not a predicate at all
  but the screener owning the memo (`BoxFor(task)` — cached under geometry-only, fresh under D-aware).
- **The weight is NOT \f$D_{ij}\f$**, so the parameter is not named that.  It is
  \f$c_{ij}=\mathrm{fold}\cdot\mathrm{Re}[D_{ij}e^{-ikR}]\f$ collocating (fold doubling × orbit
  projection × phase contraction) and \f$\mathrm{fold}\cdot|D_{ij}|\f$ gathering — see §4's last note and
  the asymmetry warning below.

### Why `std::vector<char>` and not `std::vector<bool>`

Raised as a style question (*"bool is the right modern expression"*) and worth recording, because the
intuition inverts.  `std::vector<bool>` is a bit-packed **partial specialization that is not a container of
bool**: `keep[k]` yields a proxy rather than a `bool&` (so `auto` silently captures the proxy), it has no
`data()`, each access costs a shift and a mask where `char` costs a byte load, and neighbouring elements
share a word so they **cannot be written elementwise from parallel workers** — and this walk runs under
OpenMP.  Alignment gets *worse*, not better: elements are not individually addressable at all.
⇒ `char` is the storage; `ScreenPlan::Keeps()` is the boolean *expression*, so no call site names the type.

### ⚠ The two directions still screen on different quantities

Deliberately kept (user, 2026-09-04: *"we can keep the asymmetry for now"*).  Since
\f$|D|\ge|\mathrm{Re}[De^{-ikR}]|\f$, the gather's tolerance is the tighter one and its kept set is a
strict **superset** of the scatter's — so the "identical active set / exact adjoint" claim that the
`LatticeSum1E` face and the walk comments used to make flatly is exact only up to that asymmetry at general
\f$k\f$.  Both faces now say so.  The magnitude form on the gather side is the *fix* for the 4.1 Ha
shifted-MP defect, so unifying the two means changing the scatter, which moves numbers — a separate, gated
change, and explicitly NOT what this seam did.

## 4. ✅ SETTLED (2026-09-04): **IT DOESN'T.**  THE SCREENER HOLDS NO DENSITY AT ALL

The question below was posed as a choice between three ways of giving a `DAwareScreener` access to \f$D\f$.
The answer turned out to be that it needs none, and the two findings that got there are worth keeping
because both were invisible from the spec and only appeared on reading the call sites:

1. ⛔ **OPTION B AS WRITTEN IS NOT IMPLEMENTABLE.**  A D-aware screener must bind to the *orbit-projected*
   density — `FoldProjectedD` collocating, `FoldScreenMax` gathering — and both are **private members of
   the MnD evaluator**, reading `itsStreamFold`.  The SCF caller cannot construct such a screener; it has
   no access to the projection.  So B needs the two-level `Policy::Bind` split after all — forced by
   ENCAPSULATION, not by the run-lifetime state §4 anticipated.
2. ★ **THE WALK ALREADY COMPUTES THE WEIGHT.**  `scatterShell` computes \f$c_{ij}\f$ because it is the
   scatter weight, needed whatever the screening rule; `integrateShell` computes \f$|D_{ij}|\f$ for its
   screen.  Hand that scalar *in* and the screener needs no density, no binding, and no lifetime.

⇒ **A screener is a STATELESS, RUN-LIFETIME rule**: chosen once, shared by every k-block, every SCF
iteration and every OpenMP worker.  `Tolerance` is a pure function of its arguments; there is nothing to
reseat and nothing that can go stale.  The staleness class §4's Option A table lists as this tree's most
expensive recent defect shape is **unrepresentable**, rather than guarded against by a runtime check.
The A-vs-B-vs-C argument below is kept as the record of how the question was framed.

⚠ **ONE CLAIM IN THE B COLUMN DID NOT SURVIVE**, and it is worth naming because it was the argument that
"decided it": *"`IntegratePotential`'s `screenD` parameter **becomes** the screener: same arity."*  It does
not.  A stateless screener needs the caller to supply weights, and the gather direction has no density
parameter other than `screenD` — so `screenD` stays as the WEIGHT SOURCE and the screener is an additional
parameter carrying the RULE.  `screenD == nullptr` now means "no weight information", which the walk
presents as a unit weight and every rule answers at its floor — exactly what the static sharp-field sweeps
(local PP, `V_loc`-long) have always got.  That part of the simplification was real; the arity claim was not.

### The question as it was posed

**The constraint.**  `Tolerance()` needs \f$c_{ij}=\mathrm{fold}\cdot\mathrm{Re}[D_{ij}e^{-ikR_n}]\f$.
\f$D\f$ changes every SCF iteration.  The screener is chosen once, at the top, and must be reachable from
the innermost box walk — which runs under OpenMP over shell pairs.

### Option A — the reseatable proxy (what the first cut proposed)

`DAwareScreener` holds a `shared_ptr<const DensityHandle>` at construction; the SCF calls
`handle->Reseat(&D)` once per iteration.  The screener object outlives the densities it reads.

| for | against |
|---|---|
| screener built once; nothing per-iteration to allocate | ⛔ **`Tolerance(i,j,n)` stops being a pure function** — same arguments, different answers depending on *when* you ask |
| the density never appears in the interface | ⛔ a **raw reseatable `DM*`**, i.e. a global mutable with a nicer name.  Who guarantees it outlives the read? |
| | ⛔ **thread safety becomes an unwritten invariant** (reseat must happen outside every parallel region) |
| | ⛔ **staleness is now possible and only a runtime check can catch it** — and this tree already carries three such checks for exactly this shape: the `CollocMemo` depth-1 thrash, the XC `itsSrcVersion` "STALE" warning, and the cross-run `Projector3` pollution that had to be fixed by instance-scoping |

### Option B — construct the screener per iteration and pass it by `const&`  ★ RECOMMENDED

The SCF (or the closure that already holds \f$D\f$) makes a fresh screener bound to this iteration's
density and passes `const LatticeScreener&` down.  Immutable once constructed.

★ **AND THE DENSITY IS ALREADY A PER-CALL PARAMETER ON BOTH FACES, WHICH IS THE ARGUMENT THAT DECIDES IT
FOR ME.**  This is not a new dependency being threaded through — it is a dependency that is already there,
untyped:

```cpp
// today
std::vector<rvec_t> CollocateDensity(const chmat_t& D, ...) const;
chmat_t IntegratePotential(..., const chmat_t* screenD, ...) const;   // <- the density, as a bare matrix

// option B
std::vector<rvec_t> CollocateDensity(const chmat_t& D, const LatticeScreener&, ...) const;
chmat_t IntegratePotential(..., const LatticeScreener&, ...) const;   // <- SAME ARITY, better type
```

`IntegratePotential`'s `screenD` parameter **becomes** the screener: same information, with the policy
attached instead of implied.  That is a strict simplification of a signature, not an addition to one.

| for | against |
|---|---|
| ✅ `Tolerance()` is a **pure function of its construction arguments** — the staleness class cannot exist | one object constructed per iteration (it holds a reference and a phase, not a copy — nothing to fear) |
| ✅ trivially thread-safe: immutable, shared by every worker | the two faces gain a parameter (see above: one of them LOSES one) |
| ✅ you cannot ask a screener about the wrong \f$D\f$, because there is no other \f$D\f$ it could answer for | |

### Option C — hand down the tolerances as data, not the rule

The SCF computes \f$\varepsilon_{ij}(n)\f$ into an array once per iteration; the walk just reads it.
⛔ Materialises a per-iteration array over the task list (~35 k entries on MnO) and computes tolerances for
tasks that are then killed anyway.  It also puts the RULE back at the call site, which is what the seam
exists to remove.  Mentioned for completeness; I do not recommend it.

### ▶ MY RECOMMENDATION — B, with a two-level split if state ever justifies it

Take **Option B**.  The decisive point is that the density is already flowing through both faces per call,
so B types an existing dependency rather than introducing one, while A introduces mutable shared state to
hide a dependency that was never hidden in the first place.  This codebase's own bug record is the
tiebreaker: its most expensive recent defects have all been *"an object answered for the wrong iteration"*,
and A is that shape by construction while B makes it unrepresentable.

⚠ **The one real argument for A is if the screener needs run-lifetime STATE** — a precomputed \f$|D|\f$
magnitude matrix, `FoldScreenMax`'s orbit maxima, or anything else too dear to rebuild per iteration.  If
that appears, do NOT reach for the proxy; split instead:

```cpp
class LatticeScreenerPolicy      // RUN lifetime, immutable, chosen in SolidCalculation, captured by closures
{ public: virtual ScreenerPtr Bind(const chmat_t& D, const cellphase_t&) const = 0; ... };

class LatticeScreener            // ITERATION lifetime, immutable, passed by const& into the walk
{ public: virtual double Tolerance(...) const = 0; ... };
```

That keeps *"the decision is made at a very high level"* (the POLICY travels and is captured), keeps the
hot object immutable, and still has nowhere for a stale pointer to live.  ⚠ Do not build the two-level
form up front — one level is enough for two implementations, and the split is a refactor away if a
screener ever earns it.

⚠ **A NOTE ON THE FRAMEWORK-CACHED CLOSURES**, since it is the first objection someone will raise:
`MakeCollocator` / `MakeIntegrator` build closures ONCE and call them per iteration, which sounds like it
forces A.  It does not — the closures already receive \f$D\f$ as a call argument and merely need to capture
the *policy* (immutable) and bind per call.  Checked before recommending B.

## 5. ⛔ THE PRIZE IS NOT THERE — MEASURED 2026-09-04, AND THIS INCREMENT IS CLOSED

**Do not build the geometry hoist.**  Priced at the cube by
`M_PG_BoxWalk.WhatTheGeometryHoistWouldBuy` (seconds, no SCF), on the MnO Mn–O contact at the raster the
`CP2K_COMPAT=1` row actually runs (N=40³):

| | |
|---|---|
| `ContractCube`, whole | 133.1 µs |
| the chord arithmetic ALONE | 18.5 µs — **13.9%** |
| the same loop reading precomputed bounds | 1.0 µs |
| ⇒ **HOIST CEILING** | **17.5 µs = 13.2% of the kernel, and not more** |
| RAM it would cost | 21.9 kB/task × 6152 tasks = **132 MB** |

And one box is not an argument, so the share was swept against box size — it goes like \f$1/n\f$ (the chord
is per *line*, \f$O(n^2)\f$; the point work inside it is \f$O(n^3)\f$), meaning a SMALL box is the hoist's
best case:

| box \f$n\f$ | 99 | 71 | 53 |
|---|---|---|---|
| chord share | 9.9% | 13.5% | **17.3%** |

⇒ **10–17%, not the ~40% this section claimed.**  Set that against the measured price of the screener that
makes the hoist legal at all — **+14.5% wall** (§6) — and the best case is a WASH, before paying 132 MB
against a run whose whole peak RSS is 110 MB.  The estimate was wrong in the one direction that mattered.

⚠ **WHERE THE ~40% CAME FROM, so it is not re-derived**: the per-LINE work really is a large share, but most
of it is the \f$e_2\f$ fold, and the fold reads the polynomial coefficients — which carry the density
weights.  It is not hoistable under ANY screener.  Only the chord *bounds* are iteration-invariant, and they
are the small half.

★ **THIS IS THE exp-RECURRENCE PATTERN, CAUGHT EARLY THIS TIME** (doc/OpenWork.md): that idea removed 20% of
the profile on paper and returned 3% after a fully gated implementation.  Here a cube-level number cost one
test and closed the increment before a line of it was written.  ⇒ The seam (§1–§4) still stands on its own
merits — the defect record in §2, CP2K parity, and the fact that the tolerance rule is now switchable — but
it should NOT be sold as an enabler for this.

### The original argument, kept as the record of what was believed

**⛔ `IsGeometryOnly()` WAS NOT SHIPPED AND MUST NOT COME BACK** (user, 2026-09-04).  It was left out
because it had no consumer yet; that is the weaker reason.  The real one:

> *"It is like adding `IsSlaterBasisSet()`, `IsGaussianBasisSet()`, `IsPolarizedGaussianBasisSet()`,
> `IsPWBasisSet()` to `IrrepBasisSet`.  It goes counter to most of the SOLID principles."*

Exactly so — and the timing objection would have evaporated the moment the consumer appeared, while this
one does not.  The predicate asks the object **which subclass it is** so the caller can branch; the answer
is constant per type, so Open/Closed dies at once (a third screener = a third branch at every asking site).

**THE DISTINCTION TO HOLD ON TO** — it is not about returning a `bool`:

| ❌ an IDENTITY question | ✅ a question about THIS REQUEST |
|---|---|
| *"Are you the geometry-only kind?"* | *"Can you size **this** box without a density?"* |
| one answer per class, forever | answered per task, on the merits |
| caller branches on TYPE | caller branches on an ANSWER |
| a new screener breaks callers | a new screener just answers |

The second form also degrades gracefully: a rule that is density-independent for *some* pairs and not
others can answer truthfully per task, where the predicate would have to lie.  It is the same difference as
`IsGaussianBasisSet()` versus a cross-cast that returns null — the idiom this tree already uses.

### ⚠ AND THE MEMO MUST NOT LIVE IN THE SCREENER EITHER

An earlier sketch in this file proposed the screener owning the geometry cache (`BoxFor(task)` — memoized
under geometry-only, fresh under D-aware).  **Retracted.**  A caching screener is not STATELESS, and
statelessness is the property that made §4 collapse: no reseat, no staleness, safe to share across every
OpenMP worker.  Handing that back — to the very object whose defect record §2 catalogues — to save a
recompute would be a bad trade.  It would also drag `BoxGeom`, an evaluator internal, into the screening
module.

⇒ **The cache belongs in `BoxTasks()`**, where the geometry already lives.  The screener contributes only
whether it can answer yet:

```cpp
//! Size this task from GEOMETRY ALONE.  false = this rule needs the per-term weights, so the geometry
//! cannot be hoisted and the caller must screen per iteration.  Asked ONCE PER TASK at task-list build
//! time, never on the hot path.
virtual bool ScreenStatic(const ivec3_t& n, double pf, size_t nI, size_t nJ, ScreenPlan&) const = 0;
```

`GeometryOnlyScreener` fills the plan and returns true; `DAwareScreener` returns false without inspecting
anything.  The screener stays immutable, no call site learns its type, and the RAM question below is then
the only real decision left — which is where this increment's effort should go.

The reasoning the predicate stood for is unchanged:  \f$\varepsilon_{eff}\f$
is the ONLY \f$D\f$-dependent input to `BoxGeom` — the metric, the step vectors, the product centre, the
half-widths and the chord quadratic are all pure geometry (user, same day, and it is right).  So with a
geometry-only screener **the entire box geometry becomes iteration-invariant**, and the per-line chord —
measured at ~40% of the contraction kernel and previously written off as "close to irreducible", because
it carries a `sqrt` — becomes SETUP that the task list can hold.

⚠ Not free, and the cost side is already measured (§6).  The RAM question is the one to design against:
storing chord bounds per task is \f$O(n^2)\f$ per task (~15 kB on MnO's boxes, ~144 MB over its 9760
tasks).  That is 27× smaller than the value cache step 7 deleted, but it is the same SHAPE of trade and
must be argued, not assumed.

## 6. Measured so far (the knob, not the design)

`GPW_DAWARE_SCREEN=0`, which only removes the screening — it does NOT yet hoist anything, so it prices the
COST side alone:

★ **RE-MEASURED 2026-09-04 THROUGH THE SEAM** (`CP2K_COMPAT=1 GPW_MNO_NMAX=10 MNO_SKIP_FM=1`, serial,
both arms in one session on one binary — the two arms produce IDENTICAL call counts (65 gathers / 22
collocations) and \f$E_{tot}\f$ agreeing to 2.3e-7, so the per-call column is a clean A/B):

| MnO, `CP2K_COMPAT=1` | D-aware | geometry-only | price |
|---|---|---|---|
| collocate, per call | **1.266 s** | **1.520 s** | **+20.1%** |
| gather, per call | **1.254 s** | **1.421 s** | **+13.3%** |
| wall (10 iters) | **1:54.6** | **2:11.2** | **+14.5%** |
| peak RSS | 111 MB | 110 MB | — |

⚠ **AND THE 08-28 ROW BELOW IS STALE — the box walk got ~1.45× faster after it was taken** (the run-length
stream geometry).  Quote the table above, not the one below, and note the price of geometry-only is HIGHER
than the old estimate: **+14.5% wall, not ~10%**.

The superseded 08-28 numbers, kept so the delta is visible:

| per call | D-aware (today) | geometry-only | cost |
|---|---|---|---|
| Si Γ, collocate | 0.003296 s | 0.003567 s | **+8%** |
| **MnO parity probe, collocate** | **1.839 s** | **2.100 s** | **+14.2%** |
| **MnO parity probe, gather** | **1.876 s** | **1.982 s** | **+5.7%** |

Trajectories move (Si Γ: 11 iters CONVERGED → 17 OSCILLATING, `Etot` agreeing to **2e-9**; the MnO probe's
call counts go 68/132 → 82/139), so read the PER-CALL column — whole-run CPU (386 → 464 s) mixes the cost
with the extra calls.

⚠ **ONE SEMANTIC DIFFERENCE vs THE OLD KNOB, on the geometry-only arm only** (2026-09-04).  The seam's
uniform convention is *0 weight = the term is absent*, so the gather now drops \f$|D_{ij}|=0\f$ terms under
BOTH rules.  The old `GPW_DAWARE_SCREEN=0` kept them, because the `Dmag==0` skip lived INSIDE the D-aware
branch — while the collocate side dropped them on both arms (`want` filters `D==0`, and `if (c==0.0)
continue` ran regardless of the knob).  ⇒ The old geometry-only arm was INCONSISTENT BETWEEN THE TWO
DIRECTIONS, which is exactly what the shared-active-set claim forbids; the seam fixes that rather than
preserving it.  The D-aware default arm is untouched.  Checked, not argued: all 29 `GPW.*` analytic gates
pass under `GPW_DAWARE_SCREEN=0`, `AnalyticIntegrateBackAdjoint` included.

⇒ **THE D-AWARE SCREEN BUYS ~10% ON THE BOX WALK, WEIGHTED.**  That is the honest price of dropping it, and
it is the first time it has been paid knowingly.  Set against it: the per-line term this would make
hoistable is **~40% of the kernel** (§5), the defect list in §2, and CP2K parity.  ⇒ The trade looks
clearly favourable — **pay ~10%, plausibly save ~30%** — but the ~30% is an ESTIMATE and the §5 RAM
question is unresolved, which is exactly what the queued session is for.  **Do not judge the idea on the
knob**: the knob only removes the screening, it hoists nothing.

## 7. Where the choice is made — AS BUILT

> *"The decision on which screener to construct and pass into the code is made at a very high level like in
> `SolidCalculation`."*

⇒ **`qchem::RunPolicy::DAwareScreen`** — the declared-deviation table, which is what `SolidCalculation`
itself reads for the imposition veto and the XC-mesh policy, and which the run banner prints.  The
`GPW_DAWARE_SCREEN` knob name is unchanged, so every measurement banked against `GPW_DAWARE_SCREEN=0`
still reproduces — but it now selects an OBJECT, not a branch.  It is CP2K-parity-`false`, so
`CP2K_COMPAT=1` picks up the geometry-only rule automatically: the deviation we had been taking silently
is now declared and reported.

⚠ **The screener is CONSTRUCTED one layer lower than the quote asks** — in `GPW_Evaluator`'s constructor
(`BasisSet/Lattice_3D/Evaluators/GPW/Imp/Evaluator.C`) rather than in `SolidCalculation`.  Reason: the
evaluator is built deep inside the IBS factory chain, which `SolidCalculation` does not thread options
through, and the same layer already reads `theRunPolicy().StreamFold()` for exactly this kind of choice
(`BasisSet/Lattice_3D/Imp/BasisSet.C`).  The *decision* is still made in one high-level named place; only
the `new` is local.  Moving it to `SolidCalculation` is a plumbing change, not a design change, if the
caller ever needs a per-run screener that `RunPolicy` cannot express.

## 8. Not now, but the framework should not preclude it

> *"Once this framework exists we will then have the flexibility to easily add other screening algos
> (maybe a special cubic version for speed?) if needed."*

⚠ **A `CubicScreener` would NOT fit this seam**, and that is worth saying now rather than discovering it
later: the box SHAPE is not the box TOLERANCE.  `LatticeScreener` answers "how tight", and the shape is
`BoxGeom`'s business — so an axis-aligned variant is a second seam (a `BoxShape` policy on `MakeBoxGeom`),
not another concrete here.  The note below stands on its own merits.

A `CubicScreener` (an axis-aligned box rather than a reach sphere) would drop the chord `sqrt` entirely at
the price of walking ~1/(π/6) ≈ 1.9× more points.  It is a bad trade TODAY — the chord/interval skip was
worth 1.40× on MnO when it landed — but on a machine with different arithmetic/memory balance, or for a
GPU port where branch-free beats fewer points, it may not be.  The seam is what makes that an experiment
rather than a rewrite.
