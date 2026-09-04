# The lattice screener — a DIP seam for the collocation tolerance

**STATUS: SPECCED, NOT BUILT (2026-08-28).  Queued for its own session.**  This file exists because the
decision was taken and the reasoning is worth not re-deriving; the code is one experimental `bool` today
(`GPW_DAWARE_SCREEN`), which is explicitly NOT the design.

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

## 3. The interface

```cpp
//! The tolerance policy for one (shell pair, offset) collocation term.  GEOMETRY IN, TOLERANCE OUT.
class LatticeScreener
{
public:
    virtual ~LatticeScreener() = default;
    //! The box tolerance for component pair (i,j) at cell offset n.  >= the global floor, always.
    virtual double Tolerance(size_t i, size_t j, const ivec3_t& n) const = 0;
    //! May this (pair, offset) be dropped whole?  \a pf is PairPrefactorExp (screen 1).
    virtual bool   Kill(size_t i, size_t j, const ivec3_t& n, double pf) const = 0;
    //! Is the answer INDEPENDENT of the density?  True ⇒ every BoxGeom, chord and bound is
    //! iteration-invariant and may be hoisted into the task list (see §5 — this is the point).
    virtual bool   IsGeometryOnly() const = 0;
};
```

Two concretes:
- **`GeometryOnlyScreener`** — returns \f$\varepsilon\f$ and never kills.  CP2K's rule.  Stateless.
- **`DAwareScreener`** — today's rule.

## 4. ⚠ OPEN DESIGN QUESTION — HOW DOES `DAwareScreener` REACH THE DENSITY?  **SETTLE THIS FIRST**

**Do not start coding until this is decided** (user, 2026-08-28, having second thoughts about the proxy
sketched in the first cut of this file).  It is the only genuinely contentious part of the design, and it
is the part that decides whether the seam removes a bug class or relocates one.

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

## 5. ★ THE PRIZE, and it is not just tidiness

`IsGeometryOnly()` is the reason this is worth doing now rather than as hygiene.  \f$\varepsilon_{eff}\f$
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

| per call | D-aware (today) | geometry-only | cost |
|---|---|---|---|
| Si Γ, collocate | 0.003296 s | 0.003567 s | **+8%** |
| **MnO parity probe, collocate** | **1.839 s** | **2.100 s** | **+14.2%** |
| **MnO parity probe, gather** | **1.876 s** | **1.982 s** | **+5.7%** |

Trajectories move (Si Γ: 11 iters CONVERGED → 17 OSCILLATING, `Etot` agreeing to **2e-9**; the MnO probe's
call counts go 68/132 → 82/139), so read the PER-CALL column — whole-run CPU (386 → 464 s) mixes the cost
with the extra calls.

⇒ **THE D-AWARE SCREEN BUYS ~10% ON THE BOX WALK, WEIGHTED.**  That is the honest price of dropping it, and
it is the first time it has been paid knowingly.  Set against it: the per-line term this would make
hoistable is **~40% of the kernel** (§5), the defect list in §2, and CP2K parity.  ⇒ The trade looks
clearly favourable — **pay ~10%, plausibly save ~30%** — but the ~30% is an ESTIMATE and the §5 RAM
question is unresolved, which is exactly what the queued session is for.  **Do not judge the idea on the
knob**: the knob only removes the screening, it hoists nothing.

## 7. Where the choice is made

> *"The decision on which screener to construct and pass into the code is made at a very high level like in
> `SolidCalculation`."*

Beside the imposition veto and the XC-mesh policy, which are already resolved there from
`qchem::RunPolicy`.  ⇒ It becomes a declared `CP2K_COMPAT` deviation like the others — and it is one we
have been taking silently, since CP2K screens on geometry alone.

## 8. Not now, but the framework should not preclude it

> *"Once this framework exists we will then have the flexibility to easily add other screening algos
> (maybe a special cubic version for speed?) if needed."*

A `CubicScreener` (an axis-aligned box rather than a reach sphere) would drop the chord `sqrt` entirely at
the price of walking ~1/(π/6) ≈ 1.9× more points.  It is a bad trade TODAY — the chord/interval skip was
worth 1.40× on MnO when it landed — but on a machine with different arithmetic/memory balance, or for a
GPU port where branch-free beats fewer points, it may not be.  The seam is what makes that an experiment
rather than a rewrite.
