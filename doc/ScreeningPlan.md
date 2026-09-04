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

## 4. ⚠ HOW `DAwareScreener` REACHES THE DENSITY — through a proxy, NOT through the interface (user)

> *"the `DAwareScreener` needs to get the density matrix through a proxy that holds a `DM*` at
> construction (the pointer gets reassigned at each iteration), not through the interface."*

This is the load-bearing half of the design.  Putting \f$D\f$ in `Tolerance()`'s signature would push the
density into every screener that does not want it — an ISP violation, and it would make
`GeometryOnlyScreener` take an argument it must ignore.  Instead:

```cpp
//! The density the screener reads, indirected so the SCF can re-seat it per iteration without the
//! screener (or anything holding one) being rebuilt.  ONE owner writes it; the screener only reads.
class DensityHandle
{
public:
    void         Reseat(const chmat_t* d) { itsD=d; }     //!< the iterator, once per iteration
    const chmat_t* Get() const            { return itsD; }
private:
    const chmat_t* itsD=nullptr;
};
```

`DAwareScreener` holds a `std::shared_ptr<const DensityHandle>` **at construction**.  The SCF reseats the
handle each iteration; the screener object, the task list and every box that references them are untouched.
⇒ The screener stays a POLICY object with a lifetime of the run, not a per-iteration allocation.

⚠ Two invariants to gate:
1. **Reseat happens before the first collocation of an iteration**, or the screen reads last iteration's
   \f$D\f$ — the same staleness class the XC `itsSrcVersion` guard already watches for, and it should carry
   the same kind of loud check (a serial number on the handle, compared).
2. **`Kill` and `Tolerance` must agree with each other and with the INTEGRATE side**, or collocate and
   integrate stop being adjoints.  That is what the shifted-MP defect was.  ⇒ ONE screener instance serves
   both directions; it is not two objects that happen to be configured alike.

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

| | D-aware (today) | geometry-only |
|---|---|---|
| Si Γ, collocate | 0.003296 s/call | 0.003567 s/call (**+8%**) |
| Si Γ, trajectory | 11 iters, CONVERGED | 17 iters, OSCILLATING; `Etot` agrees to **2e-9** |

⇒ The cost is single-digit percent per call.  The benefit is a large fraction of the per-line term, but it
only materialises once §5 is built.  **Do not judge the idea on the knob.**

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
