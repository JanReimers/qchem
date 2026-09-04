// File: BasisSet/Molecule/LatticeScreener.C  THE TOLERANCE POLICY OF THE GPW COLLOCATION BOX WALK.
//
// WHY THIS EXISTS (doc/ScreeningPlan.md; user, 2026-08-28: "I much prefer optionality to work through
// virtual dispatch rather than if statements in the low level code.  This makes the big code base much
// easier to work with in the long run").  Which tolerance sizes a collocation box, and kills a
// (pair, offset) term, used to be one experimental bool (GPW_DAWARE_SCREEN) consulted TWICE inside the
// innermost geometry code.  This is the DEPENDENCY INVERSION that replaces it: the box walk asks an
// abstract LatticeScreener for a tolerance and does not know -- must not branch on -- which rule answers.
//
// ⛔ AND THE FLAG HAD A RECORD.  The D-aware tolerance is the most defect-dense idea in the collocation
// path (ScreeningPlan.md §2 has the bug list: the 4.1 Ha shifted-MP defect that survived in a SECOND
// unreachable copy, FoldScreenMax's 2.3 Ha variational collapse, the union-vs-per-component box, the
// reduced-vs-full tier split).  A seam does not fix those, but it makes the alternative reachable,
// testable and switchable -- and it stops the next screening idea being a fifth branch in one function.
//
// STATELESS BY DESIGN, and that is the whole answer to "how does a D-aware rule reach the density?"
// (ScreeningPlan.md §4, settled 2026-09-04).  It does not: the walk already computes each term's
// screening weight for its own arithmetic and hands the block in.  So a screener is a RUN-LIFETIME
// immutable rule -- chosen once, shared by every k-block, every SCF iteration and every OpenMP worker --
// with no density reference to reseat, nothing that can go stale, and no thread-safety invariant to
// state.  The staleness class that produced this tree's three most expensive recent defects (the
// CollocMemo depth-1 thrash, the XC itsSrcVersion "STALE" warning, the cross-run Projector3 pollution)
// is UNREPRESENTABLE here rather than guarded against.
module;
#include <cassert>
#include <cstdlib>   // std::getenv/std::atof (the GPW_DENSITY_EPS floor)
#include <vector>
export module qchem.BasisSet.Molecule.LatticeScreener;
import qchem.Types;      // ivec3_t, rmat_t
import qchem.Blaze;      // the rmat_t block's element access
import qchem.Math;       // log/fabs/max/min (the umbrella re-exports them UNQUALIFIED -- CLAUDE.md)

export namespace qchem::BasisSet::Molecule
{

//! \brief THE COLLOCATION TOLERANCE FLOOR, in ONE place.
//!
//! Every screener is built at it, no screener may answer below it, and the (shell pair, offset) task list
//! is enumerated at it -- which is what makes that list a strict SUPERSET of any iteration's live set.
//! It lives here, with the screening policy, rather than in the walk: two copies of this number is exactly
//! the drift the seam exists to prevent.  (The walk's \c kDensityEps() delegates to it.)
//! \note DECOUPLED from the analytic 1E screen \c GPW_SCREEN_EPS -- a collocated density and an analytic
//! lattice sum are not converged by the same tolerance.
inline double CollocationEps()
{ static const double e=[]{const char* s=std::getenv("GPW_DENSITY_EPS"); return s?std::atof(s):1e-10;}(); return e; }

//! \brief The screen's verdict for ONE (shell pair, offset) task: how big the shared box must be, and
//! which of the shell pair's component terms ride it.
//!
//! Caller-owned and REUSED across the offsets of a shell pair, so screening a task allocates nothing
//! after the first.  The mask reads as boolean through \c Keeps(); \c char is only how it is stored --
//! see the note on \c itsKeep.
class ScreenPlan
{
public:
    //! The tolerance the ONE shared box must honour.  0.0 = the whole task is dead; skip it.
    double Tolerance()                    const {return itsEps;}
    bool   Dead()                         const {return itsEps==0.0;}
    //! Does component \a (a,b) of the \f$n_I\times n_J\f$ block ride the shared box?
    bool   Keeps(size_t a, size_t b)      const {assert(a*itsNJ+b<itsKeep.size()); return itsKeep[a*itsNJ+b]!=0;}

    //! \name The screener's side of the contract.  Not for callers.
    //!@{
    //! Size the mask to an \a nI x \a nJ block and clear it.  Reuses the buffer when the shape repeats.
    void Reset(size_t nI, size_t nJ)
    {
        itsNJ=nJ; itsEps=0.0;
        itsKeep.assign(nI*nJ, char(0));
    }
    void Keep(size_t a, size_t b)               {assert(a*itsNJ+b<itsKeep.size()); itsKeep[a*itsNJ+b]=char(1);}
    void SetTolerance(double eps)               {itsEps=eps;}
    //!@}

private:
    //! ⚠ NOT \c std::vector<bool>.  That specialization is bit-packed and is not a container of bool: it
    //! hands out a PROXY rather than a \c bool&, so \c auto silently captures the proxy; it has no
    //! \c data(); and neighbouring elements share a word, so it cannot be written elementwise from
    //! parallel workers.  This walk runs under OpenMP.  The mask reads as boolean through \c Keeps() --
    //! \c char is storage, not interface.
    std::vector<char> itsKeep;
    size_t            itsNJ=0;
    double            itsEps=0.0;
};

//! \brief The tolerance rule for one (shell pair, offset) collocation task.  STATELESS and RUN-LIFETIME;
//! see the file header for why it holds no density.
class LatticeScreener
{
public:
    virtual ~LatticeScreener() = default;

    //! \brief Screen one task's component terms and size the ONE box they must share.
    //!
    //! \param n     the integer cell offset -- a key; neither of today's two rules reads it.
    //! \param pf    the pair prefactor exponent (screen 1), a shell-pair property of this offset: the
    //!              whole box is below \f$\varepsilon\f$ when \f$pf\ge-\ln\varepsilon\f$.
    //! \param cij   the \f$n_I\times n_J\f$ block of screening weights over the shell pair's components.
    //!              ⚠ \b NOT \f$D_{ij}\f$.  It is \f$c_{ij}=\mathrm{fold}\cdot\mathrm{Re}[D_{ij}e^{-ikR}]\f$
    //!              collocating and \f$\mathrm{fold}\cdot|D_{ij}|\f$ gathering (ScreeningPlan.md §1) --
    //!              the two directions screen on DIFFERENT quantities, deliberately and documented, see
    //!              \c LatticeSum1E::IntegratePotential.  0 = the term is absent (dead under the stream
    //!              fold, or a zero density element).
    //! \param plan  OUT.  \c Tolerance() is the UNION tolerance -- the TIGHTEST survivor's, i.e. the
    //!              LARGEST box -- or 0.0 when nothing survives.
    //!
    //! \note THE UNION IS THE SCREENER'S DECISION, not the caller's.  One box must serve a whole shell
    //! pair, so it is sized by the tightest \f$\varepsilon\f$ present and each survivor rides a box at
    //! least as large as its own rule asked for.  That only ever ADDS sub-eps terms, never drops one, so
    //! it keeps the no-cut discipline (doc/GPWPlan.md pin).  It was open-coded identically at both call
    //! sites before this seam existed -- which is exactly the kind of duplicated decision §2's bug list
    //! is made of.
    virtual void Screen(const ivec3_t& n, double pf, const rmat_t& cij, ScreenPlan& plan) const = 0;

    //! \brief The global floor: no \c Screen may return a tolerance below it, so the STATICALLY screened
    //! (shell pair, offset) task list stays a strict SUPERSET of every iteration's live set.  Asserted at
    //! the walk (see \c NR_Evaluator::BoxTasks, whose list is built at this floor).
    virtual double Floor() const = 0;
};

//! \brief CP2K's rule: a flat \f$\varepsilon\f$, the weights ignored.
//!
//! Every task's answer is the same one, so the whole box geometry is ITERATION-INVARIANT and hoistable
//! into the task list (ScreeningPlan.md §5 -- the prize this seam exists to make reachable).  CP2K screens
//! this way: \c task_list_methods.F takes no density argument, its eps is the global \c eps_rho_rspace and
//! \c radius_list is fixed task-list data (checked in the source, 2026-08-28).
//!
//! \note IT STILL KILLS.  \a pf against the flat \f$\varepsilon\f$ IS the geometry screen -- CP2K's fixed
//! radius test.  "Geometry only" is not "no screening".
class GeometryOnlyScreener : public LatticeScreener
{
public:
    explicit GeometryOnlyScreener(double eps) : itsEps(eps), itsLnEps(-log(eps)) {assert(eps>0.0);}

    virtual void Screen(const ivec3_t&, double pf, const rmat_t& cij, ScreenPlan& plan) const
    {
        const size_t nI=cij.rows(), nJ=cij.columns();
        plan.Reset(nI,nJ);
        if (pf>=itsLnEps) return;                      // the whole box is sub-eps for EVERY term: dead
        bool any=false;
        for (size_t a=0;a<nI;a++)
            for (size_t b=0;b<nJ;b++)
                if (cij(a,b)!=0.0) {plan.Keep(a,b); any=true;}
        if (any) plan.SetTolerance(itsEps);
    }
    virtual double Floor() const {return itsEps;}

private:
    const double itsEps, itsLnEps;
};

//! \brief The D-aware rule (this tree's default, and a declared CP2K deviation -- \c RunPolicy).
//!
//! What lands on the grid is \f$c_{ij}\chi_i\chi_j\f$ and the accuracy target is ABSOLUTE on the density,
//! so the tolerance a term's box must honour is \f$\varepsilon/|c_{ij}|\f$, not \f$\varepsilon\f$: a
//! small-\f$|c|\f$ term keeps a SMALLER box (its radius shrinks with \f$\sqrt{\ln}\f$, its work with the
//! cube), and a term whose whole box is below its own tolerance is dropped.  Clamped at the floor so a
//! \f$|c|>1\f$ never grows a box BEYOND the geometry screen -- the task list is built at the floor and
//! replay is capped by what was built.
class DAwareScreener : public LatticeScreener
{
public:
    explicit DAwareScreener(double eps) : itsEps(eps) {assert(eps>0.0);}

    virtual void Screen(const ivec3_t&, double pf, const rmat_t& cij, ScreenPlan& plan) const
    {
        const size_t nI=cij.rows(), nJ=cij.columns();
        plan.Reset(nI,nJ);
        double epsUnion=0.0;
        for (size_t a=0;a<nI;a++)
            for (size_t b=0;b<nJ;b++)
            {
                const double c=cij(a,b);
                if (c==0.0) continue;                                   // absent: no weight, no term
                const double e=max(itsEps, itsEps/fabs(c));
                if (pf>=-log(e)) continue;                         // whole box below THIS term's eps
                plan.Keep(a,b);
                epsUnion = epsUnion==0.0 ? e : min(epsUnion,e);
            }
        plan.SetTolerance(epsUnion);
    }
    virtual double Floor() const {return itsEps;}

private:
    const double itsEps;
};

}
