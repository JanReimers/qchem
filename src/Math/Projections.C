// File: Projections.C  A density's projections <c|rho> onto ONE basis's functions: a one-of-N record.
//
// ONE operation, "this density's projections onto that basis", has two realizations chosen by the PRODUCER (a
// MatrixForward): an rvec_t (<f_a|rho> on a Gaussian basis; rho(r_g) on a delta basis or a mesh; occupation blocks for
// a Hubbard projector) or a ΔG_Map (rho-tilde(G) on a plane-wave basis).  The house shape for that (CLAUDE.md "One
// result, several realizations", user 2026-11-07) is a record holding every realization with exactly one populated --
// not a class hierarchy, not a std::variant with a visit in every client.
//
// ★ PROJECTIONS, NOT COEFFICIENTS (user rulings 2026-11-08, doc/Records/VOneOfN_D2Design.md (a), (b)).  The record is
// the forward's OUTPUT only.  A fitter's coefficients are its own private type, and the metric solve (S^-1, G^-1,
// Dunlap, identity) is inside the fitter; the record never carries them.
//
// ★ TWO LEVELS OF CLIENT (user ruling 2026-11-08, doc/Records/VOneOfNDesign.md §3a).
//   - Everything that only hands projections on (a density summing its blocks, a mixer) uses the DISPATCHING
//     operations only: +=, -=, *=, MaxAbs.  It never asks which realization it holds.
//   - The fitter that consumes its own projections (it built the projector, so it knows the realization) reads them
//     with the typed accessors Vector() / GMap().  A wrong one THROWS.
//   So `grep '\.Vector()\|\.GMap()'` lists exactly the realization-aware code.
//
// It lives in qcMath, below qchem.Mesh.Integrator, because MatrixForward returns it and PairProjection (qcBasisSet)
// imports that.  ΔG_Map and its field algebra live here with it (qchem.BasisSet.PairProjection re-exports this module).
module;
#include <algorithm>   // std::max
#include <cmath>       // std::abs
#include <complex>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
export module qchem.Projections;
export import qchem.Types;   // ivec3_t, dcmplx
import qchem.Blaze;          // rvec_t + its arithmetic

export namespace qchem {

//! Lexicographic comparator so a reciprocal-index triple \f$\Delta m\f$ can key a map.
struct IVec3Less
{
    bool operator()(const ivec3_t& a, const ivec3_t& b) const
    { if (a.x!=b.x) return a.x<b.x; if (a.y!=b.y) return a.y<b.y; return a.z<b.z; }
};

//! G-space components (density \f$\tilde\rho\f$ or potential \f$\tilde V\f$) keyed by the reciprocal-index
//! difference \f$\Delta m\f$ (\f$\Delta G = B\,\Delta m\f$).
using ΔG_Map = std::map<ivec3_t, dcmplx, IVec3Less>;

//---------------------------------------------------------------------------------------------------------
//  THE FIELD ALGEBRA of a ΔG_Map -- a sparse G-space field with the usual linear-space operations.  Found by
//  ADL through IVec3Less (namespace qchem), so a caller need only import this module.  A key absent from one
//  operand reads as 0, so `a - b` spans the UNION of the two index sets (what a residual ‖out−in‖ needs).
//  (V1.18: these were seven Map* helpers private to the density mixers; a mixer's arithmetic should READ as
//  arithmetic, and the type's algebra belongs with the type.)
//---------------------------------------------------------------------------------------------------------
inline ΔG_Map& operator+=(ΔG_Map& a, const ΔG_Map& b) { for (const auto& [k,v] : b) a[k]+=v; return a; }
inline ΔG_Map& operator-=(ΔG_Map& a, const ΔG_Map& b) { for (const auto& [k,v] : b) a[k]-=v; return a; }
inline ΔG_Map& operator*=(ΔG_Map& a, double s)         { for (auto& [k,v] : a) v=s*v; return a; }
inline ΔG_Map  operator+ (ΔG_Map a, const ΔG_Map& b)   { return a+=b; }
inline ΔG_Map  operator- (ΔG_Map a, const ΔG_Map& b)   { return a-=b; }
inline ΔG_Map  operator* (double s, ΔG_Map a)          { return a*=s; }
inline ΔG_Map  operator* (ΔG_Map a, double s)          { return a*=s; }
//! \f$\max_G|f(G)|\f$ -- the ∞-norm (a mixer's convergence gate).
inline double  MaxAbs(const ΔG_Map& m) { double x=0.0; for (const auto& [k,v] : m) x=std::max(x,std::abs(v)); return x; }
//! \f$\sum_i c_i f_i\f$ over any range of fields (a DIIS/Pulay extrapolation).
template <class Range> inline ΔG_Map LinearCombination(const Range& fields, const rvec_t& c)
{
    ΔG_Map r; size_t i=0;
    for (const auto& f : fields) { for (const auto& [k,v] : f) r[k]+=c[i]*v; ++i; }
    return r;
}

//! \brief A density's projections onto ONE basis's functions.  Exactly one realization is populated,
//! chosen by the PRODUCER (a \c MatrixForward, a fitter); see the file header for who may read which.
//!
//! A default-constructed record is EMPTY: the additive identity, holding no realization yet.  `+=` adopts
//! the other operand's realization, so a sum over blocks starts from `{}` without the summing code knowing
//! which realization its blocks produce.
//!
//! NON-template (user ruling 2026-11-07): an empty member is cheap, and one type serves every run scalar.
class Projections
{
public:
    Projections() = default;                                        //!< empty: the additive identity
    explicit Projections(rvec_t a) : itsKind(Kind::Vector), itsA(std::move(a)) {} //!< a Gaussian fit's c_a, a δ fit's values at its points, a mesh's values
    explicit Projections(ΔG_Map g) : itsKind(Kind::GMap),   itsG(std::move(g)) {} //!< a plane-wave basis's projections ρ̃(G)

    //! \name The PAIRED implementations' typed reads -- each THROWS on a mismatch (see the file header).
    //! The rvalue overloads move the member out, so `const rvec_t& v=f.Forward(D).Vector();` binds to a
    //! value instead of dangling into a dead temporary.
    //@{
    const rvec_t& Vector() const & { RequireKind(Kind::Vector, "Vector()"); return itsA; }
    rvec_t        Vector()      && { RequireKind(Kind::Vector, "Vector()"); return std::move(itsA); }
    const ΔG_Map& GMap()   const & { RequireKind(Kind::GMap,   "GMap()");   return itsG; }
    ΔG_Map        GMap()        && { RequireKind(Kind::GMap,   "GMap()");   return std::move(itsG); }
    //@}

    //! \name Dispatching operations: every client may use these.  A realization mismatch THROWS.
    //@{
    Projections& operator+=(const Projections& b) { return Accumulate(b, +1.0, "+="); }
    Projections& operator-=(const Projections& b) { return Accumulate(b, -1.0, "-="); }
    Projections& operator*=(double s)
    {
        if (itsKind==Kind::Vector) itsA*=s;
        else if (itsKind==Kind::GMap) itsG*=s;
        return *this;
    }
    //! The ∞-norm, \f$\max_a|c_a|\f$ (0 for an empty record).
    double MaxAbs() const
    {
        if (itsKind==Kind::GMap) return qchem::MaxAbs(itsG);
        double x=0.0;
        for (double v : itsA) x=std::max(x,std::abs(v));
        return x;
    }
    //@}

private:
    //! Which member is populated.  A private tag rather than "which member is non-empty", because an
    //! EMPTY vector (a block with no fit functions) is still the vector realization.
    enum class Kind { None, Vector, GMap };

    void RequireKind(Kind k, const char* what) const
    {
        if (itsKind==k) return;
        throw std::logic_error(std::string("Projections::")+what+": this record holds "+Name(itsKind)+
            ", not "+Name(k)+".  Only the paired implementation, which knows its realization, reads a member; "
            "everything else uses the dispatching operations.");
    }
    static const char* Name(Kind k)
    {
        return k==Kind::Vector ? "a coefficient vector (rvec_t)" : k==Kind::GMap ? "a G-space map (ΔG_Map)" : "nothing (empty)";
    }
    Projections& Accumulate(const Projections& b, double sign, const char* what)
    {
        if (b.itsKind==Kind::None) return *this;
        if (itsKind==Kind::None) { *this=b; if (sign<0) *this*=-1.0; return *this; }
        if (itsKind!=b.itsKind)
            throw std::logic_error(std::string("Projections::operator")+what+": cannot combine "+Name(itsKind)+
                                   " with "+Name(b.itsKind));
        if (itsKind==Kind::Vector)
        {
            if (itsA.size()!=b.itsA.size())
                throw std::logic_error(std::string("Projections::operator")+what+": coefficient vectors of different lengths");
            if (sign>0) itsA+=b.itsA; else itsA-=b.itsA;
        }
        else
        {
            if (sign>0) itsG+=b.itsG; else itsG-=b.itsG;
        }
        return *this;
    }

    Kind   itsKind=Kind::None;
    rvec_t itsA;   //!< populated iff itsKind==Vector
    ΔG_Map itsG;   //!< populated iff itsKind==GMap
};

} // namespace qchem
