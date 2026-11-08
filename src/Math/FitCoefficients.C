// File: FitCoefficients.C  A density's (or a field's) coefficients in ONE fit basis: a one-of-N record.
//
// ONE operation, "these coefficients in that fit basis", has two realizations chosen by the PRODUCER: an
// rvec_t (a Gaussian fit's c_a; a delta fit's values at its points; a mesh's values) or a ΔG_Map (a
// plane-wave fit's rho-tilde(G)).  The house shape for that (CLAUDE.md "One result, several realizations",
// user 2026-11-07; doc/Records/VOneOfNDesign.md) is a record holding every realization with exactly one
// populated -- not a class hierarchy, not a std::variant with a visit in every client.
//
// ★ TWO LEVELS OF CLIENT (user ruling 2026-11-08, VOneOfNDesign.md §3a).
//   - Everything ABOVE the paired implementations (the SCF, the mixers, the Hamiltonian terms, a density
//     that only hands coefficients on) uses the DISPATCHING operations only: +=, -=, *=, MaxAbs.  It never
//     asks which realization it holds.
//   - The PAIRED implementations -- the code that built (or was handed) the producer and so already knows
//     the realization: the FunctionFitter derived classes, their samplers, Hubbard_U's integrator, the
//     integrators' own adjoints -- read it with the typed accessors Vector() / GMap().  A wrong one THROWS.
//   So `grep '\.Vector()\|\.GMap()'` lists exactly the realization-aware code.
//
// It lives in qcMath, below qchem.Mesh.Integrator, because MatrixForward / MatrixAdjoint carry it and
// PairProjection (qcBasisSet) imports that.  ΔG_Map and its field algebra moved here with it (they were in
// qchem.BasisSet.PairProjection, which re-exports this module, so their spelling is unchanged).
module;
#include <algorithm>   // std::max
#include <cmath>       // std::abs
#include <complex>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
export module qchem.FitCoefficients;
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

//! \brief A density's (or a field's) coefficients in ONE fit basis.  Exactly one realization is populated,
//! chosen by the PRODUCER (a \c MatrixForward, a fitter); see the file header for who may read which.
//!
//! A default-constructed record is EMPTY: the additive identity, holding no realization yet.  `+=` adopts
//! the other operand's realization, so a sum over blocks starts from `{}` without the summing code knowing
//! which realization its blocks produce.
//!
//! NON-template (user ruling 2026-11-07): an empty member is cheap, and one type serves every run scalar.
class FitCoefficients
{
public:
    FitCoefficients() = default;                                        //!< empty: the additive identity
    explicit FitCoefficients(rvec_t a) : itsKind(Kind::Vector), itsA(std::move(a)) {} //!< a Gaussian fit's c_a, a δ fit's values at its points, a mesh's values
    explicit FitCoefficients(ΔG_Map g) : itsKind(Kind::GMap),   itsG(std::move(g)) {} //!< a plane-wave fit's coefficients ρ̃(G)

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
    FitCoefficients& operator+=(const FitCoefficients& b) { return Accumulate(b, +1.0, "+="); }
    FitCoefficients& operator-=(const FitCoefficients& b) { return Accumulate(b, -1.0, "-="); }
    FitCoefficients& operator*=(double s)
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
        throw std::logic_error(std::string("FitCoefficients::")+what+": this record holds "+Name(itsKind)+
            ", not "+Name(k)+".  Only the paired implementation, which knows its realization, reads a member; "
            "everything else uses the dispatching operations.");
    }
    static const char* Name(Kind k)
    {
        return k==Kind::Vector ? "a coefficient vector (rvec_t)" : k==Kind::GMap ? "a G-space map (ΔG_Map)" : "nothing (empty)";
    }
    FitCoefficients& Accumulate(const FitCoefficients& b, double sign, const char* what)
    {
        if (b.itsKind==Kind::None) return *this;
        if (itsKind==Kind::None) { *this=b; if (sign<0) *this*=-1.0; return *this; }
        if (itsKind!=b.itsKind)
            throw std::logic_error(std::string("FitCoefficients::operator")+what+": cannot combine "+Name(itsKind)+
                                   " with "+Name(b.itsKind));
        if (itsKind==Kind::Vector)
        {
            if (itsA.size()!=b.itsA.size())
                throw std::logic_error(std::string("FitCoefficients::operator")+what+": coefficient vectors of different lengths");
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
