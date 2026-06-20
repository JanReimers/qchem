// File: BasisSet/Molecule/Evaluators/Evaluator.C
//
// Molecular evaluator interfaces.  COPIED from BasisSet/Atom/Evaluators/Evaluator.C and left free to
// diverge (the "copy, don't share" approach): once the molecular evaluators work end to end we can see
// what genuinely remains common with the atom ones and lift that out -- guessing the shared abstraction
// up front is what we are deliberately avoiding.
//
// Known divergences from the atom version so far:
//  - NO Getl():  an atom radial carries a single l; a molecular Cartesian shell spreads over components
//    each with its own (lx,ly,lz), so "the l" is not a property of the evaluator.
//  - NO Angular / DirectAk / ExchangeAk:  that is the atomic angular (Ak) factorization of the ERI.
//    Molecular HF uses the 4-centre McMurchie-Davidson integrals instead.
//  - Kinetic:  the building block is the FULL Cartesian \f$-\nabla^2\f$ (see Grad2 below) -- there is no
//    separate centrifugal \f$l(l+1)\langle r^{-2}\rangle\f$ term (that is an atom radial/angular split).
//  - NOT (yet) a VectorFunction<double>:  the atom Evaluator is one so it can be evaluated on the DFT
//    mesh via operator()(r).  Molecular Cartesian grid-eval belongs with the DFT (3-centre) increment,
//    so for the 1E-first cut the base is just Streamable; VectorFunction + operator()(r) come back then.
//
// As in the atom design: virtual functions for the cold paths, C++20 concepts for the inline hot-loop
// kernels.  Start with 1E (Overlap/Kinetic/Nuclear); DFT (3-centre) and HF (4-centre) concepts come
// later and will diverge further (no Cache4/Grouper -- molecular caching is Cache2/Cache3).
module;
#include <iosfwd>
#include <concepts>
export module qchem.BasisSet.Molecule.Evaluators;
export import qchem.Streamable;
import qchem.Types;
import qchem.Blaze;   // rsmat_t construction/indexing in the generic matrix builders below

export namespace BasisSet::Molecule::Evaluators
{

// Abstract interface common to all molecular evaluators.  The virtuals here are only touched outside
// hot loops (sizing, identity, streaming); the per-element integral kernels are the concepts below.
class Evaluator
    : public virtual Streamable
{
public:
    virtual size_t  size    () const = 0;                  // number of basis functions (components)
    virtual size_t  maxSpan () const {return size();}      // assume no overlap for |i-j| > maxSpan
    virtual rvec_t  Norm    () const = 0;                  // per-component normalization constants
    virtual std::ostream& Write    (std::ostream&) const=0;
    virtual std::string   RadialID () const=0;             // key fragment for the integral cache
    virtual std::string   Name     () const=0;
    // Helpers for range-based loops over the i,j matrix indices.
    iv_t indices(             ) const {return iv_t(size_t(0),size());}
    iv_t indices(size_t start ) const {return iv_t(start    ,size());}  // 2nd index of symmetric matrices
};

// We probably want every evaluator to support this (plotting + numerical-integration unit tests); for a
// ground-state calculation it is only required by the DFT/Fit evaluators.
template <class E> concept isOpr_Evaluator = requires (E e, const rvec3_t& r)
{
    {e.operator()(r)} -> std::same_as<rvec_t>;
    {e.Gradient  (r)} -> std::same_as<rvec3vec_t>;
};

// One-electron inline integrals common to all Hamiltonians (the hot loops).
//
// Grad2(i,j) is the KINETIC BUILDING BLOCK \f$\langle p^2\rangle = \langle i|-\nabla^2|j\rangle\f$ --
// the FULL Cartesian negative-Laplacian (positive-definite), with NO 1/2 and NO centrifugal term.
// The 1/2 of \f$T=-\tfrac12\nabla^2\f$ is applied at the Hamiltonian boundary, NOT here
// (see BasisSet/Orbital_1E_IBS.C).  Nuclear(i,j) is the full multi-centre attraction
// \f$\sum_c -Z_c\langle i|\,|r-R_c|^{-1}|j\rangle\f$ (the evaluator carries the cluster).
template <class E> concept is1E_Evaluator = std::derived_from<E, Evaluator> && requires (E e,size_t i,size_t j)
{
    {e.Norm    (i)  } -> std::same_as<double>;
    {e.Overlap (i,j)} -> std::same_as<double>;
    {e.Grad2   (i,j)} -> std::same_as<double>;   // <p^2> = <-nabla^2>, full Cartesian, no 1/2
    {e.Nuclear (i,j)} -> std::same_as<double>;   // multi-centre sum_c -Z_c/|r-R_c|
};

// --- Generic 1E matrix builders -------------------------------------------------------------------
// The basis-set-agnostic i,j matrix-build loops, driven purely by the evaluator's inline kernels.
// (Mirrors the atom Integrals_Overlap<E>/Integrals_Kinetic<E> mixins; these are the natural candidates
// to hoist to the generic BasisSet level later -- plan Goal C.)  They fill the upper triangle of a
// symmetric matrix; rsmat_t mirrors the lower.
//
// KineticMatrix returns the kinetic BUILDING BLOCK \f$\langle p^2\rangle=\langle-\nabla^2\rangle\f$,
// i.e. just the sum of Grad2 -- NO 1/2 (the Hamiltonian applies it) and NO centrifugal term (the
// molecular Grad2 is already the full Cartesian \f$-\nabla^2\f$).  See BasisSet/Orbital_1E_IBS.C.
template <is1E_Evaluator E> rsmat_t OverlapMatrix(const E& e)
{
    rsmat_t S(e.size());
    for (auto i:e.indices()) for (auto j:e.indices(i)) S(i,j)=e.Overlap(i,j);
    return S;
}
template <is1E_Evaluator E> rsmat_t KineticMatrix(const E& e)   // <p^2>=<-nabla^2> block, no 1/2
{
    rsmat_t S(e.size());
    for (auto i:e.indices()) for (auto j:e.indices(i)) S(i,j)=e.Grad2(i,j);
    return S;
}
template <is1E_Evaluator E> rsmat_t NuclearMatrix(const E& e)
{
    rsmat_t S(e.size());
    for (auto i:e.indices()) for (auto j:e.indices(i)) S(i,j)=e.Nuclear(i,j);
    return S;
}

// TODO (later increments, will diverge from atom's Cache4/Angular design):
//  - isFit_Evaluator  : Charge/Overlap/Repulsion for auxiliary fit basis sets.
//  - isDFT_Evaluator  : 3-centre Overlap/Repulsion (DFT), via the molecular Cache2/Cache3.
//  - isHF_Evaluator   : 4-centre Direct/Exchange via McMurchie-Davidson (no Cache4/Grouper).

} //namespace
