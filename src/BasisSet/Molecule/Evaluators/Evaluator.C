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
import qchem.Cluster;  // Cluster* threaded through the Nuclear kernel

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
template <class E> concept is1E_Evaluator = std::derived_from<E, Evaluator> && requires (E e,size_t i,size_t j,const Cluster* cl)
{
    {e.Norm    (i)     } -> std::same_as<double>;
    {e.Overlap (i,j)   } -> std::same_as<double>;
    {e.Grad2   (i,j)   } -> std::same_as<double>;   // <p^2> = <-nabla^2>, full Cartesian, no 1/2
    {e.Nuclear (i,j,cl)} -> std::same_as<double>;   // multi-centre sum_c -Z_c/|r-R_c|; cluster passed per call
};

// --- 3-centre (DFT) and 4-centre (HF) inline kernels ----------------------------------------------
// These mirror the atom isDFT_Evaluator / isHF_Evaluator concepts, but the molecular kernels are
// *multi-evaluator*: each (evaluator, index) pair names one basis component, so the same kernel serves
// Coulomb and Exchange (the caller maps its loop indices into the slots).  `e` is the A slot; the
// remaining evaluators are passed as same-type references (a member may read another evaluator's data).
//
// Divergence from atom: NO isOpr_Evaluator requirement here.  The molecular DFT 3-centre integrals are
// analytic (M&D ThreeC); numerical grid-eval for the XC potential still lives on the IBS (operator()(r)),
// not the evaluator -- so requiring isOpr would (correctly) reject PG_Evaluator today.  When grid-eval
// moves onto the evaluator, fold isOpr<E> back into isDFT_Evaluator to match the atom side.
template <class E> concept isDFT_Evaluator = std::derived_from<E, Evaluator>
    && requires (const E e, size_t iA, size_t iB, size_t iC)
{
    {e.OverlapThreeC  (iA, e, iB, e, iC)} -> std::same_as<double>;   // <ab|c> overlap,   M&D 3-centre
    {e.RepulsionThreeC(iA, e, iB, e, iC)} -> std::same_as<double>;   // <ab|c> repulsion, M&D 3-centre
};

// 4-centre electron-repulsion (ab|cd), M&D 4-centre.  No Cache4/Grouper (that is the atom Ak design);
// molecular caching is the Cache2/Cache3 work of plan Stages 2-3, layered on later.
template <class E> concept isHF_Evaluator = std::derived_from<E, Evaluator>
    && requires (const E e, size_t iA, size_t iB, size_t iC, size_t iD)
{
    {e.FourC(iA, e, iB, e, iC, e, iD)} -> std::same_as<double>;   // (ab|cd)
};

// Composite the molecular Orbital_IBS will instantiate against (it IS-A 1E + DFT + HF).
template <class E> concept is1E_DFT_HF_Evaluator = is1E_Evaluator<E> && isDFT_Evaluator<E> && isHF_Evaluator<E>;

// The basis-agnostic i,j / 3C / 4C matrix-build loops that consume these kernels live in the templated
// IBS mixins (qchem.BasisSet.Molecule.IBS), mirroring the atom Integrals_*<E> / Orbital_*_IBS<E> pattern.
//
// TODO (later increment): isFit_Evaluator (Charge/Overlap/Repulsion for auxiliary fit basis sets).

} //namespace
