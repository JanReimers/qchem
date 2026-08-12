// File: BasisSet/Orbital_HF_IBS.C  Interface for a Hartree-Fock (HF) Orbital Irrep Basis Set.
module;
export module qchem.BasisSet.Orbital_HF_IBS;
export import qchem.BasisSet.Orbital_1E_IBS;

export namespace qchem::BasisSet
{

//--------------------------------------------------------------------------------
//
//! \brief The CONTRACTION face of a Hartree-Fock orbital irrep basis set -- everything a Fock-matrix
//! client ever needs, and nothing about how the numbers are produced.
//!
//! Client code rarely needs ERIs directly, it only needs the CONTRACTION over a density matrix, and
//! then only accumulated over many irreps.  That -- and only that -- is what this interface promises.
//! How a basis PRODUCES the contraction is its own business, and the two implementations in the tree
//! do it in genuinely different ways:
//!   - a concrete AO basis builds 4-index ERI blocks; that lineage is the substrate face
//!     \c Orbital_ERI4_IBS (BasisSet/Internal/Orbital_ERI4_IBS.C), which implements these four
//!     methods once, for every basis that has ERI4 blocks;
//!   - the SALC decorator \c SymmetryAdapted_IBS builds the whole-molecule AO Fock and SLICES its
//!     irrep block out of it -- it has no per-irrep-pair ERI4 blocks at all.
//!
//! \par Why the split (doc/CleanupCandidates.md R1.7)
//! \c MakeDirect / \c MakeExchange used to live here, so the SALC decorator -- which cannot answer
//! them -- had to supply dummy bodies returning an empty \c ERI4.  A generic ERI4 caller would then
//! get a ZERO Fock contribution with no diagnostic.  With the substrate on its own face, the SALC
//! decorator simply does not have that face, so the question cannot be asked of it: the failure moves
//! from a silent wrong answer to a compile-time (or, across an abstract cross-cast, a loud run-time)
//! error.  See also the ISP note in doc/CleanupCandidates.md: an interface should promise what its
//! clients consume, not what its most capable implementation happens to own.
//
template <class T> class Orbital_HF_IBS
    : public virtual IrrepBasisSet_IDs   // Name / BasisSetID identity face
{
public:
    //! Coulomb contraction of THIS irrep block against the \a bs_cd block's density:
    //! \f$S_{ab} \mathrel{+}= \sum_{cd} (ab|cd)\,D_{cd}\f$.  \a Dcd must be pre-screened for zero.
    virtual void AccumulateDirect  (rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* bs_cd) const=0;
    //! Exchange counterpart of AccumulateDirect: \f$S_{ab} \mathrel{+}= \sum_{cd} (ad|cb)\,D_{cd}\f$
    //! (same-spin densities only -- the caller supplies \a Dcd from one spin channel).
    virtual void AccumulateExchange(rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* bs_cd) const=0;
    //! Fused bra-ket contraction for an irrep PAIR: \f$J_i \mathrel{+}= J(i,j)\cdot D_j\f$ AND
    //! \f$J_j \mathrel{+}= J(j,i)\cdot D_i\f$ from ONE traversal of the pair.  Equivalent to the two
    //! independent AccumulateDirect calls; an implementation that exploits \f$J(j,i)=J(i,j)^\mathsf{T}\f$
    //! never builds the partner orientation (halving ERI RAM+build across the pair -- see
    //! \c Orbital_ERI4_IBS and doc/ERI4Rework.md \S4).  \a Di / \a Dj may be zero (an empty irrep).
    virtual void AccumulateDirectBoth(rsmat_t& Ji, rsmat_t& Jj, const smat_t<T>& Di, const smat_t<T>& Dj,
                                      const Orbital_HF_IBS<T>* cd) const=0;
    //! Exchange analogue of AccumulateDirectBoth (\f$K(i,j)=K(j,i)^\mathsf{T}\f$).  Same-spin densities
    //! only -- the caller (a single-spin density) supplies \a Di, \a Dj from one spin channel.
    virtual void AccumulateExchangeBoth(rsmat_t& Ki, rsmat_t& Kj, const smat_t<T>& Di, const smat_t<T>& Dj,
                                        const Orbital_HF_IBS<T>* cd) const=0;
};

//--------------------------------------------------------------------------------
//
//! \brief CAPABILITY FACE (cross-cast; absence means "not me") for a basis whose Fock is built ONCE from
//! the WHOLE-SYSTEM density and then SLICED per irrep -- instead of contracted per irrep PAIR.
//!
//! The pair loop in \c tComposite_CD::Accumulate*All exists for bases with real per-irrep-pair ERI4 blocks.
//! A SALC decorator has none (R1.7): every one of its \c AccumulateDirect calls rebuilds the SAME
//! whole-molecule AO Fock from one block's density and slices it, so the loop costs ~N whole-AO builds where
//! ONE would do -- \f$J\f$ is LINEAR in \f$D\f$, hence
//! \f$\sum_\Gamma J_{AO}(O_\Gamma D_\Gamma O_\Gamma^\mathsf{T}) =
//!    J_{AO}(\sum_\Gamma O_\Gamma D_\Gamma O_\Gamma^\mathsf{T})\f$.
//!
//! V1.31: the previous fix for that was \c SymFockCache, a memo hung off the basis that kept a copy of every
//! block's density and compared it ELEMENTWISE to decide freshness.  It was caching a PARTIAL AO Fock, at the
//! BASIS level, inside a loop that should not have been iterating -- and the elementwise compare was forced
//! by that position (down there the only thing in hand is a matrix, so no density Version() is reachable).
//! With the whole-system route the memo has nothing to memoize and is gone, along with its incomplete key.
//!
//! Three primitives, no stubs and no default bodies: a basis either has this face or does not have it.
template <class T> class WholeSystemFock_IBS
{
public:
    virtual ~WholeSystemFock_IBS() = default;
    //! Dimension of the underlying whole-system (AO) space -- the size of the matrices below.
    virtual size_t    AODimension() const=0;
    //! \f$D_{AO} \mathrel{+}= O\,D\,O^\mathsf{T}\f$: fold THIS block's density into the AO total.
    virtual void      AddAODensity(hmat_t<T>& Dao, const hmat_t<T>& D) const=0;
    //! The ONE whole-system build from the assembled AO density (\a exchange selects K over J).
    virtual hmat_t<T> MakeAOFock(const hmat_t<T>& Dao, bool exchange) const=0;
    //! \f$F_{ab} \mathrel{+}= O^\mathsf{T} F_{AO} O\f$: slice the whole-system Fock into THIS block.
    virtual void      SliceAOFock(hmat_t<T>& Fab, const hmat_t<T>& Fao) const=0;
};

typedef Orbital_HF_IBS<double>    Real_HF_OIBS;
typedef Orbital_HF_IBS<dcmplx> Complex_HF_OIBS;

} //namespace
