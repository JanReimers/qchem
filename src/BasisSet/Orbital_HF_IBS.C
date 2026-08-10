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

typedef Orbital_HF_IBS<double>    Real_HF_OIBS;
typedef Orbital_HF_IBS<dcmplx> Complex_HF_OIBS;

} //namespace
