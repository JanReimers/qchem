// File: BasisSet/Orbital_HF_IBS.C  Interface for a Hartree-Fock (HF) Orbital Irrep Basis Set.
module;
export module qchem.BasisSet.Orbital_HF_IBS;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Internal.ERI4;

export namespace qchem::BasisSet
{

template <class T> class Orbital_HF_IBS
    : public virtual IrrepBasisSet_IDs //avoid using statements for RadialID,AngularID
{
public:
    virtual ERI4       MakeDirect  (const Orbital_HF_IBS<T>& c) const=0; //Only called once for a given {radial,angular} ID pair.
    virtual ERI4       MakeExchange(const Orbital_HF_IBS<T>& c) const=0; //Only called once for a given {radial,angular} ID pair.
    const   ERI4&          Direct  (const Orbital_HF_IBS<T>& c) const;
    const   ERI4&          Exchange(const Orbital_HF_IBS<T>& c) const;
    // Client code rarely need ERIs directly, they only need the contraction over a density matrix.
    // In addition they only need the accumulation of these contractions over many Irreps.
    // Virtual so a SALC decorator can build these in the AO basis and slice to its irrep block.
    virtual void AccumulateDirect  (rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* bs_cd) const;
    virtual void AccumulateExchange(rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* bs_cd) const;
    //! Fused bra-ket scatter (doc/ERI4Rework.md \S4): fetch the SINGLE canonical Direct block J(this,cd)
    //! and feed BOTH Fock sub-blocks in one pass -- \f$J_i \mathrel{+}= J\cdot D_j\f$ (localized) and
    //! \f$J_j \mathrel{+}= J^\mathsf{T}\cdot D_i\f$ (scaled whole-block add) -- so J(cd,this) is never
    //! built or stored.  Equivalent to the two independent AccumulateDirect contractions, halving ERI
    //! RAM+build across an irrep pair.  All the symmetry bookkeeping lives in ERI4::ScatterBoth.
    virtual void AccumulateDirectBoth(rsmat_t& Ji, rsmat_t& Jj, const smat_t<T>& Di, const smat_t<T>& Dj,
                                      const Orbital_HF_IBS<T>* cd) const;
};

typedef Orbital_HF_IBS<double>    Real_HF_OIBS;
typedef Orbital_HF_IBS<dcmplx> Complex_HF_OIBS;

} //namespace
