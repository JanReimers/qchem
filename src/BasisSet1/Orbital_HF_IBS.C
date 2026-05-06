// File: BasisSet/Orbital_HF_IBS.C  Interface for a Hartree-Fock (HF) Orbital Irrep Basis Set.
module;
export module qchem.BasisSet1.Orbital_HF_IBS;
export import qchem.BasisSet1.Orbital_1E_IBS;
export import qchem.BasisSet.Internal.ERI4;

export namespace BasisSet1
{

template <class T> class Orbital_HF_IBS
    : public virtual Orbital_1E_IBS<T>
    , public virtual IrrepBasisSet_IDs //avoid using statements for RadialID,AngularID
{
public:
    virtual ERI4       MakeDirect  (const Orbital_HF_IBS<T>& c) const=0; //Only called once for a given {radial,angular} ID pair.
    virtual ERI4       MakeExchange(const Orbital_HF_IBS<T>& c) const=0; //Only called once for a given {radial,angular} ID pair.
    const   ERI4&          Direct  (const Orbital_HF_IBS<T>& c) const;
    const   ERI4&          Exchange(const Orbital_HF_IBS<T>& c) const;
    // Client code rarely need ERIs directly, they only need the contraction over a density matrix.
    // In addition they only need the accumulation of these contractions over many Irreps.
    void AccumulateDirect  (rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* bs_cd) const;
    void AccumulateExchange(rsmat_t& Sab, const smat_t<T>& Dcd, const Orbital_HF_IBS<T>* bs_cd) const;
};

typedef Orbital_HF_IBS<double>    Real_HF_OIBS;
typedef Orbital_HF_IBS<dcmplx> Complex_HF_OIBS;

} //namespace
