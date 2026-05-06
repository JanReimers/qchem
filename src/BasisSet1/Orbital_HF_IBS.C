// File: BasisSet/Orbital_HF_IBS.C  Interface for a Hartree-Fock (HF) Orbital Irrep Basis Set.
module;
#include <cassert>
export module qchem.BasisSet1.Orbital_HF_IBS;
export import qchem.BasisSet1.Orbital_1E_IBS;
import qchem.BasisSet1.DB_Cache;
export import qchem.BasisSet.Internal.ERI4;

export namespace BasisSet1
{



template <class T> class Orbital_HF_IBS
    : public virtual IrrepBasisSet_IDs
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

//! \brief Interface for 4-center ERI integrals used in HF calculations.
//! This particular interface if for serving up ERIs between two Irrep Basis Sets (IRBs).
//! 
// template <class T> class Integrals_BS_HF
// {
// public:
//     virtual ERI4 Direct  (const Orbital_HF_IBS<T>& a,const Orbital_HF_IBS<T>& c) const=0;
//     virtual ERI4 Exchange(const Orbital_HF_IBS<T>& a,const Orbital_HF_IBS<T>& b) const=0;
// };

} //namespace
