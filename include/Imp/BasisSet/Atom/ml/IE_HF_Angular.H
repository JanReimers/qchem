// File: Atom_ml_IE_HF_Angular.H  Angular 2e-Integrals for atoml HF basis sets.
#ifndef _Atom_ml_IE_HF_Angular_H_
#define _Atom_ml_IE_HF_Angular_H_

#include "../IE_HF_Angular.H"

namespace Atom_ml
{

class IE_BS_2E_Angular : public virtual ::AtomIE_BS_2E_Angular
{
public:
    typedef ::AtomIE_BS_2E_Angular::RVec RVec;
    virtual RVec Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
    virtual RVec ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const;
};

} //namespace

#endif //_Atom_ml_IE_HF_Angular_H_