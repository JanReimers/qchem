// File: Atoml_IE_HF_Angular.H  Angular 2e-Integrals for atoml HF basis sets.
#ifndef _Atoml_IE_HF_Angular_H_
#define _Atoml_IE_HF_Angular_H_

import qchem.BasisSet.Atom.IE;
namespace Atoml
{

class IE_BS_2E_Angular : public virtual ::AtomIE_BS_2E_Angular
{
public:
    virtual RVec Coulomb_AngularIntegrals(const iec_t* a,const iec_t* c) const;
    virtual RVec ExchangeAngularIntegrals(const iec_t* a,const iec_t* b) const;
};

} //namespace

#endif //_Atoml_IE_HF_Angular_H_