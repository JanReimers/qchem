// File: Atom_ml_IE_HF_Angular.H  Angular 2e-Integrals for atoml HF basis sets.
#ifndef _Atom_ml_IE_HF_Angular_H_
#define _Atom_ml_IE_HF_Angular_H_

import qchem.BasisSet.Atom.IE;

namespace Atom_ml
{

class IE_BS_2E_Angular : public virtual ::AtomIE_BS_2E_Angular
{
public:
    virtual RVec Coulomb_AngularIntegrals(const iec_t* a,const iec_t* c) const;
    virtual RVec ExchangeAngularIntegrals(const iec_t* a,const iec_t* b) const;
};

} //namespace

#endif //_Atom_ml_IE_HF_Angular_H_