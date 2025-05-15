// File: Atom_IE_HF_Angular.H  Angular 2e-Integrals for atom HF basis sets.
#ifndef _Atom_IE_HF_Angular_H_
#define _Atom_IE_HF_Angular_H_

#include "oml/vector.h"

class AtomIE_BS_2E_Angular
{
public:
    typedef Vector<double> RVec;
    virtual RVec Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const=0;
    virtual RVec ExchangeAngularIntegrals(size_t la, size_t lc, int ma, int mc) const=0;
};

#endif //_Atom_IE_HF_Angular_H_
