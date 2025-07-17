// File: Atom_IE_HF_Angular.H  Angular 2e-Integrals for atom HF basis sets.
#ifndef _Atom_IE_HF_Angular_H_
#define _Atom_IE_HF_Angular_H_

import oml;

class AtomIrrepIEClient;

class AtomIE_BS_2E_Angular
{
public:
    typedef Vector<double> RVec;
    typedef AtomIrrepIEClient iec_t;
    virtual RVec Coulomb_AngularIntegrals(const iec_t* a,const iec_t* c) const=0;
    virtual RVec ExchangeAngularIntegrals(const iec_t* a,const iec_t* c) const=0;
};

#endif //_Atom_IE_HF_Angular_H_
