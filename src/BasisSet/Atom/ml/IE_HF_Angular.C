// File: Atom_ml_IE_HF_Angular.H  Angular 2e-Integrals for atom-ml HF basis sets.

#include "Imp/BasisSet/Atom/ml/IE_HF_Angular.H"
#include "Imp/Integrals/AngularIntegrals.H"
#include "Imp/BasisSet/Atom/IEC.H"

namespace Atom_ml
{

IE_BS_2E_Angular::RVec IE_BS_2E_Angular::Coulomb_AngularIntegrals(const iec_t* a,const iec_t* c) const
 {
    return AngularIntegrals::Coulomb(a->l,c->l,a->m,c->m);
 }
IE_BS_2E_Angular::RVec IE_BS_2E_Angular::ExchangeAngularIntegrals(const iec_t* a,const iec_t* b) const
{
    return AngularIntegrals::Exchange(a->l,b->l,a->m,b->m);
}

}

