// File: Atoml_IE_HF_Angular.H  Angular 2e-Integrals for atoml HF basis sets.

#include "Atom/l/IE_HF_Angular.H"
#include "Atom/AngularIntegrals.H"
#include "Atom/IEC.H"

namespace Atoml
{

 IE_BS_2E_Angular::RVec IE_BS_2E_Angular::Coulomb_AngularIntegrals(const iec_t* a,const iec_t* c) const
 {
    return AngularIntegrals::Coulomb(a->l,c->l);
 }
IE_BS_2E_Angular::RVec IE_BS_2E_Angular::ExchangeAngularIntegrals(const iec_t* a,const iec_t* b) const
{
    return AngularIntegrals::Exchange(a->l,b->l);
}

}