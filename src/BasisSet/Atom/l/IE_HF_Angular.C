// File: Atoml_IE_HF_Angular.H  Angular 2e-Integrals for atoml HF basis sets.

#include <cmath>
#include "AngularIntegrals.H"
#include "l/IE_HF_Angular.H"
import qchem.BasisSet.Atom.IEClient;

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