// File: Atoml_IE_HF_Angular.H  Angular 2e-Integrals for atoml HF basis sets.

#include "Imp/BasisSet/Atom/l/IE_HF_Angular.H"
#include "Imp/Integrals/AngularIntegrals.H"

namespace Atoml
{
IE_BS_2E_Angular::RVec IE_BS_2E_Angular::Coulomb_AngularIntegrals(size_t la, size_t lc, int, int) const
{
    return AngularIntegrals::Coulomb(la,lc);
}
IE_BS_2E_Angular::RVec IE_BS_2E_Angular::ExchangeAngularIntegrals(size_t la, size_t lb, int, int) const
{
    return AngularIntegrals::Exchange(la,lb);
}

}