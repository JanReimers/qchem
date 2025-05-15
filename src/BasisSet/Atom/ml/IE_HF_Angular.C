// File: Atom_ml_IE_HF_Angular.H  Angular 2e-Integrals for atom-ml HF basis sets.

#include "Imp/BasisSet/Atom/ml/IE_HF_Angular.H"
#include "Imp/Integrals/AngularIntegrals.H"

namespace Atom_ml
{

IE_BS_2E_Angular::RVec IE_BS_2E_Angular::Coulomb_AngularIntegrals(size_t la, size_t lc, int ma, int mc) const
{
    return AngularIntegrals::Coulomb(la,lc,ma,mc);
}
IE_BS_2E_Angular::RVec IE_BS_2E_Angular::ExchangeAngularIntegrals(size_t la, size_t lb, int ma, int mb) const
{
    return AngularIntegrals::Exchange(la,lb,ma,mb);
}

}

