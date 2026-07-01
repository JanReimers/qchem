// File: Fitted.C  Fitted polarized exchange potential.
module;
#include <cassert>
#include <iostream>
#include <memory>
#include <vector>
module qchem.Hamiltonian.Internal.Terms;
import qchem.ChargeDensity;
import qchem.Energy;
import qchem.Symmetry.Spin;

namespace qchem::Hamiltonian
{


FittedVxcPol::FittedVxcPol(fbs_t& bs, ex_t& lda)
    : itsUpVxc               (new FittedVxc(bs,lda))
    , itsDownVxc             (new FittedVxc(bs,lda))
{
    assert(itsUpVxc);
    assert(itsDownVxc);
};

FittedVxcPol::~FittedVxcPol()
{
    delete itsUpVxc;
    delete itsDownVxc;
}


//########################################################################
//
//  This is where we calculate the overlap of the fit basis functions with
//  the real exchange potential,  Vxc(ro(r)), where ro is the charge density.
//
// The Hamiltonain matrix elements are calculated
//             /
//  Vxc(i,j) = | dr Vxcfit(ro(r)) Oi(r) Oj(r) .
//             /
//
//           = Sum  { Ck <Oi|Vk|Oj> } .
//
//  This last part is carried out by the base class FitImplementation.
rsmat_t FittedVxcPol::CalcMatrix(const obs_t* bs,const Spin& s,const rChargeDensity* cd) const
{
    assert(itsUpVxc);
    assert(itsDownVxc);

    if  (s==Spin::None)
    {
        std::cerr << "PolarizedFittedVxc::GetMatrix Asking for unpolarized result in Polarized Vxc" << std::endl;
        exit(-1);
    }
    const Polarized_CD* pol_cd =  dynamic_cast<const Polarized_CD*>(cd);
    if (!pol_cd)
    {
        // Spin-unpolarized SEED density (e.g. the SAD total rho, implicitly rho_up = rho_down = rho/2).
        // For Slater exchange Vx^sigma(rho/2) == Vx_unpolarized(rho_total), so each spin channel's
        // iteration-0 Vxc is just the unpolarized Vxc of the total seed density.  (Only the seed is
        // spin-agnostic; once the SCF builds orbitals the density is a Polarized_CD and the per-spin
        // branch below runs.)  Previously the dynamic_cast yielded null and -- with the assert compiled
        // out in Release -- segfaulted on the SAD-seeded polarized DFT path.
        return (s==Spin::Up ? itsUpVxc : itsDownVxc)->GetMatrix(bs, s, cd);
    }

    const DM_CD* ucd = pol_cd->GetChargeDensity(Spin::Up  );
    const DM_CD* dcd = pol_cd->GetChargeDensity(Spin::Down);

    rsmat_t Kab= s==Spin::Up ? itsUpVxc  ->GetMatrix(bs,s,ucd) : itsDownVxc->GetMatrix(bs,s,dcd);
    return Kab;
}


void FittedVxcPol::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    assert(itsUpVxc);
    assert(itsDownVxc);
    const Polarized_CD* pol_cd =  dynamic_cast<const Polarized_CD*>(cd);
    assert(pol_cd);

    const DM_CD* ucd = pol_cd->GetChargeDensity(Spin::Up  );
    const DM_CD* dcd = pol_cd->GetChargeDensity(Spin::Down);
    te.Exc = 0.0;
    itsUpVxc  ->GetEnergy(te,ucd);
    itsDownVxc->GetEnergy(te,dcd);
}

std::ostream& FittedVxcPol::Write(std::ostream& os) const
{
    assert(itsUpVxc);
    assert(itsDownVxc);

    return os << itsUpVxc << itsDownVxc;
}

} //namespace
