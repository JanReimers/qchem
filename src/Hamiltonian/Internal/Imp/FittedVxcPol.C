// File: Fitted.C  Fitted polarized exchange potential.
module;
#include <cassert>
#include <iostream>
#include <memory>
#include <stdexcept>
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
// R2.19: a pure FORWARDER -- return the CHILD'S cached reference.  It used to return by value, copying a
// matrix the child's own per-Irrep cache already held, which the tDynamic_HT_Imp_NoCache base then stored
// in scratch purely to have something to hand a reference to: one full matrix copy per call, per spin, per
// irrep, per SCF iteration, to satisfy a signature.  VxcPol -- the Hartree-Fock twin of this class -- has
// always done it this way (Imp/VxcPol.C).
const rsmat_t& FittedVxcPol::GetMatrix(const robs_t* bs,const Spin& s,const rChargeDensity* cd) const
{
    assert(itsUpVxc);
    assert(itsDownVxc);

    // A polarized term has no Spin::None block to hand back -- and the caller has to be able to SEE that
    // (R2.5: exit(-1) killed the pybind GUI and the test runner outright, taking the diagnostic with it).
    if (s==Spin::None)
        throw std::runtime_error("FittedVxcPol::GetMatrix: asked for the Spin::None (unpolarized) block of "
                                 "a polarized Vxc term -- a polarized term has an Up and a Down block, no "
                                 "total.");
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

    const rDM_CD* ucd = pol_cd->GetChargeDensity(Spin::Up  );
    const rDM_CD* dcd = pol_cd->GetChargeDensity(Spin::Down);

    return s==Spin::Up ? itsUpVxc->GetMatrix(bs,s,ucd) : itsDownVxc->GetMatrix(bs,s,dcd);
}


void FittedVxcPol::GetEnergy(EnergyBreakdown& te,const rDM_CD* cd) const
{
    assert(itsUpVxc);
    assert(itsDownVxc);
    const Polarized_CD* pol_cd =  dynamic_cast<const Polarized_CD*>(cd);
    assert(pol_cd);

    const rDM_CD* ucd = pol_cd->GetChargeDensity(Spin::Up  );
    const rDM_CD* dcd = pol_cd->GetChargeDensity(Spin::Down);
    itsUpVxc  ->GetEnergy(te,ucd);
    itsDownVxc->GetEnergy(te,dcd);
}

std::ostream& FittedVxcPol::Write(std::ostream& os) const
{
    assert(itsUpVxc);
    assert(itsDownVxc);

    return os << *itsUpVxc << *itsDownVxc;
}

} //namespace
