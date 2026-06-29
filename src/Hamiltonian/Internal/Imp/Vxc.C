// File: Vxc.C  Hartree-Fock exchange potential
module;
#include <iostream>
#include <cassert>
#include <memory>
#include <vector>
module qchem.Hamiltonian.Internal.Terms;
import qchem.Hamiltonian.Types;
import qchem.ChargeDensity;
import qchem.Energy;
import qchem.Blaze;

namespace qchem::Hamiltonian
{

//########################################################################
//
//  Let the charge density do the work.
//

rsmat_t Vxc::CalcMatrix(const obs_t* bs,const Spin&,const rChargeDensity* cd) const
{
    newCD(cd); //Set H matrix cache to dirty if cd really is new.
    auto hf_bs = dynamic_cast<const ohfbs_t*>(bs);
    assert(hf_bs);
    const DM_CD* dm = dynamic_cast<const DM_CD*>(cd);   // HF K needs the density matrix (not a fit seed)
    assert(dm && "Vxc (HF exchange): density must be a DM_CD");
    rsmat_t Kab=blazem::zero<double>(bs->GetNumFunctions());
    dm->AccumulateExchange(Kab,hf_bs);
    return Kab*-0.5;
}
void Vxc::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    newCD(cd); //Set H matrix cache to dirty if cd really is new.
    te.Exc+=0.5*cd->DM_Contract(this,cd);
}

std::ostream& Vxc::Write(std::ostream& os) const
{
    os << "    Hartee-Fock exchange potential phi(r_1)*phi(r_2)/r_12" << std::endl;
    return os;
}

} //namespace
