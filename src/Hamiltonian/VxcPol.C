// File: VxcPol.C  Polarized HF exchange potential

#include <cassert>
#include <iostream>
#include <memory>
#include <vector>


#include "VxcPol.H"
#include "Vxc.H"
import qchem.HF_IBS;
#include <Hamiltonian/TotalEnergy.H>
#include <ChargeDensity/ChargeDensity.H>
import qchem.Symmetry.Spin;

VxcPol::VxcPol()
{
};

VxcPol::~VxcPol()
{
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
Static_HT::SMat VxcPol::CalcMatrix(const ibs_t* bs,const Spin& s,const DM_CD* cd) const
{
    if  (s==Spin::None)
    {
        std::cerr << "PolarizedHartreeFockVxc::GetMatrix Asking for unpolarized result in Polarized Vxc" << std::endl;
        exit(-1);
    }
    auto hf_bs = dynamic_cast<const TOrbital_HF_IBS<double>*>(bs);
    assert(hf_bs);

    const Polarized_CD* PolExactCD =  dynamic_cast<const Polarized_CD*>(cd);
    assert(PolExactCD);
    const DM_CD* SpinCD   = PolExactCD->GetChargeDensity(s); //Get CD for this spin direction
    SMat Kab=SpinCD->GetExchange(hf_bs)*-1.0;
    return Kab;
}
void VxcPol::GetEnergy(EnergyBreakdown& te,const DM_CD* cd) const
{
    te.Exc += 0.5*cd->DM_Contract(this,cd); //This should sum K^alpha and K^beta.
}

std::ostream& VxcPol::Write(std::ostream& os) const
{
    os << "    Polarized Hartee-Fock exchange potential." << std::endl;
    return os;
}


