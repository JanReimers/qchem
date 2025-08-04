// File: VxcPol.C  Polarized HF exchange potential
module;
#include <cassert>
#include <iostream>
#include <memory>
#include <vector>
module qchem.Hamiltonian.Internal.Terms;
import qchem.Orbital_HF_IBS;
import qchem.Energy;
import qchem.ChargeDensity;
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
 SMatrix<double>  VxcPol::CalcMatrix(const ibs_t* bs,const Spin& s,const DM_CD* cd) const
{
    if  (s==Spin::None)
    {
        std::cerr << "PolarizedHartreeFockVxc::GetMatrix Asking for unpolarized result in Polarized Vxc" << std::endl;
        exit(-1);
    }
    auto hf_bs = dynamic_cast<const Orbital_HF_IBS<double>*>(bs);
    assert(hf_bs);

    const Polarized_CD* PolExactCD =  dynamic_cast<const Polarized_CD*>(cd);
    assert(PolExactCD);
    const DM_CD* SpinCD   = PolExactCD->GetChargeDensity(s); //Get CD for this spin direction
    SMatrix<double> Kab=SpinCD->GetExchange(hf_bs)*-1.0;
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


