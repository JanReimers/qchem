// File: PolarizedHartreeFockVxc.C  Fitted exchange potential.



#include "HamiltonianImplementation/PolarizedHartreeFockVxc.H"
#include "HamiltonianImplementation/HartreeFockVxc.H"
#include "TotalEnergy.H"
#include "ChargeDensity.H"
#include "BasisSet.H"
#include "Misc/Spin.H"
#include "oml/smatrix.h"
#include "oml/vector3d.h"
#include <cassert>
#include <iostream>
#include <stdlib.h>

PolarizedHartreeFockVxc::PolarizedHartreeFockVxc()
    : HamiltonianTermImplementation     ( )
{
};

PolarizedHartreeFockVxc::~PolarizedHartreeFockVxc()
{
}

bool PolarizedHartreeFockVxc::IsPolarized() const
{
    return true;
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
HamiltonianTerm::SMat PolarizedHartreeFockVxc::CalculateHamiltonianMatrix(const IrrepBasisSet* bs,const Spin& s) const
{
    if  (s.itsState==Spin::None)
    {
        std::cerr << "PolarizedHartreeFockVxc::GetMatrix Asking for unpolarized result in Polarized Vxc" << std::endl;
        exit(-1);
    }
    const PolarizedCD* PolExactCD =  dynamic_cast<const PolarizedCD*>(itsExactCD);
    assert(PolExactCD);
    const ChargeDensity* SpinCD   = PolExactCD->GetChargeDensity(s); //Get CD for this spin direction
    SMat Kab=SpinCD->GetExchange(bs)*-1.0;
    return Kab;
}
void PolarizedHartreeFockVxc::GetEnergy(TotalEnergy& te) const
{
    te.Exc += 0.5*CalculateEnergy(); //This should sum K^alpha and K^beta.
}

std::ostream& PolarizedHartreeFockVxc::Write(std::ostream& os) const
{
    if (StreamableObject::Pretty())
        os << "    Polarized Hartee-Fock exchange potential." << std::endl;
    return os;
}


