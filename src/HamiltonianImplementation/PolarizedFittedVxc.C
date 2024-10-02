// File: PolarizedFittedVxc.C  Fitted exchange potential.



#include "HamiltonianImplementation/PolarizedFittedVxc.H"
#include "HamiltonianImplementation/FittedVxc.H"
#include "ChargeDensity.H"
#include "TotalEnergy.H"
#include "Misc/Spin.H"
#include "oml/smatrix.h"
#include "oml/vector3d.h"
#include <cassert>
#include <iostream>
#include <stdlib.h>

PolarizedFittedVxc::PolarizedFittedVxc()
    : HamiltonianTermImplementation     ( )
    , itsUpVxc                    (0)
    , itsDownVxc                  (0)
{};

PolarizedFittedVxc::PolarizedFittedVxc(const rc_ptr<BasisSet>& bs, const rc_ptr<ExchangeFunctional>& lda)
    : HamiltonianTermImplementation(   )
    , itsUpVxc               (new FittedVxc(bs,lda))
    , itsDownVxc             (new FittedVxc(bs,lda))
{
    assert(itsUpVxc);
    assert(itsDownVxc);
};

PolarizedFittedVxc::~PolarizedFittedVxc()
{
    delete itsUpVxc;
    delete itsDownVxc;
}

void PolarizedFittedVxc::UseChargeDensity(const ChargeDensity* exactCD)
{
    assert(itsUpVxc);
    assert(itsDownVxc);

    HamiltonianTermImplementation::UseChargeDensity(exactCD);

    const PolarizedCD* PolExactCD =  dynamic_cast<const PolarizedCD*>(exactCD);
    assert(PolExactCD);

    const ChargeDensity* upExactCD   = PolExactCD->GetChargeDensity(Spin::Up  );
    const ChargeDensity* downExactCD = PolExactCD->GetChargeDensity(Spin::Down);
    itsUpVxc  ->UseChargeDensity(upExactCD  );
    itsDownVxc->UseChargeDensity(downExactCD);
}

bool PolarizedFittedVxc::IsPolarized() const
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
HamiltonianTerm::SMat PolarizedFittedVxc::CalculateHamiltonianMatrix(const BasisSet* bs,const Spin& s) const
{
    assert(itsUpVxc);
    assert(itsDownVxc);

    if  (s.itsState==Spin::None)
    {
        std::cerr << "PolarizedFittedVxc::GetMatrix Asking for unpolarized result in Polarized Vxc" << std::endl;
        exit(-1);
    }
    SMat Kab= s.itsState==Spin::Up ? itsUpVxc  ->BuildHamiltonian(bs,s) : itsDownVxc->BuildHamiltonian(bs,s);
    return Kab;
}


void PolarizedFittedVxc::GetEnergy(TotalEnergy& te) const
{
    assert(itsUpVxc);
    assert(itsDownVxc);

    te.Exc = 0.0;
    itsUpVxc  ->GetEnergy(te);
    itsDownVxc->GetEnergy(te);
}

std::ostream& PolarizedFittedVxc::Write(std::ostream& os) const
{
    assert(itsUpVxc);
    assert(itsDownVxc);

    return os << itsUpVxc << itsDownVxc;
}

std::istream& PolarizedFittedVxc::Read (std::istream& is)
{
    delete itsUpVxc;
    itsUpVxc = HamiltonianTerm::Factory(is);
    is >> *itsUpVxc;

    delete itsDownVxc;
    itsDownVxc = HamiltonianTerm::Factory(is);
    is >> *itsDownVxc;

    return is;
}

