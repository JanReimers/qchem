// File: Fitted.C  Fitted polarized exchange potential.



#include "Imp/Hamiltonian/FittedVxc.H"
#include "Imp/Hamiltonian/FittedVxcPol.H"
#include <ChargeDensity.H>
#include <TotalEnergy.H>
#include <Spin.H>
#include "oml/smatrix.h"
#include "oml/vector3d.h"
#include <cassert>
#include <iostream>
#include <stdlib.h>

FittedVxcPol::FittedVxcPol()
    : HamiltonianTermImp     ( )
    , itsUpVxc                    (0)
    , itsDownVxc                  (0)
{};

FittedVxcPol::FittedVxcPol(bs_t& bs, ex_t& lda, mesh_t& m)
    : HamiltonianTermImp(   )
    , itsUpVxc               (new FittedVxc(bs,lda,m))
    , itsDownVxc             (new FittedVxc(bs,lda,m))
{
    assert(itsUpVxc);
    assert(itsDownVxc);
};

FittedVxcPol::~FittedVxcPol()
{
    delete itsUpVxc;
    delete itsDownVxc;
}

void FittedVxcPol::UseChargeDensity(const ChargeDensity* exactCD)
{
    assert(itsUpVxc);
    assert(itsDownVxc);

    HamiltonianTermImp::UseChargeDensity(exactCD);

    const Polarized_CD* PolExactCD =  dynamic_cast<const Polarized_CD*>(exactCD);
    assert(PolExactCD);

    const ChargeDensity* upExactCD   = PolExactCD->GetChargeDensity(Spin::Up  );
    const ChargeDensity* downExactCD = PolExactCD->GetChargeDensity(Spin::Down);
    itsUpVxc  ->UseChargeDensity(upExactCD  );
    itsDownVxc->UseChargeDensity(downExactCD);
}

bool FittedVxcPol::IsPolarized() const
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
HamiltonianTerm::SMat FittedVxcPol::CalculateHamiltonianMatrix(const TOrbital_IBS<double>* bs,const Spin& s) const
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


void FittedVxcPol::GetEnergy(TotalEnergy& te) const
{
    assert(itsUpVxc);
    assert(itsDownVxc);

    te.Exc = 0.0;
    itsUpVxc  ->GetEnergy(te);
    itsDownVxc->GetEnergy(te);
}

std::ostream& FittedVxcPol::Write(std::ostream& os) const
{
    assert(itsUpVxc);
    assert(itsDownVxc);

    return os << itsUpVxc << itsDownVxc;
}

std::istream& FittedVxcPol::Read (std::istream& is)
{
    delete itsUpVxc;
    itsUpVxc = HamiltonianTerm::Factory(is);
    is >> *itsUpVxc;

    delete itsDownVxc;
    itsDownVxc = HamiltonianTerm::Factory(is);
    is >> *itsDownVxc;

    return is;
}

