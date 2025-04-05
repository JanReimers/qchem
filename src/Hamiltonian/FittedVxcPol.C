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
    : itsUpVxc                    (0)
    , itsDownVxc                  (0)
{};

FittedVxcPol::FittedVxcPol(bs_t& bs, ex_t& lda, mesh_t& m)
    : itsUpVxc               (new FittedVxc(bs,lda,m))
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

void FittedVxcPol::UseChargeDensity(const Exact_CD* cd)
{
    assert(itsUpVxc);
    assert(itsDownVxc);

    Dynamic_HT_Imp::UseChargeDensity(cd);

    const Polarized_CD* pol_cd =  dynamic_cast<const Polarized_CD*>(cd);
    assert(pol_cd);

    const Exact_CD* ucd = pol_cd->GetChargeDensity(Spin::Up  );
    const Exact_CD* dcd = pol_cd->GetChargeDensity(Spin::Down);
    itsUpVxc  ->UseChargeDensity(ucd  );
    itsDownVxc->UseChargeDensity(dcd);
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
Static_HT::SMat FittedVxcPol::CalculateHamiltonianMatrix(const TOrbital_IBS<double>* bs,const Spin& s) const
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
   
    return is;
}

