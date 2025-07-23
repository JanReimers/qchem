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
Static_HT::SMat FittedVxcPol::CalcMatrix(const ibs_t* bs,const Spin& s,const DM_CD* cd) const
{
    assert(itsUpVxc);
    assert(itsDownVxc);

    if  (s==Spin::None)
    {
        std::cerr << "PolarizedFittedVxc::GetMatrix Asking for unpolarized result in Polarized Vxc" << std::endl;
        exit(-1);
    }
    const Polarized_CD* pol_cd =  dynamic_cast<const Polarized_CD*>(cd);
    assert(pol_cd);

    const DM_CD* ucd = pol_cd->GetChargeDensity(Spin::Up  );
    const DM_CD* dcd = pol_cd->GetChargeDensity(Spin::Down);

    SMatrix<double> Kab= s==Spin::Up ? itsUpVxc  ->GetMatrix(bs,s,ucd) : itsDownVxc->GetMatrix(bs,s,dcd);
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

std::istream& FittedVxcPol::Read (std::istream& is)
{
   
    return is;
}

