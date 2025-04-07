// File: FittedVxc.C  Fitted exchange potential.

  

#include "Imp/Hamiltonian/FittedVxc.H"
#include "Imp/Hamiltonian/LDAVxc.H"
#include <DFT_IBS.H>
#include <TotalEnergy.H>
#include <ChargeDensity.H>
#include "oml/smatrix.h"
#include "oml/imp/binio.h"

FittedVxc::FittedVxc()
    : FittedFunctionImp<double>( )
    , itsLDAVxc                 (0)
{};


FittedVxc::FittedVxc(bs_t& bs, ex_t& lda,mesh_t& m)
    : FittedFunctionImp<double>(bs,m) //Use regular overlap for fitting.
    , itsLDAVxc                (new LDAVxc(lda))
{};

FittedVxc::~FittedVxc()
{
    delete itsLDAVxc;
}

void FittedVxc::UseChargeDensity(const DM_CD* cd)
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

Static_HT::SMat FittedVxc::CalculateHamiltonianMatrix(const ibs_t* bs,const Spin&,const DM_CD* cd) const
{
    if (newCD(cd))
    {
        itsLDAVxc->UseChargeDensity(cd);
        FittedVxc* cfvxc=const_cast<FittedVxc*>(this);
        cfvxc->DoFit(*itsLDAVxc); //use the callback GetFunctionOverlap
    }
    auto dftbs=dynamic_cast<const TOrbital_DFT_IBS<double>*>(bs);
    return FitGet3CenterOverlap(dftbs);
}

void FittedVxc::GetEnergy(TotalEnergy& te,const DM_CD* cd) const
{
    if (newCD(cd))
    {
        itsLDAVxc->UseChargeDensity(cd);
        FittedVxc* cfvxc=const_cast<FittedVxc*>(this);
        cfvxc->DoFit(*itsLDAVxc); //use the callback GetFunctionOverlap
    }
    te.Exc += 3.0/4.0 *cd->GetEnergy(this,cd);

//    double HFExc=-0.25*itsCD->GetExchangeEnergy();
//    std::cout.precision(4);
//    std::cout.width(7);
//    std::cout.setf(std::ios::fixed,std::ios::floatfield);
//    std::cout << "  Exc DFT,HF=" << te.Exc << "," << HFExc << std::endl;
}

std::ostream& FittedVxc::Write(std::ostream& os) const
{
    FittedFunctionImp<double>::Write(os);
    os << itsLDAVxc;
    return os;
}

std::istream& FittedVxc::Read (std::istream& is)
{
    FittedFunctionImp<double>::Read(is);
    delete itsLDAVxc;
    itsLDAVxc = FittablePotential::Factory(is);
    is >> itsLDAVxc;
    return is;
}

