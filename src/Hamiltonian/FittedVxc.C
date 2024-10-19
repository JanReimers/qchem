// File: FittedVxc.C  Fitted exchange potential.

 

#include "Imp/Hamiltonian/FittedVxc.H"
#include "Imp/Hamiltonian/LDAVxc.H"
#include <TotalEnergy.H>
#include <ChargeDensity.H>
#include "oml/smatrix.h"
#include "oml/imp/binio.h"

FittedVxc::FittedVxc()
    : HamiltonianTermImp     ( )
    , FittedFunctionImplementation<double>( )
    , itsLDAVxc                 (0)
{};


FittedVxc::FittedVxc(const rc_ptr<IrrepBasisSet>& bs, const rc_ptr<ExFunctional>& lda,const rc_ptr<Mesh>& m)
    : HamiltonianTermImp     (   )
    , FittedFunctionImplementation<double>(bs,m) //Use regular overlap for fitting.
    , itsLDAVxc                 (new LDAVxc(lda))
{};

FittedVxc::~FittedVxc()
{
    delete itsLDAVxc;
}

void FittedVxc::UseChargeDensity(const ChargeDensity* exactCD)
{
    HamiltonianTermImp::UseChargeDensity(exactCD);
    itsLDAVxc->UseChargeDensity(exactCD);

    DoFit(*itsLDAVxc); //use the callback GetFunctionOverlap
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

HamiltonianTerm::SMat FittedVxc::CalculateHamiltonianMatrix(const IrrepBasisSet* bs,const Spin&) const
{
    return FitGet3CenterOverlap(bs);
}

void FittedVxc::GetEnergy(TotalEnergy& te) const
{
    te.Exc += 3.0/4.0 *CalculateEnergy();

//    double HFExc=-0.25*itsExactCD->GetExchangeEnergy();
//    std::cout.precision(4);
//    std::cout.width(7);
//    std::cout.setf(std::ios::fixed,std::ios::floatfield);
//    std::cout << "  Exc DFT,HF=" << te.Exc << "," << HFExc << std::endl;
}

std::ostream& FittedVxc::Write(std::ostream& os) const
{
    FittedFunctionImplementation<double>::Write(os);
    os << itsLDAVxc;
    return os;
}

std::istream& FittedVxc::Read (std::istream& is)
{
    FittedFunctionImplementation<double>::Read(is);
    delete itsLDAVxc;
    itsLDAVxc = FittablePotential::Factory(is);
    is >> itsLDAVxc;
    return is;
}

