// File: FittedVxc.C  Fitted exchange potential.



#include "HamiltonianImplementation/FittedVxc.H"
#include "HamiltonianImplementation/LDAVxc.H"
#include "Hamiltonian/TotalEnergy.H"
#include "ChargeDensity/ChargeDensity.H"
#include "oml/smatrix.h"
#include "oml/imp/binio.h"

FittedVxc::FittedVxc()
    : HamiltonianTermImplementation     ( )
    , FittedFunctionImplementation<double>( )
    , itsLDAVxc                 (0)
{};


FittedVxc::FittedVxc(const rc_ptr<BasisSet>& bs, const rc_ptr<ExchangeFunctional>& lda)
    : HamiltonianTermImplementation     (   )
    , FittedFunctionImplementation<double>(bs,false)
    , itsLDAVxc                 (new LDAVxc(lda))
{};

FittedVxc::~FittedVxc()
{
    delete itsLDAVxc;
}

void FittedVxc::UseChargeDensity(const ChargeDensity* exactCD)
{
    HamiltonianTermImplementation::UseChargeDensity(exactCD);
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
//#include "BasisSet/QuantumNumber.H"

HamiltonianTerm::SMat FittedVxc::CalculateHamiltonianMatrix(const BasisSet* bs,const Spin&) const
{

    SMat Kab=bs->GetOverlap(this);
//    std::cout.precision(4);
//    std::cout.width(7);
//    std::cout.setf(std::ios::fixed,std::ios::floatfield);
//    std::cout << "  Basis Set QN=" << bs->GetQuantumNumber() << std::endl;
//    std::cout << "  Fitted Kab=:" << Kab << std::endl;
//    SMat Kabhf=itsExactCD->GetExchange(bs)*-0.5;
//    std::cout << "  HF Kab=:" << Kabhf << std::endl;
//    SMat delta=DirectDivide((Kabhf-Kab),Kab)*100;
//    std::cout << "  % Delta Kab=:" << delta << std::endl;
    return Kab;
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

