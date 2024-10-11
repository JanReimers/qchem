// File: FittedVxc.C  Fitted exchange potential.



#include "HamiltonianImplementation/FittedVxc.H"
#include "HamiltonianImplementation/LDAVxc.H"
#include "TotalEnergy.H"
#include "ChargeDensity.H"
#include "oml/smatrix.h"
#include "oml/imp/binio.h"

FittedVxc::FittedVxc()
    : HamiltonianTermImplementation     ( )
    , FittedFunctionImplementation<double>( )
    , itsLDAVxc                 (0)
{};


FittedVxc::FittedVxc(const rc_ptr<IrrepBasisSet>& bs, const rc_ptr<ExchangeFunctional>& lda, Mesh* m)
    : HamiltonianTermImplementation     (   )
    , FittedFunctionImplementation<double>(bs,m,false) //Use regular overlap for fitting.
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

HamiltonianTerm::SMat FittedVxc::CalculateHamiltonianMatrix(const IrrepBasisSet* bs,const Spin&) const
{
    const FittedFunctionImplementation<double>* ffi=dynamic_cast<const FittedFunctionImplementation<double>*>(this);
    assert(ffi);
    const std::vector<SMat>& overlap=bs->GetOverlap3C(ffi->itsBasisSet.get());
    int n=bs->GetNumFunctions();
    SMat Kab(n,n);
    Fill(Kab,0.0);
    size_t i=0;
    for (auto c:ffi->itsFitCoeff) Kab+=SMat(c*overlap[i++]);
    assert(!isnan(Kab));
    
    //SMat Kab=bs->GetOverlap(this);
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

