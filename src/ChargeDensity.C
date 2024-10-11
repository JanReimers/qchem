// File: ChargeDensity.C  Interface for the charge density category.
#include "ChargeDensity.H"

double ChargeDensity::FitGetConstraint  () const
{
    return  GetTotalCharge();
}

bool ChargeDensity::IsPolarized() const
{
    return false;
}

#include "oml/smatrix.h"
#include "Misc/Spin.H"
#include "oml/vector.h"
//----------------------------------------------------------------------------
//
//  Various integrals.
//
ChargeDensity::SMat PolarizedCD::GetOverlap  (const IrrepBasisSet* bs) const
{
    return
        GetChargeDensity(Spin::Up  )->GetOverlap(bs) +
        GetChargeDensity(Spin::Down)->GetOverlap(bs);
}

ChargeDensity::SMat PolarizedCD::GetRepulsion(const IrrepBasisSet* bs) const
{
//    std::cout.precision(4);
//    std::cout.width(7);
//    std::cout.setf(std::ios::fixed,std::ios::floatfield);
//    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
//    std::cout << "Spin up:" << std::endl;
    SMat Jab_up=GetChargeDensity(Spin::Up  )->GetRepulsion(bs);
//    std::cout << "  Jab_up=:" << Jab_up << std::endl;

//    std::cout << "Spin down:" << std::endl;
    SMat Jab_down=GetChargeDensity(Spin::Down)->GetRepulsion(bs);
//    std::cout << "  Jab_down=:" << Jab_down << std::endl;
//    std::cout << "Jab=:" << Jab_up + Jab_down << std::endl;
    return Jab_up + Jab_down;
}

ChargeDensity::SMat PolarizedCD::GetExchange(const IrrepBasisSet* bs) const
{
//    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
//    std::cout << "Spin up:" << std::endl;
    SMat Kab_up=GetChargeDensity(Spin::Up  )->GetExchange(bs);
//    std::cout << "  Kab_up=:" << Kab_up << std::endl;
//
//    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
//    std::cout << "Spin down:" << std::endl;
    SMat Kab_down=GetChargeDensity(Spin::Down)->GetExchange(bs);
//    std::cout << "  Kab_down=:" << Kab_up << std::endl;
    return Kab_up + Kab_down;

//    return
//        GetChargeDensity(Spin::Up  )->GetExchange(bs) +
//        GetChargeDensity(Spin::Down)->GetExchange(bs);
}

double PolarizedCD::GetEnergy(const HamiltonianTerm* v) const
{
    return GetChargeDensity(Spin::Up  )->GetEnergy(v)+GetChargeDensity(Spin::Down)->GetEnergy(v);
}

double PolarizedCD::GetTotalCharge() const
{
    return GetChargeDensity(Spin::Up)->GetTotalCharge() + GetChargeDensity(Spin::Down)->GetTotalCharge() ;
}

double PolarizedCD::GetTotalSpin() const
{
    return GetChargeDensity(Spin::Up)->GetTotalCharge() - GetChargeDensity(Spin::Down)->GetTotalCharge() ;
}

Vector<double> PolarizedCD::GetRepulsions(const IrrepBasisSet* theFitBasisSet) const
{
    return GetChargeDensity(Spin::Up  )->GetRepulsions(theFitBasisSet)
        +  GetChargeDensity(Spin::Down)->GetRepulsions(theFitBasisSet);
}

//-----------------------------------------------------------------------
//
//  Convergence and origin shifting.
//
void   PolarizedCD::ShiftOrigin(const RVec3& newcenter)
{
    GetChargeDensity(Spin::Up)  ->ShiftOrigin(newcenter) ;
    GetChargeDensity(Spin::Down)->ShiftOrigin(newcenter) ;
}

void PolarizedCD::MixIn(const ChargeDensity& cd,double c)
{
    const PolarizedCD* pcd = dynamic_cast<const PolarizedCD*>(&cd);
    if (!pcd)
    {
        std::cerr << "PolarizedCD::MixIn could not cast cd" << std::endl;
        exit(-1);
    }
    GetChargeDensity(Spin::Up)  -> MixIn(*pcd->GetChargeDensity(Spin::Up  ),c);
    GetChargeDensity(Spin::Down)-> MixIn(*pcd->GetChargeDensity(Spin::Down),c);
}

double PolarizedCD::GetChangeFrom(const ChargeDensity& cd) const
{
    const PolarizedCD* pcd = dynamic_cast<const PolarizedCD*>(&cd);
    if (!pcd)
    {
        std::cerr << "PolarizedCD::GetChangeFrom could not cast cd" << std::endl;
        exit(-1);
    }
    return GetChargeDensity(Spin::Up)  ->GetChangeFrom(*pcd->GetChargeDensity(Spin::Up  ))
           + GetChargeDensity(Spin::Down)->GetChangeFrom(*pcd->GetChargeDensity(Spin::Down)) ;
}

void PolarizedCD::ReScale(double factor)
{
    GetChargeDensity(Spin::Up)  ->ReScale(factor);
    GetChargeDensity(Spin::Down)->ReScale(factor);
}


bool PolarizedCD::IsPolarized() const
{
    return true;
}

//----------------------------------------------------------------------------------
//
//  Real space function stuff.
//
double PolarizedCD::operator()(const RVec3& r) const
{
    return (*GetChargeDensity(Spin::Up))(r) + (*GetChargeDensity(Spin::Down))(r);
}

ChargeDensity::RVec3 PolarizedCD::Gradient  (const RVec3& r) const
{
    return GetChargeDensity(Spin::Up)->Gradient(r) + GetChargeDensity(Spin::Down)->Gradient(r);
}

void PolarizedCD::Eval(const Mesh& m, Vec& v) const
{
    GetChargeDensity(Spin::Up)  ->Eval(m,v);
    GetChargeDensity(Spin::Down)->Eval(m,v);
}

//---------------------------------------------------------------------------------
//
//  Construction zone.
//
FittedPolarizedCD::FittedPolarizedCD()
    : itsSpinUpCD  (0)
    , itsSpinDownCD(0)
{};

FittedPolarizedCD::FittedPolarizedCD(const ChargeDensity* unpolcd, double Stotal)
    : itsSpinUpCD  (0)
    , itsSpinDownCD(0)
{
    const FittedCD* fcd=dynamic_cast<const FittedCD*>(unpolcd);
    if (!fcd)
    {
        std::cerr << "FittedPolarizedCD::FittedPolarizedCD could cast input charge density" << std::endl;
        exit(-1);
    }
    itsSpinUpCD  =fcd->Clone();
    itsSpinDownCD=fcd->Clone();

    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    double totalCharge=unpolcd->GetTotalCharge();
    itsSpinUpCD  ->ReScale((totalCharge+Stotal)/(2*totalCharge));
    itsSpinDownCD->ReScale((totalCharge-Stotal)/(2*totalCharge));
};

FittedPolarizedCD::FittedPolarizedCD(ChargeDensity* up, ChargeDensity* down)
    : itsSpinUpCD  (dynamic_cast<FittedCD*>(up  ))
    , itsSpinDownCD(dynamic_cast<FittedCD*>(down))
{
    if (!itsSpinUpCD)
    {
        std::cerr << "FittedPolarizedCD::FittedPolarizedCD could cast up charge density" << std::endl;
        exit(-1);
    }
    if (!itsSpinDownCD)
    {
        std::cerr << "FittedPolarizedCD::FittedPolarizedCD could cast down charge density" << std::endl;
        exit(-1);
    }
};

FittedPolarizedCD::FittedPolarizedCD(const FittedPolarizedCD& pcd)
    : itsSpinUpCD  (pcd.itsSpinUpCD  ->Clone())
    , itsSpinDownCD(pcd.itsSpinDownCD->Clone())
{
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
};

FittedPolarizedCD::~FittedPolarizedCD()
{
    delete itsSpinUpCD;
    delete itsSpinDownCD;
}


//-------------------------------------------------------------------
//
//  Access to individual components.
//
ChargeDensity* FittedPolarizedCD::GetChargeDensity(const Spin& S)
{
    assert(S.itsState!=Spin::None);
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    ChargeDensity* ret=0;
    if(S.itsState==Spin::Up  ) ret=itsSpinUpCD  ;
    if(S.itsState==Spin::Down) ret=itsSpinDownCD;
    return ret;
}

const ChargeDensity* FittedPolarizedCD::GetChargeDensity(const Spin& S) const
{
    assert(S.itsState!=Spin::None);
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    const ChargeDensity* ret=0;
    if(S.itsState==Spin::Up  ) ret=itsSpinUpCD  ;
    if(S.itsState==Spin::Down) ret=itsSpinDownCD;
    return ret;
}

double FittedPolarizedCD::GetSelfRepulsion() const
{
    return GetRepulsion(this);
}

// <ro(1) | 1/r12 | ff(2)>
double FittedPolarizedCD::GetRepulsion(const FittedFunction* ff) const
{
    double Jup  =itsSpinUpCD  ->GetRepulsion(ff);
    double Jdown=itsSpinDownCD->GetRepulsion(ff);

    return Jup + Jdown;
}

//-------------------------------------------------------------------------
//
//  Fitted function stuff.
//
double FittedPolarizedCD::DoFit(const DensityFFClient& ffc)
{
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    const PolarizedCD* polcd=dynamic_cast<const PolarizedCD*>(&ffc);
    assert(polcd);
    double lam_bar=0;
    lam_bar += itsSpinUpCD  ->DoFit(*polcd->GetChargeDensity(Spin::Up  ));
    lam_bar += itsSpinDownCD->DoFit(*polcd->GetChargeDensity(Spin::Down));
    return lam_bar/2.0;
}

double FittedPolarizedCD::DoFit(const ScalarFFClient& ffc)
{
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    double lam_bar=0;
    lam_bar += itsSpinUpCD  ->DoFit(ffc);
    lam_bar += itsSpinDownCD->DoFit(ffc);
    return lam_bar/2.0;
}

Vector<double> FittedPolarizedCD::GetRepulsions(const IrrepBasisSet* theFitBasisSet) const
{
    return GetChargeDensity(Spin::Up  )->GetRepulsions(theFitBasisSet)
        +  GetChargeDensity(Spin::Down)->GetRepulsions(theFitBasisSet);
    
}


void FittedPolarizedCD::ReScale(double factor) //Fit *= factor
{
    itsSpinUpCD   -> ReScale(factor);
    itsSpinDownCD -> ReScale(factor);
}

void FittedPolarizedCD::ShiftOrigin(const RVec3& newCenter)
{
    itsSpinUpCD   -> ShiftOrigin(newCenter);
    itsSpinDownCD -> ShiftOrigin(newCenter);
}

void FittedPolarizedCD::FitMixIn(const FittedFunction& ff,double c) // this = this*(1-c) + that*c.
{
    const FittedPolarizedCD* polcd=dynamic_cast<const FittedPolarizedCD*>(&ff);
    assert(polcd);
    itsSpinUpCD   -> FitMixIn(*polcd->itsSpinUpCD  ,c);
    itsSpinDownCD -> FitMixIn(*polcd->itsSpinDownCD,c);
}

double FittedPolarizedCD::FitGetChangeFrom(const FittedFunction& ff) const
{
    return itsSpinUpCD   -> FitGetChangeFrom(ff)
           + itsSpinDownCD -> FitGetChangeFrom(ff);
}

//--------------------------------------------------------------------------
//
//  Streamable stuff.
//
std::ostream& FittedPolarizedCD::Write(std::ostream& os) const
{
    assert(itsSpinUpCD);
    assert(itsSpinDownCD);
    return os << itsSpinUpCD << itsSpinDownCD;
}

std::istream& FittedPolarizedCD::Read (std::istream& is)
{
    delete itsSpinUpCD;
    itsSpinUpCD=FittedCD::Factory(is);
    assert(itsSpinUpCD);
    is >> *itsSpinUpCD;

    delete itsSpinDownCD;
    itsSpinDownCD=FittedCD::Factory(is);
    assert(itsSpinDownCD);
    is >> *itsSpinDownCD;

    return is;
}

FittedCD* FittedPolarizedCD::Clone() const
{
    return new FittedPolarizedCD(*this);
}


