// File: CompositeCD.C  Exact implementation of the charged density.



#include "ChargeDensityImplementation/CompositeCD/CompositeCD.H"
#include "BasisSet.H"
#include "Imp/Containers/ptr_vector_io.h"
#include "oml/smatrix.h"
#include <cassert>

//------------------------------------------------------------------------------------
//
//  Construction zone.
//
CompositeCD::CompositeCD()
{};

void CompositeCD::Insert(ChargeDensity* cd)
{
    itsCDs.push_back(cd);
}

typedef optr_vector1<ChargeDensity*>::iterator ITER;
typedef optr_vector1<ChargeDensity*>::const_iterator CITER;
//-----------------------------------------------------------------------------
//
//  Totale energy terms for a charge density.
//
ChargeDensity::SMat CompositeCD::GetOverlap  (const IrrepBasisSet* bs) const
{
    // No UT coverage
    int n=bs->GetNumFunctions();
    SMat S(n,n);
    Fill(S,0.0);
    for (auto c:itsCDs) S+=c->GetOverlap(bs);
    return S;
}

ChargeDensity::SMat CompositeCD::GetRepulsion(const IrrepBasisSet* bs_ab) const
{
    int n=bs_ab->GetNumFunctions();
    SMat J(n,n);
    Fill(J,0.0);
    for (auto c:itsCDs) J+=c->GetRepulsion(bs_ab);
    return J;
}

ChargeDensity::SMat CompositeCD::GetExchange(const IrrepBasisSet* bs_ab) const
{
    int n=bs_ab->GetNumFunctions();
    SMat K(n,n);
    Fill(K,0.0);
    for (auto c:itsCDs) K+=c->GetExchange(bs_ab);
    return K;
}

double CompositeCD::GetEnergy(const HamiltonianTerm* v) const
{
    double ret=0.0;
    for (auto c:itsCDs) ret+=c->GetEnergy(v);
    return ret;
}

double CompositeCD::GetTotalCharge() const
{
    double ret=0.0;
    for (auto c:itsCDs) ret+=c->GetTotalCharge();
    return ret;
}

//------------------------------------------------------------------------------
//
//  Required by fitting routines.
//
void CompositeCD::InjectOverlaps(FittedFunction* ff, const IrrepBasisSet* fbs) const
{
    // No UT coverage
    for (auto c:itsCDs) c->InjectOverlaps(ff,fbs);
}

void CompositeCD::InjectRepulsions(FittedFunction* ff, const IrrepBasisSet* fbs) const
{
    for (auto c:itsCDs) c->InjectRepulsions(ff,fbs);
}

//-------------------------------------------------------------------------
//
//  SCF convergence stuff.
//
void CompositeCD::ReScale(double factor)
{
    // No UT coverage
    for (auto c:itsCDs) c->ReScale(factor);
}

void CompositeCD::ShiftOrigin(const RVec3& newCenter)
{
    // No UT coverage
    for (auto c:itsCDs) c->ShiftOrigin(newCenter);
}

void CompositeCD::MixIn(const ChargeDensity& cd,double f)
{
    const CompositeCD* ecd = dynamic_cast<const CompositeCD*>(&cd);
    assert(ecd);
    CITER  b(ecd->itsCDs.begin());
    for (auto c:itsCDs)
    {
        c->MixIn(**b,f);
        b++;
    }
}

double CompositeCD::GetChangeFrom(const ChargeDensity& cd) const
{
    const CompositeCD* ecd = dynamic_cast<const CompositeCD*>(&cd);
    assert(ecd);
    assert(itsCDs.size()==ecd->itsCDs.size());
    CITER b(ecd->itsCDs.begin());
    double ret=0;
    for (auto c:itsCDs)
    {
        double t = c->GetChangeFrom(**b);
        ret = t > ret ? t : ret;
        b++;
    }
    return ret;
}

//-------------------------------------------------------------------------
//
//  Real space function stuff.
//
double CompositeCD::operator()(const RVec3& r) const
{
    double ret=0.0;
    for (auto c:itsCDs) ret+=c->operator()(r);
    return ret;
}

void  CompositeCD::Eval(const Mesh& mesh, Vec& v) const
{
    // No UT coverage
    for (auto c:itsCDs) c->Eval(mesh,v); //Each  does v+= operation.
}

ChargeDensity::Vec3 CompositeCD::Gradient  (const RVec3& r) const
{
    // No UT coverage
    Vec3 ret(0,0,0);
    for (auto c:itsCDs) ret+=c->Gradient(r);
    return ret;
}


//-----------------------------------------------------------------------
//
//  Streamable stuff.
//
std::ostream& CompositeCD::Write(std::ostream& os) const
{
    os << itsCDs;
    return os;
}

std::istream& CompositeCD::Read (std::istream& is)
{
    is >> itsCDs;
    return is;
}


