// File: CompositeCD.C  Exact implementation of the charged density.



#include "ChargeDensityImplementation/CompositeCD/CompositeCD.H"
#include "BasisSet/BasisSet.H"
#include "Misc/ptrvector_io.h"
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

typedef optr_vector<ChargeDensity*>::iterator ITER;
typedef optr_vector<ChargeDensity*>::const_iterator CITER;
//-----------------------------------------------------------------------------
//
//  Totale energy terms for a charge density.
//
ChargeDensity::SMat CompositeCD::GetOverlap  (const BasisSet* bs) const
{
    int n=bs->GetNumFunctions();
    SMat S(n,n);
    Fill(S,0.0);
    for (CITER b(itsCDs.begin()); b!=itsCDs.end(); b++) S+=b->GetOverlap(bs);
    return S;
}

ChargeDensity::SMat CompositeCD::GetRepulsion(const BasisSet* bs_ab) const
{
    int n=bs_ab->GetNumFunctions();
    SMat J(n,n);
    Fill(J,0.0);
    for (CITER c(itsCDs.begin()); c!=itsCDs.end(); c++)
        J+=c->GetRepulsion(bs_ab);
    return J;
}

ChargeDensity::SMat CompositeCD::GetExchange(const BasisSet* bs_ab) const
{
    int n=bs_ab->GetNumFunctions();
    SMat K(n,n);
    Fill(K,0.0);
    for (CITER b(itsCDs.begin()); b!=itsCDs.end(); b++)
        K+=b->GetExchange(bs_ab);
    return K;
}

double CompositeCD::GetEnergy(const HamiltonianTerm* v) const
{
    double ret=0.0;
    for (CITER b(itsCDs.begin()); b!=itsCDs.end(); b++) ret+=b->GetEnergy(v);
    return ret;
}

double CompositeCD::GetTotalCharge() const
{
    double ret=0.0;
    for (CITER b(itsCDs.begin()); b!=itsCDs.end(); b++) ret+=b->GetTotalCharge();
    return ret;
}

//------------------------------------------------------------------------------
//
//  Required by fitting routines.
//
void CompositeCD::InjectOverlaps(FittedFunction* ff, const BasisSet* fbs) const
{
    for (CITER b(itsCDs.begin()); b!=itsCDs.end(); b++) b->InjectOverlaps(ff,fbs);
}

void CompositeCD::InjectRepulsions(FittedFunction* ff, const BasisSet* fbs) const
{
    for (CITER b(itsCDs.begin()); b!=itsCDs.end(); b++) b->InjectRepulsions(ff,fbs);
}

//-------------------------------------------------------------------------
//
//  SCF convergence stuff.
//
void CompositeCD::ReScale(double factor)
{
    for (ITER i(itsCDs.begin()); i!=itsCDs.end(); i++) i->ReScale(factor);
}

void CompositeCD::ShiftOrigin(const RVec3& newCenter)
{
    for (ITER i(itsCDs.begin()); i!=itsCDs.end(); i++) i->ShiftOrigin(newCenter);
}

void CompositeCD::MixIn(const ChargeDensity& cd,double c)
{
    const CompositeCD* ecd = dynamic_cast<const CompositeCD*>(&cd);
    assert(ecd);
    ITER i(itsCDs.begin());
    CITER  b(ecd->itsCDs.begin());
    for (; i!=itsCDs.end()&&b!=ecd->itsCDs.end(); i++,b++) i->MixIn(*b,c);
}

double CompositeCD::GetChangeFrom(const ChargeDensity& cd) const
{
    const CompositeCD* ecd = dynamic_cast<const CompositeCD*>(&cd);
    assert(ecd);
    CITER i(itsCDs.begin());
    CITER b(ecd->itsCDs.begin());
    double ret=0;
    for (; i!=itsCDs.end()&&b!=ecd->itsCDs.end(); i++,b++)
    {
        double t = i->GetChangeFrom(*b);
        ret = t > ret ? t : ret;
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
    for (CITER b(itsCDs.begin()); b!=itsCDs.end(); b++) ret+=b->operator()(r);
    return ret;
}

void  CompositeCD::Eval(const Mesh& mesh, Vec& v) const
{
    for (CITER b(itsCDs.begin()); b!=itsCDs.end(); b++) b->Eval(mesh,v); //Each  does v+= operation.
}

ChargeDensity::Vec3 CompositeCD::Gradient  (const RVec3& r) const
{
    Vec3 ret(0,0,0);
    for (CITER b(itsCDs.begin()); b!=itsCDs.end(); b++) ret+=b->Gradient(r);
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


