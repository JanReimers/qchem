// File: CompositeCD.C  Exact implementation of the charged density.



#include "Imp/ChargeDensity/CompositeCD.H"
#include <Irrep_BS.H>
#include "Imp/Containers/ptr_vector_io.h"
#include "oml/smatrix.h"
#include "oml/vector.h"
#include <cassert>

//------------------------------------------------------------------------------------
//
//  Construction zone.
//
Composite_Exact_CD::Composite_Exact_CD()
{};

void Composite_Exact_CD::Insert(Exact_CD* cd)
{
    itsCDs.push_back(cd);
}

typedef optr_vector1<Exact_CD*>::iterator ITER;
typedef optr_vector1<Exact_CD*>::const_iterator CITER;
//-----------------------------------------------------------------------------
//
//  Total energy terms for a charge density.
//
Exact_CD::SMat Composite_Exact_CD::GetRepulsion(const TOrbital_IBS<double>* bs_ab) const
{
    int n=bs_ab->GetNumFunctions();
    SMat J(n,n);
    Fill(J,0.0);
    for (auto c:itsCDs) J+=c->GetRepulsion(bs_ab);
    return J;
}

Exact_CD::SMat Composite_Exact_CD::GetExchange(const TOrbital_IBS<double>* bs_ab) const
{
    int n=bs_ab->GetNumFunctions();
    SMat K(n,n);
    Fill(K,0.0);
    for (auto c:itsCDs) K+=c->GetExchange(bs_ab);
    return K;
}

double Composite_Exact_CD::GetEnergy(const HamiltonianTerm* v) const
{
    double ret=0.0;
    for (auto c:itsCDs) ret+=c->GetEnergy(v);
    return ret;
}

double Composite_Exact_CD::GetTotalCharge() const
{
    double ret=0.0;
    for (auto c:itsCDs) ret+=c->GetTotalCharge();
    return ret;
}

//------------------------------------------------------------------------------
//
//  Required by fitting routines.
//
Vector<double> Composite_Exact_CD::GetRepulsion3C(const Fit_IBS* fbs) const
{
    Vector<double> ret(fbs->size());
    Fill(ret,0.0);
    for (auto c:itsCDs) ret+=c->GetRepulsion3C(fbs);
    return ret;
}

//-------------------------------------------------------------------------
//
//  SCF convergence stuff.
//
void Composite_Exact_CD::ReScale(double factor)
{
    // No UT coverage
    for (auto c:itsCDs) c->ReScale(factor);
}

void Composite_Exact_CD::ShiftOrigin(const RVec3& newCenter)
{
    // No UT coverage
    for (auto c:itsCDs) c->ShiftOrigin(newCenter);
}

void Composite_Exact_CD::MixIn(const Exact_CD& cd,double f)
{
    const Composite_Exact_CD* ecd = dynamic_cast<const Composite_Exact_CD*>(&cd);
    assert(ecd);
    CITER  b(ecd->itsCDs.begin());
    for (auto c:itsCDs)
    {
        c->MixIn(**b,f);
        b++;
    }
}

double Composite_Exact_CD::GetChangeFrom(const Exact_CD& cd) const
{
    const Composite_Exact_CD* ecd = dynamic_cast<const Composite_Exact_CD*>(&cd);
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
double Composite_Exact_CD::operator()(const RVec3& r) const
{
    double ret=0.0;
    for (auto c:itsCDs) ret+=c->operator()(r);
    return ret;
}

Exact_CD::Vec3 Composite_Exact_CD::Gradient  (const RVec3& r) const
{
    // No UT coverage
    Vec3 ret(0,0,0);
    for (auto c:itsCDs) ret+=c->Gradient(r);
    return ret;
}


