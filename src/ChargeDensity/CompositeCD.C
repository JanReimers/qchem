// File: CompositeCD.C  Exact implementation of the charged density.



#include "Imp/ChargeDensity/CompositeCD.H"
#include <HF_IBS.H>
#include <Fit_IBS.H>
#include "oml/smatrix.h"
#include "oml/vector.h"
#include <cassert>

//------------------------------------------------------------------------------------
//
//  Construction zone.
//
Composite_CD::Composite_CD()
{};

void Composite_CD::Insert(DM_CD* cd)
{
    itsCDs.push_back(std::unique_ptr<DM_CD>(cd));
}

//-----------------------------------------------------------------------------
//
//  Total energy terms for a charge density.
//
DM_CD::SMat Composite_CD::GetRepulsion(const TOrbital_HF_IBS<double>* bs_ab) const
{
    int n=bs_ab->GetNumFunctions();
    SMat J(n,n);
    Fill(J,0.0);
    for (auto& c:itsCDs) J+=c->GetRepulsion(bs_ab);
    return J;
}

DM_CD::SMat Composite_CD::GetExchange(const TOrbital_HF_IBS<double>* bs_ab) const
{
    int n=bs_ab->GetNumFunctions();
    SMat K(n,n);
    Fill(K,0.0);
    for (auto& c:itsCDs) K+=c->GetExchange(bs_ab);
    return K;
}

double Composite_CD::DM_Contract(const Static_HT* v) const
{
    double ret=0.0;
    for (auto& c:itsCDs) ret+=c->DM_Contract(v);
    return ret;
}

double Composite_CD::DM_Contract(const Dynamic_HT* v,const DM_CD* cd) const
{
    double ret=0.0;
    for (auto& c:itsCDs) ret+=c->DM_Contract(v,cd);
    return ret;
}

double Composite_CD::GetTotalCharge() const
{
    double ret=0.0;
    for (auto& c:itsCDs) ret+=c->GetTotalCharge();
    return ret;
}

//------------------------------------------------------------------------------
//
//  Required by fitting routines.
//
Vector<double> Composite_CD::GetRepulsion3C(const Fit_IBS* fbs) const
{
    Vector<double> ret(fbs->size());
    Fill(ret,0.0);
    for (auto& c:itsCDs) ret+=c->GetRepulsion3C(fbs);
    return ret;
}

//-------------------------------------------------------------------------
//
//  SCF convergence stuff.
//
void Composite_CD::ReScale(double factor)
{
    // No UT coverage
    for (auto& c:itsCDs) c->ReScale(factor);
}

void Composite_CD::ShiftOrigin(const RVec3& newCenter)
{
    // No UT coverage
    for (auto& c:itsCDs) c->ShiftOrigin(newCenter);
}

void Composite_CD::MixIn(const DM_CD& cd,double f)
{
    const Composite_CD* ecd = dynamic_cast<const Composite_CD*>(&cd);
    assert(ecd);
    auto  b(ecd->itsCDs.begin());
    for (auto& c:itsCDs)
    {
        c->MixIn(**b,f);
        b++;
    }
}

double Composite_CD::GetChangeFrom(const DM_CD& cd) const
{
    const Composite_CD* ecd = dynamic_cast<const Composite_CD*>(&cd);
    assert(ecd);
    assert(itsCDs.size()==ecd->itsCDs.size());
    auto  b(ecd->itsCDs.begin());
    double ret=0;
    for (auto& c:itsCDs)
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
double Composite_CD::operator()(const RVec3& r) const
{
    double ret=0.0;
    for (auto& c:itsCDs) ret+=c->operator()(r);
    return ret;
}

DM_CD::Vec3 Composite_CD::Gradient  (const RVec3& r) const
{
    // No UT coverage
    Vec3 ret(0,0,0);
    for (auto& c:itsCDs) ret+=c->Gradient(r);
    return ret;
}


