// File: Imp/CompositeCD.C  Composite charged density, which is any array of Irrep DM_CDs.
module;
#include <cassert>
#include <vector>
#include <memory>
module qchem.CompositeCD;
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
SMatrix<double> Composite_CD::GetRepulsion(const Orbital_HF_IBS<double>* bs_ab) const
{
    int n=bs_ab->GetNumFunctions();
    SMatrix<double> J(n,n);
    Fill(J,0.0);
    for (auto& c:itsCDs) J+=c->GetRepulsion(bs_ab);
    return J;
}

SMatrix<double> Composite_CD::GetExchange(const Orbital_HF_IBS<double>* bs_ab) const
{
    int n=bs_ab->GetNumFunctions();
    SMatrix<double> K(n,n);
    Fill(K,0.0);
    for (auto& c:itsCDs) K+=c->GetExchange(bs_ab);
    return K;
}

double Composite_CD::DM_Contract(const Static_CC* v) const
{
    double ret=0.0;
    for (auto& c:itsCDs) ret+=c->DM_Contract(v);
    return ret;
}

double Composite_CD::DM_Contract(const Dynamic_CC* v,const DM_CD* cd) const
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
    Vector<double> ret(fbs->GetNumFunctions());
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
        ret += c->GetChangeFrom(**b);
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

RVec3 Composite_CD::Gradient  (const RVec3& r) const
{
    // No UT coverage
    RVec3 ret(0,0,0);
    for (auto& c:itsCDs) ret+=c->Gradient(r);
    return ret;
}


