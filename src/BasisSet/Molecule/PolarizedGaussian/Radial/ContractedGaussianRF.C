// File: ContractedGaussianRF.C  Contracted Gaussian in 3D space.


#include "Common/stl_io.h"
#include <iostream>

#include "PolarizedGaussian/Radial/ContractedGaussianH3.H"
#include "PolarizedGaussian/Radial/GaussianRF.H"
#include "PolarizedGaussian/Block.H"
#include "PolarizedGaussian/MnD/Hermite1.H"
#include "PolarizedGaussian/Radial/ContractedGaussianRF.H"

namespace PolarizedGaussian
{
//#######################################################################
//
//  Contracted Gaussian implementation
//
ContractedGaussianRF::ContractedGaussianRF()
    : RadialCommon()
    , cs          ()
    , gs          ()
{};

ContractedGaussianRF::ContractedGaussianRF(const Vector<double>& coeffs, std::vector<RadialFunction*>& rfs)
    : RadialCommon( rfs[0]->GetCenter (), rfs[0]->GetL())
    , cs            (coeffs)
    , unormalized_cs(coeffs)
    , gs() 
{
    for (auto r:rfs) gs.push_back(std::unique_ptr<RadialFunction>(r));
    assert(cs.size()==gs.size());
//
//  Absorb normalizations into the contraction coefficients.
//
    size_t i=1;
    for (auto& g:gs) cs(i++)*=g->GetNormalization(Polarization(GetL(),0,0));
    Check(); //Make sure all L's and centres are the same.
};

Hermite1* ContractedGaussianRF::MakeH1() const
{
    Hermite1* ret=new Hermite1;
    size_t i=1;
    for (auto& g:gs) ret->Add(g->GetH1(),cs(i++));
    return ret;
}

bool ContractedGaussianRF::operator==(const RadialFunction& rf) const
{
    bool ret=false;
    const ContractedGaussianRF* cg = dynamic_cast<const ContractedGaussianRF*>(&rf);
    if (cg)
    {
        bool base     = RadialCommon::operator==(rf);
        bool coeff    = cs.size()==cg->cs.size() && cs==cg->cs;
        bool index    = true;

        if (gs.size()==cg->gs.size())
            for (size_t i=0;i<gs.size();i++) 
                index = index && gs[i]->GetID()==cg->gs[i]->GetID();
        ret = base && coeff && index;
    }
    return ret;
}

double ContractedGaussianRF::GetNormalization(const Polarization& p) const
{
    // Client code should use the integral engine for this.
    assert(false);
    return 0;
}

double ContractedGaussianRF::GetCharge(const Polarization& p) const
{
    double ret=0.0;
    size_t i=1;
    for (auto& g:gs) ret += cs(i++)*g->GetCharge(p);
    return ret;
}

RadialFunction::sd_t ContractedGaussianRF::GetExponents() const
{
    sd_t ret;
    for (auto& i:gs)
        for (auto& e:i->GetExponents()) ret.insert(e);
    return ret;
}

RadialFunction::vd_t ContractedGaussianRF::GetCoeff() const
{
    vd_t ret;
//    for (auto c:unormalized_cs) ret.push_back(c);
    for (auto c:cs) ret.push_back(c);
    return ret;
}
 

double ContractedGaussianRF::Integrate(qchem::IType2C type,const RadialFunction* rb, const Polarization& pa, const Polarization& pb,CDCache& cache,const Cluster* cl) const
{
    double s=0;
    size_t i=1;
    for (auto& g:gs) s += cs(i++)*rb->Integrate(type,g.get(),pb,pa,cache,cl); //swap pols
    return s;
}

//
//  At his we know this is radial c.  Keep c at the end of the call
//
double ContractedGaussianRF::Integrate(qchem::IType3C type,const RadialFunction* ra, const RadialFunction* rb, const Polarization& pa, const Polarization& pb, const Polarization& pc,CDCache& cache) const
{
    double s=0;
    size_t i=1;
    for (auto& g:gs) 
        s += cs(i++)*ra->Integrate(type,rb,pb,pa,pc,cache,g.get()); //swap pols
    return s;
}

double ContractedGaussianRF::Integrate(qchem::IType3C type,const RadialFunction* ra, const Polarization& pa, const Polarization& pb, const Polarization& pc,CDCache& cache,const RadialFunction* rc) const
{
    double s=0;
    size_t i=1;
    for (auto& g:gs)
        s += cs(i++)*ra->Integrate(type,g.get(),pb,pa,pc,cache,rc); //swap pols
    return s;
}

// this is rd
double ContractedGaussianRF::Integrate(rf_t* ra,rf_t* rb,rf_t* rc,po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache& cache) const
{
    double s=0;
    size_t i=1;
    for (auto& g:gs)
        s += cs(i++)*rc->Integrate(ra,rb,pa,pb,pc,pd,cache,g.get()); //swap pols
    return s;
}

// this is rc
double ContractedGaussianRF::Integrate(rf_t* ra,rf_t* rb,po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache& cache, rf_t* rd) const
{
    double s=0;
    size_t i=1;
    for (auto& g:gs)
        s += cs(i++)*rb->Integrate(ra,pa,pb,pc,pd,cache,g.get(),rd); //swap pols
    return s;
}

// this is rb
double ContractedGaussianRF::Integrate(rf_t* ra,po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache& cache, rf_t* rc, rf_t* rd) const
{
    double s=0;
    size_t i=1;
    for (auto& g:gs)
        s += cs(i++)*ra->Integrate(g.get(),pb,pa,pc,pd,cache,rc,rd); //swap pols
    return s;
}



Hermite3* ContractedGaussianRF::GetH3(const RadialFunction& r1, const RadialFunction& r2) const
{
    ContractedGaussianH3* ret = new ContractedGaussianH3(cs);
    for (auto& g:gs) ret->Insert( g->GetH3(r1,r2) );
    return ret;
}




std::ostream& ContractedGaussianRF::Write(std::ostream& os) const
{
    os << "Contracted {";
    for(auto& g:gs)
    {
        std::ostringstream s;
        s << *g;
        std::string ss(s.str());
        std::string num=ss.substr(10);
        os << num << " ";
    }
    os << "}";
    return os;
}


double ContractedGaussianRF::operator()(const RVec3& r) const
{
    double ret=0;
    size_t i=1;
    for (auto& g:gs) ret+=cs(i++)* g->operator()(r);
    return ret;
}

RVec3 ContractedGaussianRF::Gradient(const RVec3& r) const
{
    RVec3 ret(0,0,0);
    size_t i=1;
    for (auto& g:gs) ret+=cs(i++)* g->Gradient(r);
    return ret;
}

void ContractedGaussianRF::Eval(const Mesh& mesh, Vector<double>& vec) const
{
    size_t ig=1;
    for (auto& g:gs)
    {
        Vector<double>::iterator i(vec.begin());
        Vector<double>::const_iterator  v(g->operator()(mesh).begin());
        for (; i!=vec.end()&&v; i++,v++) *i += cs(ig++) * (*v);
    }
}

void ContractedGaussianRF::EvalGrad(const Mesh& mesh, Vector<RVec3>& vec) const
{
//    CVITER c(cs.begin());
//    CITER p(gs.begin());
//    for (; p!=gs.end()&&c!=cs.end(); c++,p++)
    size_t ig=1;
    for (auto& g:gs)
    {
        Vector<RVec3>::iterator i(vec.begin());
        Vector<RVec3>::const_iterator  v(g->Gradient(mesh).begin());
        for (; i&&v; i++,v++) *i += cs(ig++)* (*v);
    }
}


RadialFunction* ContractedGaussianRF::Clone() const
{
    assert(false);
    return 0;
    // return new  ContractedGaussianRF(*this);
}



RadialFunction* ContractedGaussianRF::Clone(const RVec3& newCenter) const
{
    std::vector<RadialFunction*> newList;
    for (auto& g:gs) newList.push_back(g->Clone(newCenter));
    return new  ContractedGaussianRF(cs,newList);
}

//
//  Make sure all primatives are self consistent.
//
void ContractedGaussianRF::Check() const
{
    const RVec3& center=gs[0]->GetCenter();
    int          L     =gs[0]->GetL     ();
    for(unsigned i=1; i<gs.size(); i++)
    {
        if (center!=gs[i]->GetCenter() || L!=gs[i]->GetL())
        {
            std::cerr
                << "ContractedGaussianRF ID=" << GetID()
                << "Primatives have different L or Center" << std::endl
                << "center(0)=" << center << "  L(0)=" << L << std::endl
                << "center(" << i << ")=" << gs[i]->GetCenter()
                << "  L(" << i << ")=" << gs[i]->GetL() << std::endl;
        }
    }
}

} //namespace PolarizedGaussian
