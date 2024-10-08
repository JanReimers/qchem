// File: ContractedGaussianRF.C  Contracted Gaussian in 3D space.


#include "Misc/Polarization.H"
#include "BasisSetImplementation/PolarizedGaussian/Gaussian/GaussianRF.H"
#include "BasisSetImplementation/PolarizedGaussian/ContractedGaussian/ContractedGaussianRF.H"
#include "BasisSetImplementation/PolarizedGaussian/ContractedGaussian/ContractedGaussianH3.H"
#include "BasisSetImplementation/PolarizedGaussian/BasisFunctionBlock.H"
#include "BasisSetImplementation/PolarizedGaussian/Hermite/Hermite1.H"
#include "Mesh/MeshBrowser.H"
#include "Misc/ptr_vector1_io.h"
#include "oml/imp/binio.h"
#include "oml/io3d.h"
#include "oml/smatrix.h"
#include "oml/matrix.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>

//#######################################################################
//
//  Contracted Gaussian implementation
//
ContractedGaussianRF::ContractedGaussianRF()
    : RadialFunctionImplementation()
    , cs             ( )
    , gs         ( )
{};

ContractedGaussianRF::ContractedGaussianRF(const Vector<double>& coeffs, std::vector<RadialFunction*>& rfs)
    : RadialFunctionImplementation( rfs[0]->GetCenter (), rfs[0]->GetL())
    , cs    (coeffs)
    , gs(rfs   ) //Copy pointers
{
    assert(cs.size()==gs.size());
//
//  Absorb normalizations into the contraction coefficients.
//
    for (auto i:gs.indices()) cs(i+1)*=gs[i]->GetNormalization(Polarization(GetL(),0,0));
    Check(); //Make sure all L's and centres are the same.
};


Hermite1* ContractedGaussianRF::MakeH1() const
{
    Hermite1* ret=new Hermite1;
    for (auto i:gs.indices()) ret->Add(gs[i]->GetH1(),cs(i+1));
    return ret;
}

bool ContractedGaussianRF::operator==(const RadialFunction& rf) const
{
    bool ret=false;
    const ContractedGaussianRF* cg = dynamic_cast<const ContractedGaussianRF*>(&rf);
    if (cg)
    {
        bool base     = RadialFunctionImplementation::operator==(rf);
        bool coeff    = cs.size()==cg->cs.size() && cs==cg->cs;
        bool index    = true;

        for (auto i:gs.indices()) index = index && gs[i]->GetID()==cg->gs[i]->GetID();
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
    for (auto i:gs.indices()) ret += cs(i+1)*gs[i]->GetCharge(p);
    return ret;
}

double ContractedGaussianRF::Integrate(Types2C type,const RadialFunction* rb, const Polarization& pa, const Polarization& pb,CDcache& cache,const Cluster* cl) const
{
    double s=0;
    for (auto i:gs.indices()) s += cs(i+1)*rb->Integrate(type,gs[i],pb,pa,cache,cl); //swap pols
    return s;
}

//
//  At his we know this is radial c.  Keep c at the end of the call
//
double ContractedGaussianRF::Integrate(Types3C type,const RadialFunction* ra, const RadialFunction* rb, const Polarization& pa, const Polarization& pb, const Polarization& pc,CDcache& cache) const
{
    double s=0;
    for (auto i:gs.indices()) 
        s += cs(i+1)*ra->Integrate(type,rb,pb,pa,pc,cache,gs[i]); //swap pols
    return s;
}

double ContractedGaussianRF::Integrate(Types3C type,const RadialFunction* ra, const Polarization& pa, const Polarization& pb, const Polarization& pc,CDcache& cache,const RadialFunction* rc) const
{
    double s=0;
    for (auto i:gs.indices()) 
        s += cs(i+1)*ra->Integrate(type,gs[i],pb,pa,pc,cache,rc); //swap pols
    return s;
}

// this is rd
double ContractedGaussianRF::Integrate(rf_t* ra,rf_t* rb,rf_t* rc,po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDcache& cache) const
{
    double s=0;
    for (auto i:gs.indices()) 
        s += cs(i+1)*rc->Integrate(ra,rb,pa,pb,pc,pd,cache,gs[i]); //swap pols
    return s;
}

// this is rc
double ContractedGaussianRF::Integrate(rf_t* ra,rf_t* rb,po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDcache& cache, rf_t* rd) const
{
    double s=0;
    for (auto i:gs.indices()) 
        s += cs(i+1)*rb->Integrate(ra,pa,pb,pc,pd,cache,gs[i],rd); //swap pols
    return s;
}

// this is rb
double ContractedGaussianRF::Integrate(rf_t* ra,po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDcache& cache, rf_t* rc, rf_t* rd) const
{
    double s=0;
    for (auto i:gs.indices()) 
        s += cs(i+1)*ra->Integrate(gs[i],pb,pa,pc,pd,cache,rc,rd); //swap pols
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
    if (Binary()) os << cs << gs;
    if (Ascii ()) os << cs << std::endl << gs << std::endl;
    if (Pretty())
    {
        os << "Contracted {";
        for(auto& g:gs)
        {
            std::ostringstream s;
            s << g << std::ends;
            std::string ss(s.str());
            std::string num=ss.substr(10);
            os << num << " ";
        }
        os << "}";
    }
    if (!Pretty()) RadialFunctionImplementation::Write(os);
    return os;
}

std::istream& ContractedGaussianRF::Read(std::istream& is)
{
    is >> cs >> gs;
    RadialFunctionImplementation::Read(is);
    Check(); //Make all L's and centers are the same.
    return is;
}

double ContractedGaussianRF::operator()(const RVec3& r) const
{
    double ret=0;
    for (auto i:gs.indices()) ret+=cs(i+1)* gs[i]->operator()(r);
    return ret;
}

RVec3 ContractedGaussianRF::Gradient(const RVec3& r) const
{
    RVec3 ret(0,0,0);
    for (auto i:gs.indices()) ret+=cs(i+1)* gs[i]->Gradient(r);
    return ret;
}

void ContractedGaussianRF::Eval(const Mesh& mesh, Vector<double>& vec) const
{
    for (auto ig:gs.indices())
    {
        Vector<double>::iterator i(vec.begin());
        Vector<double>::const_iterator  v(gs[ig]->operator()(mesh).begin());
        for (; i!=vec.end()&&v; i++,v++) *i += cs(ig+1) * (*v);
    }
}

void ContractedGaussianRF::EvalGrad(const Mesh& mesh, Vector<RVec3>& vec) const
{
//    CVITER c(cs.begin());
//    CITER p(gs.begin());
//    for (; p!=gs.end()&&c!=cs.end(); c++,p++)
    for (auto ig:gs.indices())
    {
        Vector<RVec3>::iterator i(vec.begin());
        Vector<RVec3>::const_iterator  v(gs[ig]->Gradient(mesh).begin());
        for (; i&&v; i++,v++) *i += cs(ig+1)* (*v);
    }
}


RadialFunction* ContractedGaussianRF::Clone() const
{
    return new  ContractedGaussianRF(*this);
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




