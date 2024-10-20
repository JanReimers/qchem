// File: GaussianRF.C  Primative Gaussian in 3D space.


#include "Imp/BasisSet/PolarizedGaussian/MnD/Hermite1.H"
#include "Imp/BasisSet/PolarizedGaussian/MnD/RNLM.H"
#include "Imp/BasisSet/PolarizedGaussian/Radial/GaussianRF.H"
#include "Imp/BasisSet/PolarizedGaussian/Radial/GaussianCD.H"
#include "Imp/BasisSet/PolarizedGaussian/Radial/GaussianH3.H"
#include "Imp/BasisSet/PolarizedGaussian/CDCache.H"
#include "Imp/BasisSet/PolarizedGaussian/Block.H"
#include <Cluster.H>
#include "oml/smatrix.h"
#include "oml/imp/binio.h"
#include "oml/io3d.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <stdlib.h>

namespace PolarizedGaussian
{
//#######################################################################
//
//   Gaussian radial function implementation
//

GaussianRF::GaussianRF()
    : RadialCommon()
    , itsExponent(0)
{};

GaussianRF::GaussianRF(double theExponent, const RVec3& theCenter, int theL)
    : RadialCommon(theCenter,theL)
    , itsExponent(theExponent)
{
    if(itsExponent < 0)
    {
        std::cerr << "GaussianRF exponent < 0" << std::endl;
        exit(-1);
    }
};

Hermite1* GaussianRF::MakeH1() const
{
    return new Hermite1(itsExponent,GetL());
}


bool GaussianRF::operator==(const RadialFunction& rf) const
{
    bool ret=false;
    const GaussianRF* g = dynamic_cast<const GaussianRF*>(&rf);
    if (g)
    {
        bool base     = RadialCommon::operator==(rf);
        bool exponent = fabs((itsExponent-g->itsExponent)/(itsExponent+g->itsExponent)*2.0) < 0.001; // 0.1%.
        ret= exponent && base;
    }
    return ret;
}

double GaussianRF::Integrate(qchem::IType2C type,const RadialFunction* rb, const Polarization& pa, const Polarization& pb,CDCache& cache,const Cluster* cl) const
{
    double s=0.0;
    const GaussianRF* gb=dynamic_cast<const GaussianRF*>(rb);
    if (!gb)
        return rb->Integrate(type,this,pb,pa,cache,cl); 

    Polarization zero(0,0,0);
    const GaussianCD& ab=cache.find(this,gb);
    switch (type)
    {
        case qchem::Overlap2C :
            s=pow(Pi/ab.AlphaP,1.5)*ab.Eij*ab.H2(zero,pa,pb);
            break;
        case qchem::Repulsion2C :
            {
                auto NLMs=GaussianCD::GetNMLs(this->GetL());
                const Hermite1& H1a=this->GetH1();
                const Hermite1& H1b=gb->GetH1();
                const RNLM& R=cache.find(ab);
                 
                double factor=1.0/(ab.ab*sqrt(ab.AlphaP));
                factor = (pb.GetTotalL()%2) ? -factor : factor;

                for (auto bNLM:NLMs)
                {
                    if (bNLM> pa) continue;
                    double ha=H1a(bNLM,pa);
                    if (ha==0.0) continue;
                    double RR=0.0;
                    for (int n=0; n<=pb.n; n++)
                        for (int l=0; l<=pb.l; l++)
                            for (int m=0; m<=pb.m; m++)
                            {
                                Polarization NLMp(n,l,m);
                                double hb=H1b(NLMp,pb);
                                if (hb!=0.0)
                                    RR+=hb*R(bNLM+NLMp);
                            } //for nlm
                    
                    s += ha*RR;
                }//for (auto bNLM:NLMs)
                s*=2*Pi52*factor;
            }  // case          
            break;
        case qchem::Kinetic :
            {
                double factor=0.5*pow(Pi/ab.AlphaP,1.5)*ab.Eij;
                double h = GetKinetic(pa,pb,ab);
                if (h!=0) s=factor*h;
                break;
            }
        case qchem::Nuclear :
            {
                assert(cl);
                RNLM R; //Create and empty aux function.
                //  Loop over nuclear centres and add the RNML contribution from each nucleus.
                for (auto atom:*cl) 
                    R.Add(RNLM(ab.Ltotal,ab.AlphaP,ab.P-atom->itsR), -1.0*(atom->itsZ) );

                auto NLMs=GaussianCD::GetNMLs(ab.Ltotal);
                const Polarization Pab = pa + pb;
                for (auto bNLM:NLMs)
                {
                    if (bNLM > Pab) continue;
                    if(double h = ab.H2(bNLM,pa,pb);h!=0) 
                        s+=h*R(bNLM);
                }
                        
                s*=2*Pi/ab.AlphaP*ab.Eij;
            }
            break;
        case qchem::InvOverlap :
        case qchem::InvRepulsion :
        case qchem::Charge :
        case qchem::Normalization :
            std::cerr << "GaussianRF::Integrate switch case not handled." <<  std::endl;
    } //switch
  
    return s;
}

//---------------------------------------------------------------------------------------
//
//  Do 3 center contractions
//
//  this is the c center. <ab|c>
double GaussianRF::Integrate(qchem::IType3C type,const RadialFunction* ra, const RadialFunction* rb, const Polarization& pa, const Polarization& pb, const Polarization& pc,CDCache& cache) const
{
    const GaussianRF* ga=dynamic_cast<const GaussianRF*>(ra);
    if (!ga) 
        return ra->Integrate(type,rb,pb,pa,pc,cache,this);
    const GaussianRF* gb=dynamic_cast<const GaussianRF*>(rb);
    if (!gb) 
        return rb->Integrate(type,ra,pa,pb,pc,cache,this);
    return Integrate3C(type,ga,gb,pa,pb,pc,cache,this);
}

//this = rb
double GaussianRF::Integrate(qchem::IType3C type,const RadialFunction* ra, const Polarization& pa, const Polarization& pb, const Polarization& pc,CDCache& cache,const RadialFunction* rc) const
{
    const GaussianRF* ga=dynamic_cast<const GaussianRF*>(ra);
    if (!ga) 
        return ra->Integrate(type,this,pb,pa,pc,cache,rc);
    const GaussianRF* gc=dynamic_cast<const GaussianRF*>(rc);
    assert(gc);
    return Integrate3C(type,ga,this,pa,pb,pc,cache,gc);
    
}

double GaussianRF::Integrate3C(qchem::IType3C type,grf_t* ga,grf_t* gb, po_t& pa, po_t& pb, po_t& pc,CDCache& cache, grf_t* gc)
{
    assert(ga);
    assert(gb);
    assert(gc);
    double s=0.0;
    switch (type)
    {
        case qchem::Overlap3C :
            {
                Hermite3* H3=gc->GetH3(*ga,*gb);
                s=(*H3)(pa,pb,pc);
                delete H3;
            }
            break;
        case qchem::Repulsion3C :
            {
                const GaussianCD& ab(cache.find(ga,gb));
                const RNLM&        R(cache.find(ab,gc));

                auto  NLMs=GaussianCD::GetNMLs(ab.Ltotal);
                const Hermite1& Hc=gc->GetH1();
                assert(&Hc);
                const Polarization Pab = pa+pb;
                for (auto nlm:NLMs)
                {
                    
                    if (nlm>Pab) continue;
                    double hab = ab.H2(nlm,pa,pb);
                    if (hab==0.0) continue;
                    double Rs=0.0;
                    for (int n=0; n<=pc.n; n++)
                            for (int l=0; l<=pc.l; l++)
                                for (int m=0; m<=pc.m; m++)
                                {
                                    Polarization NLMp(n,l,m);
                                    if (double h=Hc(NLMp,pc);h!=0.0)
                                        Rs+=h*R(nlm+NLMp);
                                } //for m
                    if (Rs!=0) s+=hab*Rs;
                   
                } //for (auto nlm:NLMs)
                
                double factor=1.0/(ab.AlphaP*gc->itsExponent*sqrt(ab.AlphaP+gc->itsExponent));
                factor = (pc.GetTotalL()%2) ? -factor : factor;
                s*=2*Pi52 * ab.Eij*factor;
            }
            
            break;
    }
    return s;

}



// this is rd
double GaussianRF::Integrate(rf_t* ra,rf_t* rb,rf_t* rc,po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache& cache) const
{
    const GaussianRF* gc=dynamic_cast<const GaussianRF*>(rc);
    if (!gc) 
        return rc->Integrate(ra,rb,pa,pb,pc,pd,cache,this);
    const GaussianRF* gb=dynamic_cast<const GaussianRF*>(rb);
    if (!gb) 
        return rb->Integrate(ra,pa,pb,pc,pd,cache,rc,this);
    const GaussianRF* ga=dynamic_cast<const GaussianRF*>(ra);
    if (!ga) 
        return ra->Integrate(rb,pb,pa,pc,pd,cache,rc,this);
    return Integrate4C(ga,gb,pa,pb,pc,pd,cache,gc,this);
}

// this = rc
double GaussianRF::Integrate(rf_t* ra,rf_t* rb, po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache& cache, rf_t* rd) const
{
    const GaussianRF* gb=dynamic_cast<const GaussianRF*>(rb);
    if (!gb) 
        return rb->Integrate(ra,pa,pb,pc,pd,cache,this,rd);
    const GaussianRF* ga=dynamic_cast<const GaussianRF*>(ra);
    if (!ga) 
        return ra->Integrate(rb,pb,pa,pc,pd,cache,this,rd);
    
    const GaussianRF* gd=dynamic_cast<const GaussianRF*>(rd);
    assert(gd);
    
    return Integrate4C(ga,gb,pa,pb,pc,pd,cache,this,gd);

}
// this = rb
double GaussianRF::Integrate(rf_t* ra, po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache& cache, rf_t* rc, rf_t* rd) const
{
    const GaussianRF* ga=dynamic_cast<const GaussianRF*>(ra);
    if (!ga) 
        return ra->Integrate(this,pb,pa,pc,pd,cache,rc,rd);

    const GaussianRF* gc=dynamic_cast<const GaussianRF*>(rc);
    assert(gc);
    const GaussianRF* gd=dynamic_cast<const GaussianRF*>(rd);
    assert(gd);
    
    return Integrate4C(ga,this,pa,pb,pc,pd,cache,gc,gd);
}

double GaussianRF::Integrate4C(grf_t* ga,grf_t* gb, po_t& pa, po_t& pb, po_t& pc, po_t& pd,CDCache& cache, grf_t* gc, grf_t* gd)
{
    assert(ga);
    assert(gb);
    assert(gc);
    assert(gd);
    
    const GaussianCD& ab(cache.find(ga,gb));
    const GaussianCD& cd(cache.find(gc,gd));

//    std::cout.precision(5);
//    std::cout.width(8);

    const std::vector<Polarization>& abNLMs=GaussianCD::GetNMLs(ab.Ltotal);
    const std::vector<Polarization>& cdNLMs=GaussianCD::GetNMLs(cd.Ltotal);


    double lambda=2*Pi52/(ab.AlphaP*cd.AlphaP*sqrt(ab.AlphaP+cd.AlphaP)); //M&D 3.31
    lambda*=ab.Eij*cd.Eij; //M&D 2.25
    const RNLM& rnlm(cache.find(ab,cd)); //M&D section 4A

    double s=0.0;
    const Polarization Pab = pa + pb;
    const Polarization Pcd = pc + pd;
    // NLM loop for ab
    for (auto abNLM:abNLMs)
    {
        if (abNLM > Pab) continue;
        double hab = ab.H2(abNLM,pa,pb);
        if (hab==0.0) continue;
        for (auto cdNLM:cdNLMs)
        {
            if (cdNLM > Pcd) continue;
            double hcd = cd.H2(cdNLM,pc,pd);
            if (hcd==0) continue;

            double r=rnlm(abNLM + cdNLM);
            if(r!=0)
                s+=hab*hcd*r*cdNLM.GetSign();;
        } //cdNLM
    } //abNLM
    s=s*lambda;

    return s;
}




double GaussianRF::GetKinetic(const Polarization& p1,const Polarization& p2,const GaussianCD& ab) const
{
    static Polarization p0(0,0,0),x(1,0,0),y(0,1,0),z(0,0,1);
    const  Polarization& P1= p1;
    const  Polarization& P2= p2;

    double txx=   P1.n * P2.n * ab.H2(p0,P1-x,P2-x)
                  - 2*P1.n * ab.b    * ab.H2(p0,P1-x,P2+x)
                  - 2*P2.n * ab.a    * ab.H2(p0,P1+x,P2-x)
                  + 4*ab.ab          * ab.H2(p0,P1+x,P2+x);

    double tyy=   P1.l * P2.l * ab.H2(p0,P1-y,P2-y)
                  - 2*P1.l * ab.b    * ab.H2(p0,P1-y,P2+y)
                  - 2*P2.l * ab.a    * ab.H2(p0,P1+y,P2-y)
                  + 4*ab.ab          * ab.H2(p0,P1+y,P2+y);

    double tzz=   P1.m * P2.m * ab.H2(p0,P1-z,P2-z)
                  - 2*P1.m * ab.b    * ab.H2(p0,P1-z,P2+z)
                  - 2*P2.m * ab.a    * ab.H2(p0,P1+z,P2-z)
                  + 4*ab.ab          * ab.H2(p0,P1+z,P2+z);

    return txx+tyy+tzz;
}
//---------------------------------------------------------------------------------------
//
//  Calculate 3 center hermite functions.
//  Here *this is treated as the third argument, center C.
//
Hermite3* GaussianRF::GetH3(const RadialFunction& r1, const RadialFunction& r2) const
{
    const GaussianRF* g1 = dynamic_cast<const GaussianRF*>(&r1);
    const GaussianRF* g2 = dynamic_cast<const GaussianRF*>(&r2);
    assert(g1);
    assert(g2);

    const double& a =g1->itsExponent, b =g2->itsExponent, c =itsExponent;
    const RVec3  & A =g1->GetCenter(), B =g2->GetCenter(), C =GetCenter();
    const int   & La=g1->GetL     (), Lb=g2->GetL     () ,Lc=GetL     ();

    double alphaQ = a+b+c;

    RVec3   AB=A-B;
    RVec3   AC=A-C;
    RVec3   BC=B-C;
    RVec3   Q = (a*A+b*B+c*C)/alphaQ;

    double Eabc = pow(Pi/alphaQ,1.5)*exp( -(a*b*AB*AB + a*c*AC*AC + b*c*BC*BC) / alphaQ );

    return new GaussianH3(alphaQ,Q-A,Q-B,Q-C,La,Lb,Lc,Eabc);
}

//----------------------------------------------------------------------------
//
// double factorial table, starts at -1, so you have to add 1 to the index.
//
double DoubleFactData[14] = {1,1,1,2,3,8,15,48,105,384,945,3840,10395,46080};
inline double DoubleFact(int i)
{
    return DoubleFactData[i+1];
}


double GaussianRF::GetNormalization(const Polarization& p) const
{
    assert(2*p.n-1 <= 12);
    assert(2*p.l-1 <= 12);
    assert(2*p.m-1 <= 12);
    double s=pow(Pi/(2*itsExponent),1.5);
    double t=intpow(4*itsExponent,p.GetTotalL());
    double f=DoubleFact(2*p.n-1) * DoubleFact(2*p.l-1) * DoubleFact(2*p.m-1);
    return sqrt(t/(s*f));
}

double GaussianRF::GetCharge(const Polarization& p) const
{
    assert(2*p.n-1 <= 12);
    assert(2*p.l-1 <= 12);
    assert(2*p.m-1 <= 12);
    if ((p.n%2) || (p.l%2) || (p.m%2)) return 0;
    double s=pow(Pi/itsExponent,1.5);
    double t=intpow(2*itsExponent,p.GetTotalL()/2);
    double f=DoubleFact(p.n-1) * DoubleFact(p.l-1) * DoubleFact(p.m-1);
    return s*f/t;
}


//--------------------------------------------------------------------------------
//
//  Streamable object stuff.
//
std::ostream& GaussianRF::Write(std::ostream& os) const
{
    if (Binary())
    {
        BinaryWrite(itsExponent,os);
    }
    if (Ascii ())
        os << itsExponent << " ";

    if (Pretty())
    {
        os << "Primative  " << std::setw(8) << itsExponent;
    }

    if (!Pretty()) RadialCommon::Write(os);

    return os;
}

std::istream& GaussianRF::Read(std::istream& is)
{
    if (Binary())
    {
        BinaryRead(itsExponent,is);
    }
    else
    {
        is >> itsExponent;
    }
    RadialCommon::Read(is);
    return is;
}

RadialFunction* GaussianRF::Clone() const
{
    return new  GaussianRF(*this);
}

RadialFunction* GaussianRF::Clone(const RVec3& newCenter) const
{
    return new GaussianRF(itsExponent,newCenter,GetL());
}

//----------------------------------------------------------------
//
//  Scalar function stuff.
//
double GaussianRF::operator()(const RVec3& r) const
{
    RVec3 dr=GetCenter()-r;
    return exp(-itsExponent*dr*dr);
}

RVec3 GaussianRF::Gradient(const RVec3& r) const
{
    RVec3 dr=GetCenter()-r;
    return -2*itsExponent* (*this)(r) * dr;
}

} //namespace PolarizedGaussian
