// File: GaussianRF.C  Primative Gaussian in 3D space.


#include "oml/smatrix.h"
#include "BasisSetImplementation/PolarizedGaussian/Gaussian/GaussianRF.H"
#include "BasisSetImplementation/PolarizedGaussian/Gaussian/GaussianCD.H"
#include "BasisSetImplementation/PolarizedGaussian/Gaussian/GaussianH3.H"
#include "BasisSetImplementation/PolarizedGaussian/BasisFunctionBlock.H"
#include "Misc/Polarization.H"
#include "Misc/ERIList.H"
#include "BasisSetImplementation/PolarizedGaussian/Hermite/Hermite1.H"
#include "BasisSetImplementation/PolarizedGaussian/Auxillary/RNLM.H"
#include "Mesh/MeshBrowser.H"
#include "Cluster/ClusterBrowser.H"
#include "oml/imp/binio.h"
#include <cmath>
#include <iostream>
#include <iomanip>
#include <cassert>
#include <stdlib.h>

//#######################################################################
//
//   Gaussian radial function implementation
//

GaussianRF::GaussianRF()
    : RadialFunctionImplementation()
    , itsExponent(0)
{};

GaussianRF::GaussianRF(double theExponent, const RVec3& theCenter, int theL)
    : RadialFunctionImplementation(theCenter,theL)
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
        bool base     = RadialFunctionImplementation::operator==(rf);
        bool exponent = fabs((itsExponent-g->itsExponent)/(itsExponent+g->itsExponent)*2.0) < 0.001; // 0.1%.
        ret= exponent && base;
    }
    return ret;
}

//---------------------------------------------------------------------------------------
//
//  Do 3 center contractions
//

void GaussianRF::Get2CenterIntegrals(Types2C type, BFBP& p, SMat& ret, const Cluster* cl, double scale) const
{
    p.radials.push_back(this);
    switch (p.radials.itsIndex)
    {
    case 0:
        p.b->itsRadial->Get2CenterIntegrals(type,p,ret,cl,scale);
        break;
    case 1:
        switch (type)
        {
        case Overlap2C:
            GetOverlap2CInternal(p,ret,scale);
            break;
        case Repulsion2C:
            GetRepulsion2CInternal(p,ret,scale);
            break;
        case Kinetic:
            GetKinetic2CInternal(p,ret,scale);
            break;
        case Nuclear:
            GetNuclear2CInternal(p,ret,cl,scale);
            break;
        }
        break;
    default:
        std::cerr << "GaussianRF::Get2CenterIntegrals Unhandeled index: " << p.radials.itsIndex << std::endl;
        assert(false);
        exit(-1);
    }
    p.radials.pop_back();
}

void GaussianRF::Get2CenterIntegrals(Types2C type, BFBP& p, Mat& ret, double scale) const
{
    p.radials.push_back(this);
    switch (p.radials.itsIndex)
    {
    case 0:
        p.b->itsRadial->Get2CenterIntegrals(type,p,ret,scale);
        break;
    case 1:
        switch (type)
        {
        case Overlap2C:
            GetOverlap2CInternal(p,ret,scale);
            break;
        case Repulsion2C:
            GetRepulsion2CInternal(p,ret,scale);
            break;
        default:
            assert(false);
            break;
        }
        break;
    default:
        std::cerr << "GaussianRF::Get2CenterIntegrals Unhandeled index: " << p.radials.itsIndex << std::endl;
        assert(false);
        exit(-1);
    }
    p.radials.pop_back();
}

//---------------------------------------------------------------------------------------
//
//  Do 3 center contractions
//

void GaussianRF::Get3CenterIntegrals(Types3C type, BFBT& t, std::vector<SMat >& ret, double scale)
{
    t.radials.push_back(this);
    switch (t.radials.itsIndex)
    {
    case 0:
        t.b->itsRadial->Get3CenterIntegrals(type,t,ret,scale);
        break;
    case 1:
        t.c->itsRadial->Get3CenterIntegrals(type,t,ret,scale);
        break;
    case 2:
        switch (type)
        {
        case Overlap3C:
            GetOverlap3CInternal(t,ret,scale);
            break;
        case Repulsion3C:
            GetRepulsion3CInternal(t,ret,scale);
            break;
        }
        break;
    default:
        std::cerr << "GaussianRF::GetRepulsion3C Unhandeled index: " << t.radials.itsIndex << std::endl;
        assert(false);
        exit(-1);
    }
    t.radials.pop_back();
}

void GaussianRF::GetRepulsion4C(BFBQ& q, ERIList& eris, double scale)
{
    q.radials.push_back(this);
    switch (q.radials.itsIndex)
    {
    case 0:
        q.b->itsRadial->GetRepulsion4C(q,eris,scale);
        break;
    case 1:
        q.c->itsRadial->GetRepulsion4C(q,eris,scale);
        break;
    case 2:
        q.d->itsRadial->GetRepulsion4C(q,eris,scale);
        break;
    case 3:
        GetRepulsion4CInternal(q,eris,scale);
        break;
    default:
        std::cerr << "GaussianRF::GetRepulsion4C Unhandeled index: " << q.radials.itsIndex << std::endl;
        assert(false);
        exit(-1);
    }
    q.radials.pop_back();
}


void GaussianRF::GetOverlap2CInternal  (BFBP& p, SMat& ret, double scale) const
{
    assert(p.radials.a());
    assert(p.radials.b());

    Polarization zero(0,0,0);
    GaussianCD ab(*p.radials.a(),*p.radials.b());
    double factor=scale*pow(Pi/ab.AlphaP,1.5)*ab.Eij;
    SMatrix<double>::Subscriptor s(ret);

    std::vector<Polarization>::const_iterator bPa(p.a->itsPols.begin());
    for (int nPa=p.a->itsN;  bPa!=p.a->itsPols.end();  bPa++,nPa++)
    {
        std::vector<Polarization>::const_iterator bPb(p.b->itsPols.begin());
        for (int nPb=p.b->itsN;  bPb!=p.b->itsPols.end();  bPb++,nPb++)
        {
            double h = ab.H2(zero,*bPa,*bPb);
            if (h!=0 && nPa <= nPb) s(nPa,nPb)+=factor*h;
        }
    }

}

void GaussianRF::GetOverlap2CInternal  (BFBP& p, Mat& ret, double scale) const
{
    assert(p.radials.a());
    assert(p.radials.b());

    Polarization zero(0,0,0);
    GaussianCD ab(*p.radials.a(),*p.radials.b());
    double factor=scale*pow(Pi/ab.AlphaP,1.5)*ab.Eij;
    Matrix<double>::Subscriptor s(ret);

    std::vector<Polarization>::const_iterator bPa(p.a->itsPols.begin());
    for (int nPa=p.a->itsN;  bPa!=p.a->itsPols.end();  bPa++,nPa++)
    {
        std::vector<Polarization>::const_iterator bPb(p.b->itsPols.begin());
        for (int nPb=p.b->itsN;  bPb!=p.b->itsPols.end();  bPb++,nPb++)
        {
            double h = ab.H2(zero,*bPa,*bPb);
            if (h!=0) s(nPa,nPb)+=factor*h;
        }
    }

}

void GaussianRF::GetRepulsion2CInternal(BFBP& p, SMat& ret, double scale) const
{
    assert(p.radials.a());
    assert(p.radials.b());
    GaussianCD ab(*p.radials.a(),*p.radials.b());
    int La=p.radials.a()->GetL();
    const std::vector<Polarization>& NLMs=GaussianCD::GetNMLs(La);
    double Alpha = ab.a;
    Matrix<double> Rs=p.radials.b()->GetAux(NLMs,p.b->itsPols,La,Alpha,p.radials.a()->GetCenter());
    const Matrix<double>&         rsub(Rs);
    SMatrix<double>::Subscriptor s(ret);

    std::vector<Polarization>::const_iterator bNLM(NLMs.begin());
    for (int nNLM=1;  bNLM!=NLMs.end();  bNLM++,nNLM++)
    {
        std::vector<Polarization>::const_iterator bPa(p.a->itsPols.begin());
        for (int nPa=p.a->itsN;  bPa!=p.a->itsPols.end();  bPa++,nPa++)
        {

            if (bNLM->n <= bPa->n && bNLM->l <= bPa->l && bNLM->m <= bPa->m)
            {
                double h=p.radials.a()->GetH1()(*bNLM,*bPa);
                if (h!=0.0)
                {
                    h*=2*Pi52*scale;
                    for (int nPb=p.b->itsN; nPb<=p.b->itsN+p.b->itsPols.size()-1; nPb++)
                        if (nPa <= nPb) s(nPa,nPb) += h * rsub(nNLM,nPb-p.b->itsN+1);
                }
            }
        }
    }

}

void GaussianRF::GetRepulsion2CInternal(BFBP& p, Mat& ret, double scale) const
{
    assert(p.radials.a());
    assert(p.radials.b());
    GaussianCD ab(*p.radials.a(),*p.radials.b());
    int La=p.radials.a()->GetL();
    const std::vector<Polarization>& NLMs=GaussianCD::GetNMLs(La);
    double Alpha = ab.a;
    Matrix<double> Rs=p.radials.b()->GetAux(NLMs,p.b->itsPols,La,Alpha,p.radials.a()->GetCenter());
    const Matrix<double>&         rsub(Rs);
    Matrix<double>::Subscriptor s(ret);

    std::vector<Polarization>::const_iterator bNLM(NLMs.begin());
    for (int nNLM=1;  bNLM!=NLMs.end();  bNLM++,nNLM++)
    {
        std::vector<Polarization>::const_iterator bPa(p.a->itsPols.begin());
        for (int nPa=p.a->itsN;  bPa!=p.a->itsPols.end();  bPa++,nPa++)
        {

            if (bNLM->n <= bPa->n && bNLM->l <= bPa->l && bNLM->m <= bPa->m)
            {
                double h=p.radials.a()->GetH1()(*bNLM,*bPa);
                if (h!=0.0)
                {
                    h*=2*Pi52*scale;
                    for (int nPb=p.b->itsN; nPb<=p.b->itsN+p.b->itsPols.size()-1; nPb++)
                        s(nPa,nPb) += h * rsub(nNLM,nPb-p.b->itsN+1);
                }
            }
        }
    }

}
void GaussianRF::GetKinetic2CInternal  (BFBP& p, SMat& ret, double scale) const
{
    assert(p.radials.a());
    assert(p.radials.b());
    GaussianCD ab(*p.radials.a(),*p.radials.b());
    double factor=scale*0.5*pow(Pi/ab.AlphaP,1.5)*ab.Eij;
    SMatrix<double>::Subscriptor s(ret);
    std::vector<Polarization>::const_iterator bPa(p.a->itsPols.begin());
    for (int nPa=p.a->itsN;  bPa!=p.a->itsPols.end();  bPa++,nPa++)
    {
        std::vector<Polarization>::const_iterator bPb(p.b->itsPols.begin());
        for (int nPb=p.b->itsN;  bPb!=p.b->itsPols.end();  bPb++,nPb++)
        {
            double h = GetKinetic(*bPa,*bPb,ab);
            if (h!=0 && nPa <= nPb) s(nPa,nPb)+=factor*h;
        }
    }
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


void GaussianRF::GetNuclear2CInternal  (BFBP& p, SMat& ret, const Cluster* cl, double scale) const
{
    assert(p.radials.a());
    assert(p.radials.b());
    GaussianCD ab(*p.radials.a(),*p.radials.b());
    RNLM R; //Create and empty aux function.
    //  Loop over nuclear centers and add the RNML contribution from each nucleaus.
    for (ClusterBrowser cb(*cl); cb; cb++) R.Add(RNLM(ab.Ltotal,ab.AlphaP,ab.P-(*cb).itsR), -1.0*(cb->itsZ) );

    double factor=scale*2*Pi/ab.AlphaP*ab.Eij;
    SMatrix<double>::Subscriptor s(ret);
    const std::vector<Polarization>& NLMs=GaussianCD::GetNMLs(ab.Ltotal);

    for (std::vector<Polarization>::const_iterator bNLM(NLMs.begin());  bNLM!=NLMs.end();  bNLM++)
    {
        std::vector<Polarization>::const_iterator bPa(p.a->itsPols.begin());
        for (int nPa=p.a->itsN;  bPa!=p.a->itsPols.end();  bPa++,nPa++)
        {
            std::vector<Polarization>::const_iterator bPb(p.b->itsPols.begin());
            for (int nPb=p.b->itsN;  bPb!=p.b->itsPols.end();  bPb++,nPb++)
            {
                const Polarization Pab = *bPa + *bPb;
                if (bNLM->n <= Pab.n && bNLM->l <= Pab.l && bNLM->m <= Pab.m   && nPa <= nPb)
                {
                    double h = ab.H2(*bNLM,*bPa,*bPb);
                    if(h!=0) s(nPa,nPb)+=h*R(*bNLM)*factor;
                }
            }
        }
    }
}

void GaussianRF::GetOverlap3CInternal(BFBT& t, std::vector<SMat>& ret, double scale) const
{
    assert(t.radials.a());
    assert(t.radials.b());
    assert(t.radials.c());
    Hermite3* H3=t.radials.c()->GetH3(*t.radials.a(),*t.radials.b());
    
    std::vector<Polarization>::const_iterator bPa(t.a->itsPols.begin());
    for (int nPa=t.a->itsN;  bPa!=t.a->itsPols.end();  bPa++,nPa++)
    {
        std::vector<Polarization>::const_iterator bPb(t.b->itsPols.begin());
        for (int nPb=t.b->itsN;  bPb!=t.b->itsPols.end();  bPb++,nPb++)
        {
            if (nPa <= nPb)
            {
                std::vector<Polarization>::const_iterator bPc(t.c->itsPols.begin());
                for (int nPc=0;  bPc!=t.c->itsPols.end();  bPc++,nPc++)
                {
                    ret[nPc](nPa,nPb)+=scale*(*H3)(*bPa,*bPb,*bPc);
                }
            }
        }
    }
    delete H3;
}

//
// The t variable contains info on polarizations and the contracted radial functions.
// The radials argument contains all four primative radial functions.
// The scale argument should be a product of all contraction coefficients.
// The output gets stored in the ret matrix list.
//
void GaussianRF::GetRepulsion3CInternal(BFBT& t, std::vector<SMat>& ret, double scale) const
{
    assert(t.radials.a());
    assert(t.radials.b());
    assert(t.radials.c());

    GaussianCD ab(*t.radials.a(),*t.radials.b());

    const std::vector<Polarization>& NLMs=GaussianCD::GetNMLs(ab.Ltotal);

    Matrix<double> Rs=t.radials.c()->GetAux(NLMs,t.c->itsPols,ab.Ltotal,ab.AlphaP,ab.P);
    const Matrix<double>& rsub(Rs);

    std::vector<Polarization>::const_iterator bNLM(NLMs.begin());
    for (int nNLM=1;  bNLM!=NLMs.end();  bNLM++,nNLM++)
    {
        std::vector<Polarization>::const_iterator bPa(t.a->itsPols.begin());
        for (int nPa=t.a->itsN;  bPa!=t.a->itsPols.end();  bPa++,nPa++)
        {
            std::vector<Polarization>::const_iterator bPb(t.b->itsPols.begin());
            for (int nPb=t.b->itsN;  bPb!=t.b->itsPols.end();  bPb++,nPb++)
            {
                const Polarization Pab = *bPa + *bPb;
                if (bNLM->n <= Pab.n && bNLM->l <= Pab.l && bNLM->m <= Pab.m   &&   nPa <= nPb)
                {
//                    double hab = AmIFlipped() ? H2(*bNLM,*bPb,*bPa) : H2(*bNLM,*bPa,*bPb);
                    double hab = ab.H2(*bNLM,*bPa,*bPb);
                    hab*=2*Pi52 * ab.Eij * scale;
                    if (hab!=0)
                    {
                        std::vector<Polarization>::const_iterator bPc(t.c->itsPols.begin());
                        for (int nPc=1;  bPc!=t.c->itsPols.end();  bPc++,nPc++)
                        {
                            double r=rsub(nNLM,nPc);
                            if(r!=0)
                            {
                                ret[nPc-1](nPa,nPb)+=hab*r;
                            }
                        }
                    }
                }
            }
        }
    }
}
//
// The q variable contains info on polarizations and the contracted radial functions.
// The radials argument contains all four primative radial functions.
// The scale argument should be a product of all contraction coefficients.
// The output gets stored in the eris super matrix.
//
#include "BasisSetImplementation/SphericalGaussian/GaussianIntegrals.H"

void GaussianRF::GetRepulsion4CInternal( BFBQ& q, ERIList& eris, double scale) const
{
    assert(q.radials.a());
    assert(q.radials.b());
    assert(q.radials.c());
    assert(q.radials.d());

    GaussianCD ab(*q.radials.a(),*q.radials.b());
    GaussianCD cd(*q.radials.c(),*q.radials.d());

//    std::cout.precision(5);
//    std::cout.width(8);

    const std::vector<Polarization>& abNLMs=GaussianCD::GetNMLs(ab.Ltotal);
    const std::vector<Polarization>& cdNLMs=GaussianCD::GetNMLs(cd.Ltotal);

    typedef std::vector<Polarization>::const_iterator PolITER;

    double lambda=2*Pi52/(ab.AlphaP*cd.AlphaP*sqrt(ab.AlphaP+cd.AlphaP)); //M&D 3.31
    lambda*=ab.Eij*cd.Eij*scale; //M&D 2.25
    double alpha=ab.AlphaP*cd.AlphaP/(ab.AlphaP+cd.AlphaP); //M&D 3.32
    RVec3 PQ = ab.P-cd.P; //M&D 3.32
    RNLM rnlm(ab.Ltotal+cd.Ltotal,alpha,PQ); //M&D section 4A

    // nlm loop for a
    int nPa=q.a->itsN;
    for (PolITER bPa=q.a->itsPols.begin();  bPa!=q.a->itsPols.end();  bPa++,nPa++)
    {
        // nlm loop for b
        int nPb=q.b->itsN;
        PolITER bPb=q.b->itsPols.begin();
        // Don't go below the diagonal if the a and b radials are the same.
        if (q.a->itsRadial->GetID()==q.b->itsRadial->GetID())
        {
            bPb=bPa;
            nPb=nPa;
        }
        for (;  bPb!=q.b->itsPols.end();  bPb++,nPb++)
        {
            int nPc=q.c->itsN;
            PolITER bPc=q.c->itsPols.begin();
            // Don't go below the diagonal if the a and c radials are the same.
            if (q.a->itsRadial->GetID()==q.c->itsRadial->GetID())
            {
                bPc=bPa;
                nPc=nPa;
            }
            for (;  bPc!=q.c->itsPols.end();  bPc++,nPc++)
            {
                // nlm loop for d
                int nPd=q.d->itsN;
                PolITER bPd=q.d->itsPols.begin();
                PolITER bPd_end=q.d->itsPols.end();
                if (nPa==nPc && q.b->itsRadial->GetID()==q.d->itsRadial->GetID())
                {
                    bPd=bPb;
                    nPd=nPb;
                }
                else if (nPa==nPc && nPd<nPb && q.b->itsRadial->GetID()!=q.d->itsRadial->GetID())
                {
//                    std::cout << "Skipping a,b,c,d=" << nPa << ", " <<  nPb << ", " <<  nPc << ", " <<  nPd << std::endl;
                    continue;

                }
                // Don't go below the diagonal if the d and c radials are the same.
                else if (q.c->itsRadial->GetID()==q.d->itsRadial->GetID())
                {
                    bPd=bPc;
                    nPd=nPc;
                }


                for (;  bPd!=bPd_end;  bPd++,nPd++)
                {
//                    if (tracker_eris(nPa,nPb,nPc,nPd)==-1.0)
//                    {
//                        std::cout << "Allready assigned a,b,c,d=" << nPa << ", " <<  nPb << ", " <<  nPc << ", " <<  nPd << std::endl;
//                        assert(false);
//                    }
//                    tracker_eris(nPa,nPb,nPc,nPd)=-1.0;
                    assert(nPa<=nPb);
                    assert(nPa<=nPc);
                    if (nPa==nPc)
                        assert(nPb<=nPd);
                    else
                        assert(nPc<=nPd);

                    double ERItemp=0.0;
                    const Polarization Pab = *bPa + *bPb;
                    const Polarization Pcd = *bPc + *bPd;
                    // NLM loop for ab
                    for (PolITER abNLM=abNLMs.begin();  abNLM!=abNLMs.end();  abNLM++)
                    {
                        if (!(abNLM->n <= Pab.n && abNLM->l <= Pab.l && abNLM->m <= Pab.m   &&   nPa <= nPb)) continue;
                        double hab = ab.H2(*abNLM,*bPa,*bPb);
                        if (hab==0.0) continue;
                        //hab*=lambda;
                        // NLM loop for cd
                        for (PolITER cdNLM=cdNLMs.begin();  cdNLM!=cdNLMs.end();  cdNLM++)
                        {
                            double sign=cdNLM->GetSign();
                            if (!(cdNLM->n <= Pcd.n && cdNLM->l <= Pcd.l && cdNLM->m <= Pcd.m   &&   nPc <= nPd)) continue;

                            double hcd = cd.H2(*cdNLM,*bPc,*bPd);
                            if (hcd==0) continue;

                            double r=rnlm(*abNLM + *cdNLM);
                            if(r!=0)
                            {
                                ERItemp+=hab*hcd*r*sign;
                            }
                        } //cdNLM
                    } //abNLM
                    eris(nPa,nPb,nPc,nPd)+=ERItemp*lambda;

                } //d

            } //c

        } //b
    } //a



//
//  Make sure we got all the trackers
//
//    {
//        int nPa=q.a->itsN;
//        for (PolITER bPa=q.a->itsPols.begin();  bPa!=q.a->itsPols.end();  bPa++,nPa++)
//        {
//            int nPb=q.b->itsN;
//            for (PolITER bPb=q.b->itsPols.begin();  bPb!=q.b->itsPols.end();  bPb++,nPb++)
//            {
//                int nPc=q.c->itsN;
//                for (PolITER bPc=q.c->itsPols.begin();  bPc!=q.c->itsPols.end();  bPc++,nPc++)
//                {
//                    int nPd=q.d->itsN;
//                    for (PolITER bPd=q.d->itsPols.begin();  bPd!=q.d->itsPols.end();  bPd++,nPd++)
//                    {
//                        if (tracker_eris(nPa,nPb,nPc,nPd)!=-1.0)
//                        {
//                            std::cout << "Never assigned a,b,c,d=" << nPa << ", " <<  nPb << ", " <<  nPc << ", " <<  nPd << std::endl;
//                            assert(false);
//                        }
//                    }
//
//                }
//            }
//        }
//    }
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


//---------------------------------------------------------------------------------------
//
//  Calculate 3 center Auxillar  functions.
//
Matrix<double> GaussianRF::GetAux(const std::vector<Polarization>& NLM, const std::vector<Polarization>& Pc,
                                  int LP, double AlphaP, const RVec3& P) const
{
    double alpha =AlphaP*itsExponent/(AlphaP+itsExponent);
    RNLM R(LP+GetL(),alpha,P-GetCenter());

    Matrix<double> ret(NLM.size(),Pc.size());
    Fill(ret,0.0);
    Matrix<double>::Subscriptor s(ret);

    const Hermite1& H1=GetH1();
    assert(&H1);
    double factor=1.0/(AlphaP*itsExponent*sqrt(AlphaP+itsExponent));
    std::vector<Polarization>::const_iterator bPc(Pc.begin());
    for (int nPc=1; bPc!=Pc.end(); bPc++,nPc++)
    {
        double sign = (bPc->GetTotalL()%2) ? -factor : factor;

        for (int n=0; n<=(*bPc).n; n++)
            for (int l=0; l<=(*bPc).l; l++)
                for (int m=0; m<=(*bPc).m; m++)
                {
                    Polarization NLMp(n,l,m);
                    double h=H1(NLMp,*bPc);
                    if (h!=0.0)
                    {
                        h*=sign;
                        std::vector<Polarization>::const_iterator bNLM(NLM.begin());
                        for (int nNLM=1; bNLM!=NLM.end(); bNLM++,nNLM++) s(nNLM,nPc)+=h*R(*bNLM+NLMp);
                    }
                }
    }
    return ret;
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

    if (!Pretty()) RadialFunctionImplementation::Write(os);

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
    RadialFunctionImplementation::Read(is);
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






