// File: SphericalGaussianIE.C  Here is where all the integral get calculated.


#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianBF.H"
#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianBS.H"
#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianIE.H"
#include "BasisSetImplementation/SphericalGaussian/SphericalSymmetryQN.H"
#include "BasisSetImplementation/SphericalGaussian/GaussianIntegrals.H"
#include "BasisSetImplementation/SphericalGaussian/SlaterIntegrals.H"
#include "BasisSet/IntegralDataBase.H"
//#include "BasisSet/TBasisSetBrowser.H"
#include "BasisSet/BasisGroup.H"
#include "Cluster/Cluster.H"
#include "Misc/MatrixList.H"
#include "Misc/ERIProxy.H"
#include <cassert>
#include <iostream>
#include <stdlib.h>

double FourPi2=4*4*Pi*Pi;

//-----------------------------------------------------------------
//
//  Construction zone.
//
SphericalGaussianIE::SphericalGaussianIE()
    : IntegralEngineImplementation<double>()
    , itsL(0)
    , itsExponents()
{};

void SphericalGaussianIE::Insert(const TBasisSet<double>* theSet)
{
    IntegralEngineImplementation<double>::Insert(theSet);
    itsExponents.SetLimits(VecLimits(itsN));
    InitializeLocalStuff();
}

void SphericalGaussianIE::InitializeLocalStuff()
{
    assert(itsExponents.size()==itsBasisSet->GetNumFunctions());
    int i=1;
    for (auto b:*itsBasisSet)
    {
        auto bsg=dynamic_cast<const SphericalGaussianBF*>(b);
        assert(bsg);
        if (b==*itsBasisSet->begin()) 
            itsL=bsg->itsL;
        else
            if (bsg->itsL != itsL)
            {
                std::cerr << "SphericalGaussianIE::InitializeLocalStuff basis set has more than one L quantum number" << std::endl;
                exit(-1);
            }
        itsExponents(i++)=bsg->itsExponent;
        
    }
}

//-----------------------------------------------------------------
//
//  Streamable Object stuff
//
IntegralEngine<double>* SphericalGaussianIE::Clone() const
{
    return new SphericalGaussianIE(*this);
}
//----------------------------------------------------------------------------------------
//
//  Overlap type integrals
//
SMatrix<double> SphericalGaussianIE::MakeOverlap() const
{
    CheckBasisSet();
    SMatrix<double> ret(itsN,itsN);
    SMatrix<double>::Subscriptor s(ret);

    double a=TwoLPlusOne(itsL);
    for (index_t i=1; i<=itsN; i++)
        for (index_t j=i; j<=itsN; j++)
            s(i,j)=a*GaussianIntegral(itsExponents(i)+itsExponents(j),2*itsL);

    Normalize(ret);
    return ret;
}

Matrix<double> SphericalGaussianIE::MakeOverlap(const TBasisSet<double>& theBasisSet) const
{
    CheckBasisSet();
    const SphericalGaussianIE* OtherIE;
    const IntegralEngine<double>*      theOtherIE;

    theOtherIE=theBasisSet.GetDataBase()->GetIntegralEngine();
    OtherIE  = dynamic_cast<const SphericalGaussianIE*>(theOtherIE);
    assert(OtherIE);
    index_t otherN=OtherIE->itsN;
    int     otherL=OtherIE->itsL;

    Matrix<double> ret(itsN,otherN);
    Matrix<double>::Subscriptor s(ret);


    const Vector<double>& e1(itsExponents);
    const Vector<double>& e2(OtherIE->itsExponents);
    const Vector<double>& n1(itsNormalizations);
    const Vector<double>& n2(OtherIE->itsNormalizations);

    for (index_t i=1; i<=itsN; i++)
        for (index_t j=1; j<=otherN; j++)
            s(i,j)=GaussianIntegral(e1(i)+e2(j),itsL+otherL)*n1(i)*n2(j);

    return ret;
}


Vector<double> SphericalGaussianIE::MakeOverlap(const ScalarFunction<double>& f) const
{
    CheckBasisSet();
    Vector<double> ret(itsN);
    const SphericalGaussianBF* bfag= dynamic_cast<const SphericalGaussianBF*>(&f);
    assert(bfag);
    assert(bfag->itsL==0);
    assert(itsL==0);

    Vector<double>::Subscriptor s(ret);
    const Vector<double>& e(itsExponents);

    for (index_t i=1; i<=itsN; i++) s(i) = GaussianIntegral(e(i)+bfag->itsExponent,itsL+bfag->itsL);

    Normalize(ret);
    ret *= bfag->GetNormalization();
    return ret;
}

void SphericalGaussianIE::MakeOverlap3C(MList& mlist, const TBasisSet<double>& bs) const
{
    mlist.Empty();
    for (auto b=bs.beginT(); b!=bs.end(); b++) mlist.Add(MakeOverlap(**b));
    mlist.Clear();
}

//----------------------------------------------------------------------------------------
//
//  Repulsion type integrals
//
SMatrix<double> SphericalGaussianIE::MakeRepulsion() const
{
    CheckBasisSet();
    SMatrix<double> ret(itsN,itsN);
    SMatrix<double>::Subscriptor s(ret);

    const Vector<double>& e(itsExponents);

    double a=TwoLPlusOne(itsL);
    for (index_t i=1; i<=itsN; i++)
        for (index_t j=i; j<=itsN; j++)
            s(i,j)=a*GaussianRepulsionIntegral(e(i),e(j),itsL,itsL);

    Normalize(ret);
    return ret;
}

Matrix<double> SphericalGaussianIE::MakeRepulsion(const TBasisSet<double>& theBasisSet) const
{
    CheckBasisSet();
    const SphericalGaussianIE* OtherIE;
    const IntegralEngine<double>*      theOtherIE;

    theOtherIE=theBasisSet.GetDataBase()->GetIntegralEngine();
    OtherIE  = dynamic_cast<const SphericalGaussianIE*>(theOtherIE);
    assert(OtherIE);
    index_t otherN=OtherIE->itsN;
    int     otherL=OtherIE->itsL;

    Matrix<double> ret(itsN,otherN);
    Matrix<double>::Subscriptor s(ret);


    const Vector<double>& e1(itsExponents);
    const Vector<double>& e2(OtherIE->itsExponents);
    const Vector<double>& n1(itsNormalizations);
    const Vector<double>& n2(OtherIE->itsNormalizations);

    for (index_t i=1; i<=itsN; i++)
        for (index_t j=1; j<=otherN; j++)
            s(i,j)=GaussianRepulsionIntegral(e1(i),e2(j),itsL,otherL)*n1(i)*n2(j);

    return ret;
}



Vector<double> SphericalGaussianIE::MakeRepulsion(const ScalarFunction<double>& f) const
{
    CheckBasisSet();
    Vector<double> ret(itsN);
    const SphericalGaussianBF* bfag= dynamic_cast<const SphericalGaussianBF*>(&f);
    assert(bfag);
    assert(bfag->itsL==0);
    assert(itsL==0);

    Vector<double>::Subscriptor s(ret);
    const Vector<double>& e(itsExponents);

    for (index_t i=1; i<=itsN; i++) s(i) = GaussianRepulsionIntegral(e(i),bfag->itsExponent,itsL,bfag->itsL);

    Normalize(ret);
    ret *= bfag->GetNormalization();
    return ret;
}


void SphericalGaussianIE::MakeRepulsion3C(MList& mlist, const TBasisSet<double>& bs) const
{
    mlist.Empty();
    for (auto b=bs.beginT(); b!=bs.end(); b++) mlist.Add(MakeRepulsion(**b));
    mlist.Clear();
}

//
//  This is where we do the big double loop over basis sets.
//
void SphericalGaussianIE::MakeRepulsion4C(ERIList& Coulomb, ERIList& exchange, const BasisGroup* bg) const
{
    // TODO count integrals and % non zero.
    Coulomb.SetSize(bg->GetNumFunctions());
    exchange.SetSize(bg->GetNumFunctions());
    ERIList tracker_eris(Coulomb.GetSize());
    for (auto ba=bg->begin();ba!=bg->end();ba++)
    {
        int start_a=(*ba)->GetStartIndex();
        const SphericalGaussianBS* sg_a=dynamic_cast<const SphericalGaussianBS*>(*ba);
        assert(sg_a);
        const SphericalGaussianIE* ie_a=dynamic_cast<const SphericalGaussianIE*>(sg_a->GetDataBase()->GetIntegralEngine());
        assert(ie_a);

        const Vector<double>& ea(ie_a->itsExponents);
        const Vector<double>& na(ie_a->itsNormalizations);
        int Na=ie_a->itsN;
        int La=ie_a->itsL;
        for (auto bb=ba;bb!=bg->end();bb++)
        {
            int start_b=(*bb)->GetStartIndex();
            const SphericalGaussianBS* sg_b=dynamic_cast<const SphericalGaussianBS*>(*bb);
            assert(sg_b);
            const SphericalGaussianIE* ie_b=dynamic_cast<const SphericalGaussianIE*>(sg_b->GetDataBase()->GetIntegralEngine());
            assert(ie_b);

            const Vector<double>& eb(ie_b->itsExponents);
            const Vector<double>& nb(ie_b->itsNormalizations);
            int Nb=ie_b->itsN;
            int Lb=ie_b->itsL;

            for (auto bc=ba;bc!=bg->end();bc++)
            {
                int start_c=(*bc)->GetStartIndex();
                const SphericalGaussianBS* sg_c=dynamic_cast<const SphericalGaussianBS*>(*bc);
                assert(sg_c);
                const SphericalGaussianIE* ie_c=dynamic_cast<const SphericalGaussianIE*>(sg_c->GetDataBase()->GetIntegralEngine());
                assert(ie_c);

                const Vector<double>& ec(ie_c->itsExponents);
                const Vector<double>& nc(ie_c->itsNormalizations);
                int Nc=ie_c->itsN;
                int Lc=ie_c->itsL;
                
                for (auto bd=bc;bd!=bg->end();bd++)
                {
                    int start_d=(*bd)->GetStartIndex();
                    const SphericalGaussianBS* sg_d=dynamic_cast<const SphericalGaussianBS*>(*bd);
                    assert(sg_d);
                    const SphericalGaussianIE* ie_d=dynamic_cast<const SphericalGaussianIE*>(sg_d->GetDataBase()->GetIntegralEngine());
                    assert(ie_d);

                    const Vector<double>& ed(ie_d->itsExponents);
                    const Vector<double>& nd(ie_d->itsNormalizations);
                    int Nd=ie_d->itsN;
                    int Ld=ie_d->itsL;

                    ERIProxy cp(Coulomb     ,start_a,start_b,start_c,start_d);
                    ERIProxy ep(exchange    ,start_a,start_b,start_c,start_d);
                    ERIProxy tp(tracker_eris,start_a,start_b,start_c,start_d);
                    for (index_t ia=1; ia<=Na; ia++)
                        for (index_t ib=1; ib<=Nb; ib++)
                            for (index_t ic=1; ic<=Nc; ic++)
                                for (index_t id=1; id<=Nd; id++)
                                {
                                    double norm=na(ia)*nb(ib)*nc(ic)*nd(id);
                                    // Coulomb case
                                    double J=0.0,K=0.0;
                                    if (La==Lb && Lc==Ld) //TODO these ifs can go outside the abcd loops.
                                    {
                                        SlaterIntegrals R(ea(ia)+eb(ib),ec(ic)+ed(id));
                                        J=FourPi2*R(0,La,Lb,Lc,Ld)*norm;
                                    }
                                    // Exchange case
                                    if (La==Lc && Lb==Ld) //TODO these ifs can go outside the abcd loops.
                                    {
                                        SlaterIntegrals R(ea(ia)+eb(ib),ec(ic)+ed(id));
                                        K=FourPi2*R.DoExchangeSum(La,Lb,Lc,Ld)*norm;
                                    }

                                    if (tp(ia,ib,ic,id)!=-1.0)
                                    {
                                        tp(ia,ib,ic,id)=-1.0;
                                        cp(ia,ib,ic,id)=J;
                                        ep(ia,ib,ic,id)=K;
                                    }
                                    else if (fabs(cp(ia,ib,ic,id)-J)>1e-8)
                                    {
                                        std::cout << " Re-assign J " << ia << ","<< ib << ","<< ic << ","<< id << ", was" << ep(ia,ib,ic,id) << " new=" << J << std::endl;
                                    }
                                    else if (fabs(ep(ia,ib,ic,id)-K)>1e-8)
                                    {
                                        std::cout << " Re-assign K " << ia << ","<< ib << ","<< ic << ","<< id << ", was" << ep(ia,ib,ic,id) << " new=" << J << std::endl;
                                    }

                                }
                }
            }
        }
    }
}


//----------------------------------------------------------------------------------------
//
//  Special integrals
//
SMatrix<double> SphericalGaussianIE::MakeKinetic() const
{
    CheckBasisSet();
    SMatrix<double> ret(itsN,itsN);
    SMatrix<double>::Subscriptor s(ret);

    const Vector<double>& e(itsExponents);

    double a=TwoLPlusOne(itsL);
    for (index_t i=1; i<=itsN; i++)
        for (index_t j=i; j<=itsN; j++)
        {
            double t=e(i)+e(j);
            int L1=itsL+1;
            s(i,j)=0.5*a*
                   (
                       (L1*L1 + itsL*L1) * GaussianIntegral(t,2*itsL-2)
                       -2 *  L1 * (e(i)+e(j)) * GaussianIntegral(t,2*itsL  )
                       +4 * e(i)*e(j)         * GaussianIntegral(t,2*itsL+2)
                   );
        }

    Normalize(ret);
    return ret;
}

SMatrix<double> SphericalGaussianIE::MakeNuclear(const Cluster& theCluster) const
{
    CheckBasisSet();
    assert(theCluster.GetNumAtoms()==1);

    SMatrix<double> ret(itsN,itsN);
    SMatrix<double>::Subscriptor s(ret);

    const Vector<double>& e(itsExponents);

    double a=TwoLPlusOne(itsL);
    for (index_t i=1; i<=itsN; i++)
        for (index_t j=i; j<=itsN; j++)
            s(i,j)= a*GaussianIntegral(e(i)+e(j),2*itsL-1);

    ret *= (-(double)(theCluster.GetNuclearCharge()));
    Normalize(ret);
    return ret;
}




SMatrix<double> SphericalGaussianIE::MakeOverlap(const TBasisFunction<double>& bf) const
{
    CheckBasisSet();
    const SphericalGaussianBF* bfag= dynamic_cast<const SphericalGaussianBF*>(&bf);
    assert(bfag);

    SMatrix<double> ret(itsN,itsN);
    SMatrix<double>::Subscriptor s(ret);

    const Vector<double>& e(itsExponents);
    double ebf=bfag->itsExponent;
    int    lbf=bfag->itsL;
    assert(lbf==0);

    double a=TwoLPlusOne(itsL);
    for (index_t i=1; i<=itsN; i++)
        for (index_t j=i; j<=itsN; j++)
            s(i,j)=a*GaussianIntegral(e(i)+e(j)+ebf,2*itsL+lbf);

    Normalize(ret);
    ret*=bf.GetNormalization();
    return ret;
}

SMatrix<double> SphericalGaussianIE::MakeRepulsion(const TBasisFunction<double>& bf) const
{
    CheckBasisSet();
    const SphericalGaussianBF* bfag= dynamic_cast<const SphericalGaussianBF*>(&bf);
    assert(bfag);

    SMatrix<double> ret(itsN,itsN);
    SMatrix<double>::Subscriptor s(ret);

    const Vector<double>& e(itsExponents);
    double ebf=bfag->itsExponent;
    int    lbf=bfag->itsL;
    assert(lbf==0);

    double a=TwoLPlusOne(itsL);
    for (index_t i=1; i<=itsN; i++)
        for (index_t j=i; j<=itsN; j++)
        {
            SlaterIntegrals R(e(i)+e(j),ebf);
            s(i,j)=a*FourPi2*R(0,itsL,itsL,lbf,0);
//            s(i,j)=a*GaussianRepulsionIntegral(e(i)+e(j),ebf,2*itsL,lbf);

        }

    Normalize(ret);
    ret*=bf.GetNormalization();
    return ret;
}









