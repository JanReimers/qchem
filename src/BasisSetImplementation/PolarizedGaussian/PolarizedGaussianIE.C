// File: PolarizedGaussianIE.C  Here is where all the integral get calculated.


#include "BasisSet.H"
#include "IntegralDataBase.H"
//#include "BasisSetImplementation/PolarizedGaussian/PolarizedGaussianIE.H"
#include "BasisSetImplementation/PolarizedGaussian/PolarizedGaussianBF.H"
#include "BasisSetImplementation/PolarizedGaussian/PolarizedGaussianBS.H"
#include "BasisSetImplementation/PolarizedGaussian/BasisFunctionBlock.H"
#include "Misc/Polarization.H"
#include "BasisSetImplementation/PolarizedGaussian/RadialFunction.H"
#include "Misc/MatrixList.H"
#include "Misc/ERIList.H"
#include "BasisSetImplementation/PolarizedGaussian/Gaussian/GaussianRF.H"

#include <cassert>

//-----------------------------------------------------------------
//
//  Construction zone.
//

PolarizedGaussianIE::PolarizedGaussianIE()
    : IntegralEngineImplementation<double>()
{};

PolarizedGaussianIE::PolarizedGaussianIE(const PolarizedGaussianIE& pgie)
    : IntegralEngineImplementation<double>(pgie)
{};

//-----------------------------------------------------------------
//
//  Streamable Object stuff
//
IntegralEngine<double>* PolarizedGaussianIE::Clone() const
{
    return new PolarizedGaussianIE(*this);
}

typedef optr_vector1<BasisFunctionBlock*>::const_iterator CITER;
//----------------------------------------------------------------------------------------
//
//  Overlap type integrals
//
PolarizedGaussianIE::SMat PolarizedGaussianIE::MakeOverlap() const
{
    CheckBasisSet();
    SMat ret(itsN,itsN);
    Fill(ret,0.0);
    auto blocks=GetBlocks(*itsBasisSet);
    for (CITER a(blocks.begin()); a!=blocks.end(); a++)
        for (CITER b(a); b!=blocks.end(); b++)
        {
            BasisFunctionBlockPair p(*a,*b);
            (*a)->itsRadial->Get2CenterIntegrals(RadialFunction::Overlap2C,p,ret,NULL,1.0);
        }


    Normalize(ret);
    return ret;
}

PolarizedGaussianIE::Mat PolarizedGaussianIE::MakeOverlap(const TBasisSet<double>& theBasisSet) const
{
    // No UT coverage
    CheckBasisSet();
    const PolarizedGaussianIE*    OtherIE;
    const IntegralEngine<double>* theOtherIE;

    theOtherIE=theBasisSet.GetDataBase()->GetIntegralEngine();
    OtherIE  = dynamic_cast<const PolarizedGaussianIE*>(theOtherIE);
    assert(OtherIE);


    Mat ret(itsN,OtherIE->itsN);
    Fill(ret,0.0);
    const optr_vector1<BasisFunctionBlock*>& blocksa=GetBlocks(*itsBasisSet);
    const optr_vector1<BasisFunctionBlock*>& blocksb=GetBlocks(theBasisSet);
    for (CITER a(blocksa.begin()); a!=blocksa.end(); a++)
        for (CITER b(blocksb.begin()); b!=blocksb.end(); b++)
        {
            BasisFunctionBlockPair p(*a,*b);
            (*a)->itsRadial->Get2CenterIntegrals(RadialFunction::Overlap2C,p,ret,1.0);
        }

    Mat::Subscriptor      rs(ret);
    for (unsigned int i=1; i<=itsN; i++)
        for (unsigned int j=1; j<=OtherIE->itsN; j++)
            rs(i,j)*=itsNormalizations(i)*OtherIE->itsNormalizations(j);

    return ret;
}

PolarizedGaussianIE::Vec PolarizedGaussianIE::MakeOverlap(const ScalarFunction<double>& f) const
{
    // No UT coverage
    CheckBasisSet();
    const PolarizedGaussianBF* bfpolg= dynamic_cast<const PolarizedGaussianBF*>(&f);
    assert(bfpolg);

    SMat mat(itsN,itsN);
    Fill(mat,0.0);
    BasisFunctionBlock block(bfpolg->itsRadial->Clone(),1);
    block.Add(bfpolg->itsPol);

    const optr_vector1<BasisFunctionBlock*>& blocks=GetBlocks(*itsBasisSet);
    for (CITER a(blocks.begin()); a!=blocks.end(); a++)
    {
        BasisFunctionBlockPair p(*a,&block);
        (*a)->itsRadial->Get2CenterIntegrals(RadialFunction::Overlap2C,p,mat,NULL,1.0);
    }

    Vec ret=mat.GetRow(1);
    Normalize(ret);
    ret *= bfpolg->GetNormalization();
    return ret;
}

//------------------------------------------------------------------------------------
//
//  Calculates overlap matricies <ab|c> for a whole set basis functions c.
//  Everything gets normalized and the elements below the diagonal are defined.
//
void PolarizedGaussianIE::MakeOverlap3C(MList& ret, const TBasisSet<double>& otherBS) const
{
    ret.Empty();
    auto c=otherBS.begin();
    const optr_vector1<BasisFunctionBlock*>& blocks=GetBlocks(otherBS);
    for (CITER block(blocks.begin()); block!=blocks.end(); block++)
    {
        std::vector<SMat > list=MakeMatrixList((*block)->size());
        MakeOverlap3C(**block,list);
       for(std::vector<SMat >::iterator blist(list.begin()); blist!=list.end(); blist++,c++)
        {
            SMat& m=*blist;
            const PolarizedGaussianBF* bfpolg= dynamic_cast<const PolarizedGaussianBF*>(*c);
            assert(bfpolg);
            m *= bfpolg->GetNormalization();
            Normalize(m);
            ret.Add(m);
        }
    }

    ret.Clear();
}


void PolarizedGaussianIE::MakeOverlap3C(const BasisFunctionBlock& c,std::vector<SMat >& ret) const
{
    CheckBasisSet();

    const optr_vector1<BasisFunctionBlock*>& blocks=GetBlocks(*itsBasisSet);
    for (CITER a(blocks.begin()); a!=blocks.end(); a++)
        for (CITER b(a); b!=blocks.end(); b++)
        {
            BasisFunctionBlockTriplet t(*a,*b,&c);
            (*a)->itsRadial->Get3CenterIntegrals(RadialFunction::Overlap3C,t,ret,1.0);
        }
}




//----------------------------------------------------------------------------------------
//
//  Repulsion type integrals
//
PolarizedGaussianIE::SMat PolarizedGaussianIE::MakeRepulsion() const
{
    CheckBasisSet();
    SMat ret(itsN,itsN);
    Fill(ret,0.0);
    const optr_vector1<BasisFunctionBlock*>& blocks=GetBlocks(*itsBasisSet);
    for (CITER a(blocks.begin()); a!=blocks.end(); a++)
        for (CITER b(a); b!=blocks.end(); b++)
        {
            BasisFunctionBlockPair p(*a,*b);
            (*a)->itsRadial->Get2CenterIntegrals(RadialFunction::Repulsion2C,p,ret,NULL,1.0);
        }

    Normalize(ret);
    return ret;
}

PolarizedGaussianIE::Mat PolarizedGaussianIE::MakeRepulsion(const TBasisSet<double>& theBasisSet) const
{
    CheckBasisSet();
    const PolarizedGaussianIE*    OtherIE;
    const IntegralEngine<double>* theOtherIE;

    theOtherIE=theBasisSet.GetDataBase()->GetIntegralEngine();
    OtherIE  = dynamic_cast<const PolarizedGaussianIE*>(theOtherIE);
    assert(OtherIE);

    Mat ret(itsN,OtherIE->itsN);
    Fill(ret,0.0);
    const optr_vector1<BasisFunctionBlock*>& blocksa=GetBlocks(*itsBasisSet);
    const optr_vector1<BasisFunctionBlock*>& blocksb=GetBlocks(theBasisSet);
    for (CITER a(blocksa.begin()); a!=blocksa.end(); a++)
        for (CITER b(blocksb.begin()); b!=blocksb.end(); b++)
        {
            BasisFunctionBlockPair p(*a,*b);
            (*a)->itsRadial->Get2CenterIntegrals(RadialFunction::Repulsion2C,p,ret,1.0);
        }


    Mat            ::Subscriptor      rs(ret);
    for (unsigned int i=1; i<=itsN; i++)
        for (unsigned int j=1; j<=OtherIE->itsN; j++)
            rs(i,j)*=itsNormalizations(i)*OtherIE->itsNormalizations(j);

    return ret;
}

PolarizedGaussianIE::Vec PolarizedGaussianIE::MakeRepulsion(const ScalarFunction<double>& f) const
{
    CheckBasisSet();
    const PolarizedGaussianBF* bfpolg= dynamic_cast<const PolarizedGaussianBF*>(&f);
    assert(bfpolg);

    SMat mat(itsN,itsN);
    Fill(mat,0.0);
    BasisFunctionBlock block(bfpolg->itsRadial->Clone(),1);
    block.Add(bfpolg->itsPol);

    const optr_vector1<BasisFunctionBlock*>& blocks=GetBlocks(*itsBasisSet);
    for (CITER a(blocks.begin()); a!=blocks.end(); a++)
    {
        BasisFunctionBlockPair p(*a,&block);
        (*a)->itsRadial->Get2CenterIntegrals(RadialFunction::Repulsion2C,p,mat,NULL,1.0);
    }

    Vec ret=mat.GetRow(1);
    Normalize(ret);
    ret *= bfpolg->GetNormalization();
    return ret;
}


//------------------------------------------------------------------------------------
//
//  Calculates repulsion matricies <ab|1/r12|c> for a whole set basis functions c.
//  Everything gets normalized and the elements below the diagonal are defined.
//
#include "BasisSetImplementation/SphericalGaussian/GaussianIntegrals.H"
void PolarizedGaussianIE::MakeRepulsion3C(MList& ret, const TBasisSet<double>& otherBS) const
{
    ret.Empty();
    auto c=otherBS.begin();
    const optr_vector1<BasisFunctionBlock*>& blocks=GetBlocks(otherBS);
    for (CITER block(blocks.begin()); block!=blocks.end(); block++)
    {
        std::vector<SMat > list=MakeMatrixList((*block)->size());
        MakeRepulsion3C(**block,list);

        for(std::vector<SMat >::iterator blist(list.begin()); blist!=list.end(); blist++,c++)
        {
            SMat& m=*blist;
            const PolarizedGaussianBF* bfpolg= dynamic_cast<const PolarizedGaussianBF*>(*c);
            assert(bfpolg);
            m *= bfpolg->GetNormalization();
            Normalize(m);
            ret.Add(m);
        }
    }
    ret.Clear();
}



//------------------------------------------------------------------------------------
//
//  Calculates repulsion matricies <ab|1/r12|c> for a block of basis functions c, over
//  the whole basis set a and b.
//  *** No normalization is done at this point***
//  **** Only elements above the diagonal will be defined, everything below the diagonal
//  will be zero ***.
//

void PolarizedGaussianIE::MakeRepulsion3C(const BasisFunctionBlock& c,std::vector<SMat >& ret) const
{
    CheckBasisSet();
    const optr_vector1<BasisFunctionBlock*>& blocks=GetBlocks(*itsBasisSet);
    for (CITER a(blocks.begin()); a!=blocks.end(); a++)
        for (CITER b(a); b!=blocks.end(); b++)
        {
            BasisFunctionBlockTriplet t(*a,*b,&c);
            (*a)->itsRadial->Get3CenterIntegrals(RadialFunction::Repulsion3C,t,ret,1.0);
        }

}

//------------------------------------------------------------------------------------
//
//  Calculates repulsion integrals <ab|1/r12|cd> over the whole internal basis set.
//  *** No normalization is done at this point***
//  **** Only elements above the diagonal will be defined, everything below the diagonal
//  will be zero ***.
//
void PolarizedGaussianIE::MakeRepulsion4C(ERIList& coulomb,  ERIList& exchange, const BasisGroup* bg) const
{
    CheckBasisSet();
    assert(bg->GetNumBasisSets()==1); //TODO Lets worry about support multiple groups later.
    coulomb.SetSize(bg->GetNumFunctions());

    const BasisSet* bs1=*bg->begin(); //Get the first basis set.
    const optr_vector1<BasisFunctionBlock*>& blocks=GetBlocks(*bs1);
    for (CITER a(blocks.begin()); a!=blocks.end(); a++)
        for (CITER b(a); b!=blocks.end(); b++)
        {
            for (CITER c(a); c!=blocks.end(); c++)
                for (CITER d(c); d!=blocks.end(); d++)
                {
                    BasisFunctionBlockQuartet q(*a,*b,*c,*d);
                    (*a)->itsRadial->GetRepulsion4C(q,coulomb,1.0);
                }
        }

    // Now normalize everything.
    const TBasisSet<double>* tbs=dynamic_cast<const TBasisSet<double>*>(bs1);
    Vector<double> Norm=tbs->GetDataBase()->GetNormalization();
    int N=bs1->GetNumFunctions();
    for (int ia=1; ia<=N; ia++)
        for (int ib=ia; ib<=N; ib++)
            for (int ic=ia; ic<=N; ic++)
                for (int id= ia==ic ? ib : ic; id<=N; id++)
                {
                    coulomb(ia,ib,ic,id)*=Norm(ia)*Norm(ib)*Norm(ic)*Norm(id);
                }

    // Now verify
#if 0
    typedef List<Polarization>::const_iterator PolITER;
    for (CITER a(blocks.begin()); a!=blocks.end(); a++)
        for (CITER b(a); b!=blocks.end(); b++)
        {
            for (CITER c(a); c!=blocks.end(); c++)
                for (CITER d(c); d!=blocks.end(); d++)
                {
                    double ERItotal=0.0;
                    int nPa=a->itsN;
                    for (PolITER bPa=a->itsPols.begin();  bPa!=a->itsPols.end();  bPa++,nPa++)
                    {
                        int nPb=b->itsN;
                        for (PolITER bPb=b->itsPols.begin();  bPb!=b->itsPols.end();  bPb++,nPb++)
                        {
                            int nPc=c->itsN;
                            for (PolITER bPc=c->itsPols.begin();  bPc!=c->itsPols.end();  bPc++,nPc++)
                            {
                                int nPd=d->itsN;
                                for (PolITER bPd=d->itsPols.begin();  bPd!=d->itsPols.end();  bPd++,nPd++)
                                {
                                    ERItotal+=coulomb(nPa,nPb,nPc,nPd);
                                }
                            }
                        }
                    }
                    GaussianRF* ga=dynamic_cast<GaussianRF*>(a->itsRadial);
                    GaussianRF* gb=dynamic_cast<GaussianRF*>(b->itsRadial);
                    GaussianRF* gc=dynamic_cast<GaussianRF*>(c->itsRadial);
                    GaussianRF* gd=dynamic_cast<GaussianRF*>(d->itsRadial);
                    if (ga && gb && gc && gd)
                    {
                        double ea=ga->itsExponent;
                        double eb=gb->itsExponent;
                        double ec=gc->itsExponent;
                        double ed=gd->itsExponent;
                        int la=a->itsPols.begin()->GetTotalL();
                        int lb=b->itsPols.begin()->GetTotalL();
                        int lc=c->itsPols.begin()->GetTotalL();
                        int ld=d->itsPols.begin()->GetTotalL();
                        double norm=GaussianNorm(ea,la)*GaussianNorm(eb,lb)*GaussianNorm(ec,lc)*GaussianNorm(ed,ld);
                        StreamableObject::SetToPretty();
//                        std::cout << "Checking <" << *bPa << *bPb << "||" <<  *bPc <<  *bPd << ">" << std::endl;
//                        double sg_result=GaussianRepulsionIntegral(ea,eb,ec,ed,la,lb,lc,ld)*norm*sqrt(2*la+1)*sqrt(2*lb+1)*sqrt(2*lc+1)*sqrt(2*ld+1);
                        double sg_result=GaussianRepulsionIntegral(ea,eb,ec,ed,la,lb,lc,ld)*norm;
                        if (la==1 && lb==1 && lc==1 && ld==1)
                            sg_result+=GaussianExchangeIntegral(ea,eb,ec,ed,la,lb,lc,ld)*norm;
                        if (la==1 && lb==0 && lc==1 && ld==0)
                            sg_result=GaussianExchangeIntegral(ea,eb,ec,ed,la,lb,lc,ld)*norm;
                        if (la==0 && lb==1 && lc==0 && ld==1)
                            sg_result=GaussianExchangeIntegral(ea,eb,ec,ed,la,lb,lc,ld)*norm;

                        sg_result*=sqrt(2*la+1)*sqrt(2*lb+1)*sqrt(2*lc+1)*sqrt(2*ld+1);
//                        if (bPa->n==0 && bPb->n==1 && bPc->n==0 && bPd->n==1)
//                        {
//                           std::cout << "expected " <<  eris(nPa,nPb,nPc,nPd) << ", got " << sg_result << std::endl;
//                        }
                        if (fabs(sg_result-ERItotal)>1e-8)
                        {
                            std::cout << "Integral missmatch <" << a->itsN << b->itsN << "||" <<  c->itsN <<  d->itsN << ">, expected "
                                      <<  ERItotal << ", got " << sg_result << std::endl;

                        }
                        else
                        {
                            std::cout << "Integral confirmed <" << a->itsN << b->itsN << "||" <<  c->itsN <<  d->itsN << ">, expected " <<  ERItotal << std::endl;
                        }

                    }

                }
        }
#endif



//    eris.Dump(std::cout);

}



//----------------------------------------------------------------------------------------
//
//  Special integrals
//
PolarizedGaussianIE::SMat PolarizedGaussianIE::MakeKinetic() const
{
    CheckBasisSet();
    SMat ret(itsN,itsN);
    Fill(ret,0.0);
    const optr_vector1<BasisFunctionBlock*>& blocks=GetBlocks(*itsBasisSet);
    for (CITER a(blocks.begin()); a!=blocks.end(); a++)
        for (CITER b(a); b!=blocks.end(); b++)
        {
            BasisFunctionBlockPair p(*a,*b);
            (*a)->itsRadial->Get2CenterIntegrals(RadialFunction::Kinetic,p,ret,NULL,1.0);
        }

    Normalize(ret);
    return ret;
}

PolarizedGaussianIE::SMat PolarizedGaussianIE::MakeNuclear(const Cluster& theCluster) const
{
    CheckBasisSet();
    SMat ret(itsN,itsN);
    Fill(ret,0.0);
    const optr_vector1<BasisFunctionBlock*>& blocks=GetBlocks(*itsBasisSet);
    for (CITER a(blocks.begin()); a!=blocks.end(); a++)
        for (CITER b(a); b!=blocks.end(); b++)
        {
            BasisFunctionBlockPair p(*a,*b);
            (*a)->itsRadial->Get2CenterIntegrals(RadialFunction::Nuclear,p,ret,&theCluster,1.0);
        }

    Normalize(ret);
    return ret;
}

//--------------------------------------------------------
//
//  Private utilities.
//
const optr_vector1<BasisFunctionBlock*>& PolarizedGaussianIE::GetBlocks(const BasisSet& bs) const
{
    const PolarizedGaussianBS* bspolg= dynamic_cast<const PolarizedGaussianBS*>(&bs);
    assert(bspolg);
    return bspolg->itsBlocks;
}




