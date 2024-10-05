// File: PolarizedGaussianIE1.C  Here is where all the integral get calculated.


//#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianBF.H"
//#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianBS.H"
#include "BasisSetImplementation/PolarizedGaussian/PolarizedGaussianIE1.H"
#include "BasisSetImplementation/PolarizedGaussian/BasisFunctionBlock.H"
#include "BasisSetImplementation/PolarizedGaussian/ContractedGaussian/ContractedGaussianRF.H"
//#include "BasisSetImplementation/SphericalGaussian/SphericalSymmetryQN.H"
#include "Cluster.H"
#include "oml/matrix.h"
#include "oml/smatrix.h"
//#include "IntegralDataBase.H"
//#include "BasisSet.H"
//#include "Cluster.H"
//#include "Misc/MatrixList.H"
#include "Misc/ERIList.H"
#include "Misc/ERIProxy.H"
#include "Misc/MatrixList.H"
//#include <cassert>
//#include <iostream>
//#include <stdlib.h>

//-----------------------------------------------------------------
//
//  Construction zone.
//
PolarizedGaussianIE1::PolarizedGaussianIE1(const blocks_t& _blocks, const RVec& _ns)
    : blocks(_blocks)
    , ns(_ns)
    , ons(OuterProduct(ns))
{
    for (auto bl:blocks)
        for (auto p:bl->itsPols)
        {
            radials.push_back(bl->itsRadial);
            pols.push_back(p);
        }
};

void PolarizedGaussianIE1::Normalize(SMat& s) const
{
    s=DirectMultiply(s,ons);
}

//-----------------------------------------------------------------
//
//  Streamable Object stuff
//
IntegralEngine1<double>* PolarizedGaussianIE1::Clone() const
{
    return new PolarizedGaussianIE1(*this);
}
//----------------------------------------------------------------------------------------
//
//  Overlap type integrals
//
PolarizedGaussianIE1::SMat PolarizedGaussianIE1::MakeOverlap() const
{
    SMat s(size());
    Fill(s,0.0);
    for (auto a(blocks.begin()); a!=blocks.end(); a++)
        for (auto b(a); b!=blocks.end(); b++)
        {
            BasisFunctionBlockPair p(*a,*b);
            (*a)->itsRadial->Get2CenterIntegrals(RadialFunction::Overlap2C,p,s,NULL,1.0);
        }
    
    int N=size();
    SMat s1(size());
    for (index_t ia=0;ia<N;ia++)
        for (index_t ib=ia;ib<N;ib++)
        {
            s1(ia+1,ib+1)=radials[ia]->Integrate(RadialFunction::Overlap2C,radials[ib],pols[ia],pols[ib],cache);
            assert(fabs(s(ia+1,ib+1)-s1(ia+1,ib+1))<1e-14);
        }    
        
    cache.Report(std::cout);
    Normalize(s);
    return s;
}

PolarizedGaussianIE1::SMat PolarizedGaussianIE1::MakeRepulsion() const
{
    
    SMat s(size());
    Fill(s,0.0);
    for (auto a(blocks.begin()); a!=blocks.end(); a++)
        for (auto b(a); b!=blocks.end(); b++)
        {
            BasisFunctionBlockPair p(*a,*b);
            (*a)->itsRadial->Get2CenterIntegrals(RadialFunction::Repulsion2C,p,s,NULL,1.0);
        }
        
     int N=size();
    SMat s1(size());
    for (index_t ia=0;ia<N;ia++)
        for (index_t ib=ia;ib<N;ib++)
        {
            s1(ia+1,ib+1)=radials[ia]->Integrate(RadialFunction::Repulsion2C,radials[ib],pols[ia],pols[ib],cache);
            double err=fabs(s(ia+1,ib+1)-s1(ia+1,ib+1));
            //std::cout << "Repulsion error=" << err << std::endl;
            assert(err<1e-14);
        }    
        
    cache.Report(std::cout);

    Normalize(s);
    return s;
}
PolarizedGaussianIE1::SMat PolarizedGaussianIE1::MakeKinetic() const
{
    SMat s(size());
    Fill(s,0.0);
    for (auto a(blocks.begin()); a!=blocks.end(); a++)
        for (auto b(a); b!=blocks.end(); b++)
        {
            BasisFunctionBlockPair p(*a,*b);
            (*a)->itsRadial->Get2CenterIntegrals(RadialFunction::Kinetic,p,s,NULL,1.0);
        }
    
    int N=size();
    SMat s1(size());
    for (index_t ia=0;ia<N;ia++)
        for (index_t ib=ia;ib<N;ib++)
        {
            s1(ia+1,ib+1)=radials[ia]->Integrate(RadialFunction::Kinetic,radials[ib],pols[ia],pols[ib],cache);
            double err=fabs(s(ia+1,ib+1)-s1(ia+1,ib+1));
            //std::cout << "Repulsion error=" << err << std::endl;
            assert(err<1e-14);
        }    
        

    Normalize(s);
    return s;
}
//
PolarizedGaussianIE1::SMat PolarizedGaussianIE1::MakeNuclear(const Cluster& cl) const
{
    SMat s(size());
    Fill(s,0.0);
    for (auto a(blocks.begin()); a!=blocks.end(); a++)
        for (auto b(a); b!=blocks.end(); b++)
        {
            BasisFunctionBlockPair p(*a,*b);
            (*a)->itsRadial->Get2CenterIntegrals(RadialFunction::Nuclear,p,s,&cl,1.0);
        }
    
    int N=size();
    SMat s1(size());
    for (index_t ia=0;ia<N;ia++)
        for (index_t ib=ia;ib<N;ib++)
        {
            assert(&cl);
            s1(ia+1,ib+1)=radials[ia]->Integrate(RadialFunction::Nuclear,radials[ib],pols[ia],pols[ib],cache,&cl);
            double err=fabs(s(ia+1,ib+1)-s1(ia+1,ib+1));
            std::cout << "Nuclear error=" << err << std::endl;
            assert(err<1e-14);
        }    
        

    Normalize(s);
    return s;
}


PolarizedGaussianIE1::Mat PolarizedGaussianIE1::MakeOverlap(const IE* ie) const
{
    // No UT coverage.
    assert(false); 
    return Mat();
}

//
PolarizedGaussianIE1::RVec PolarizedGaussianIE1::MakeOverlap(const ScalarFunction<double>& f) const
{
    // No UT coverage.  Only used for numerical integrations.
    assert(false);
    return RVec();
}

//PolarizedGaussianIE1::SMat PolarizedGaussianIE1::MakeOverlap(const BasisFunctionBlock& c) const
//{    
//    SMat s(size());
//    for (auto a(bls.begin()); a!=bls.end(); a++)
//            s(i,j)=GaussianIntegral(es(i)+es(j)+eo,2*L+Lo)*ns(i)*ns(j)*no;
//    return s;
//}

void PolarizedGaussianIE1::MakeOverlap3C(MList& mlist, const IE* ie) const
{
    const PolarizedGaussianIE1* other=dynamic_cast<const PolarizedGaussianIE1*>(ie);;
    assert(other);
    assert(false);
//     mlist.Empty();
//     for (auto c:other->bls) mlist.Add(MakeOverlap(c))
//     mlist.Clear();
}
//
////----------------------------------------------------------------------------------------
////
////  Repulsion type integrals
////

//
PolarizedGaussianIE1::Mat PolarizedGaussianIE1::MakeRepulsion(const IE* ie) const
{
    assert(false);
    return SMat();

}
//
//
//
Vector<double> PolarizedGaussianIE1::MakeRepulsion(const ScalarFunction<double>& f) const
{
    assert(false);
    return RVec();
}

//
//PolarizedGaussianIE1::SMat PolarizedGaussianIE1::MakeRepulsion(const bf_tuple& bf) const
//{    
//    assert(false);
//    return SMat();
//}


void PolarizedGaussianIE1::MakeRepulsion3C(MList& mlist, const IE* ie) const
{
    const PolarizedGaussianIE1* other=dynamic_cast<const PolarizedGaussianIE1*>(ie);;
    assert(other);
    assert(false);
}
//
////
////  This is where we do the big double loop over basis sets.
////
void PolarizedGaussianIE1::MakeRepulsion4C(ERIList& Coulomb, ERIList& exchange, const iev_t& iev) const
{
    assert(false);
}

PolarizedGaussianIE1::jk_t PolarizedGaussianIE1::Make4C(const iev_t&) const
{
    ERIList1 J,K;
    return std::make_pair(J,K);
}

////
//
////----------------------------------------------------------------------------------------
////
////  Special integrals
////

PolarizedGaussianIE1::RVec PolarizedGaussianIE1::MakeNormalization() const
{
    return ns;
}

PolarizedGaussianIE1::RVec PolarizedGaussianIE1::MakeCharge() const
{
    assert(false);
    return RVec();
}

std::ostream& PolarizedGaussianIE1::Write(std::ostream& os) const
{
    return os;
}
std::istream& PolarizedGaussianIE1::Read (std::istream& is)
{
    return is;
}



