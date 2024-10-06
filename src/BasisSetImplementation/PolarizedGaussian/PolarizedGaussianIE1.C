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
#include <vector>

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
        
    //cache.Report(std::cout);
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
            if (err>0.0) std::cout << "Repulsion error=" << log10(err) << std::endl;
            assert(err<1e-14);
        }    
        
    //cache.Report(std::cout);

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
            //if (err>0.0) std::cout << "Kinetic error=" << log10(err) << std::endl;
            assert(err<1e-12);
        }    
        

    Normalize(s);
    return s;
}
//
PolarizedGaussianIE1::SMat PolarizedGaussianIE1::MakeNuclear(const Cluster& cl) const
{
//    SMat s(size());
//    Fill(s,0.0);
//    for (auto a(blocks.begin()); a!=blocks.end(); a++)
//        for (auto b(a); b!=blocks.end(); b++)
//        {
//            BasisFunctionBlockPair p(*a,*b);
//            (*a)->itsRadial->Get2CenterIntegrals(RadialFunction::Nuclear,p,s,&cl,1.0);
//        }
    
    int N=size();
    SMat s(size());
    for (index_t ia=0;ia<N;ia++)
        for (index_t ib=ia;ib<N;ib++)
        {
            assert(&cl);
            s(ia+1,ib+1)=radials[ia]->Integrate(RadialFunction::Nuclear,radials[ib],pols[ia],pols[ib],cache,&cl);
            //double err=fabs(s(ia+1,ib+1)-s1(ia+1,ib+1));
            //if (err>0.0) std::cout << "Nuclear error=" << log10(err) << std::endl;
            //assert(err<1e-12);
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


template <class M> std::vector<M> MakeMatrixList(int n, int N)
{
    std::vector<M> ret;
    for (index_t i=0; i<n; i++)
    {
        ret.push_back(M(N,N));
        Fill(ret.back(),0.0);
    }
    return ret;
}


void PolarizedGaussianIE1::MakeOverlap3C(MList& mlist, const IE* ie) const
{
    const PolarizedGaussianIE1* other=dynamic_cast<const PolarizedGaussianIE1*>(ie);;
    mlist.Empty();
    int Nc=other->size();
    for (index_t ic=0;ic<Nc;ic++)
    {
        SMat s=MakeOverlap3C(other->radials[ic],other->pols[ic]);
        s*=other->ns(ic+1);
        mlist.Add(s);
    }    
    
    mlist.Clear();
}

PolarizedGaussianIE1::SMat PolarizedGaussianIE1::MakeOverlap3C(const RadialFunction* rc, const Polarization& pc) const
{
    int N=size();
    SMat s(size());
    for (index_t ia=0;ia<N;ia++)
        for (index_t ib=ia;ib<N;ib++)
            s(ia+1,ib+1)=rc->Integrate(RadialFunction::Overlap3C,radials[ia],radials[ib],pols[ia],pols[ib],pc,cache);
        
    Normalize(s);
    return s;
}
//------------------------------------------------------------------------------------
//
//  Calculates repulsion matricies <ab|1/r12|c> for a block of basis functions c, over
//  the whole basis set a and b.
//  *** No normalization is done at this point***
//  **** Only elements above the diagonal will be defined, everything below the diagonal
//  will be zero ***.
//

void PolarizedGaussianIE1::MakeOverlap3C(const BasisFunctionBlock& c,std::vector<SMat >& ret) const
{
   for (auto a(blocks.begin()); a!=blocks.end(); a++)
        for (auto b(a); b!=blocks.end(); b++)
        {
            BasisFunctionBlockTriplet t(*a,*b,&c);
            (*a)->itsRadial->Get3CenterIntegrals(RadialFunction::Overlap3C,t,ret,1.0);
        }

}

//
////----------------------------------------------------------------------------------------
////
////  Repulsion type integrals
////

//
void PolarizedGaussianIE1::MakeRepulsion3C(MList& mlist,const IE* ie) const
{
    const PolarizedGaussianIE1* other=dynamic_cast<const PolarizedGaussianIE1*>(ie);;
    mlist.Empty();
    int Nc=other->size();
    for (index_t ic=0;ic<Nc;ic++)
    {
        SMat s=MakeRepulsion3C(other->radials[ic],other->pols[ic]);
        s*=other->ns(ic+1);
        mlist.Add(s);
    }    
    
    mlist.Clear();


}
PolarizedGaussianIE1::SMat PolarizedGaussianIE1::MakeRepulsion3C(const RadialFunction* rc, const Polarization& pc) const
{
    int N=size();
    SMat s(size());
    for (index_t ia=0;ia<N;ia++)
        for (index_t ib=ia;ib<N;ib++)
            s(ia+1,ib+1)=rc->Integrate(RadialFunction::Repulsion3C,radials[ia],radials[ib],pols[ia],pols[ib],pc,cache);
        
    Normalize(s);
    return s;
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
PolarizedGaussianIE1::Mat PolarizedGaussianIE1::MakeRepulsion(const IE*) const
{    
    assert(false);
    return SMat();
}



//
////
////  This is where we do the big double loop over basis sets.
////
void PolarizedGaussianIE1::MakeRepulsion4C(ERIList& Coulomb, ERIList& exchange, const iev_t& iev) const
{
    assert(false);
}

PolarizedGaussianIE1::jk_t PolarizedGaussianIE1::Make4C(const iev_t& ies) const
{
#ifdef DEBUG
    assert(ies.size()==1);  //Don't have multiple Irrep basis sets for Molecules yet.
    const IE* ie=ies[0];
    const PolarizedGaussianIE1* other=dynamic_cast<const PolarizedGaussianIE1*>(ie);
    assert(other);
    assert(other==this);   //Again don't have multiple Irrep basis sets for Molecules yet.
#endif
    
    int N=size();
    ERIList1 J(N,-1.0),K;
    std::cout << N << " " << J.itsData.size() <<" " << K.itsData.size() << std::endl;

    for (index_t ia:ns.indices())
        for (index_t ib:ns.indices(ia))
            for (index_t ic:ns.indices())
                for (index_t id:ns.indices(ic))
                {
                    if (J(ia,ib,ic,id)==-1.0)
                    {
                        //std::cout << "abcd=(" << ia << "," << ib << "," << ic << "," << id << ")" << std::endl;
                        double norm=ns(ia)*ns(ib)*ns(ic)*ns(id);
                        assert(radials[id-1]);
                        J(ia,ib,ic,id)=norm * radials[id-1]->Integrate(radials[ia-1],radials[ib-1],radials[ic-1],pols[ia-1],pols[ib-1],pols[ic-1],pols[id-1],cache);
                    }
                }
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



