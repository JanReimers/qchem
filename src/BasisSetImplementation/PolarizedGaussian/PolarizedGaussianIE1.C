// File: PolarizedGaussianIE1.C  Here is where all the integral get calculated.


#include "BasisSetImplementation/PolarizedGaussian/PolarizedGaussianIE1.H"
#include "BasisSetImplementation/PolarizedGaussian/BasisFunctionBlock.H"
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include "Misc/ERIList.H"
#include "Misc/ERIProxy.H"
#include "Misc/MatrixList.H"

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
AnalyticIE<double>* PolarizedGaussianIE1::Clone() const
{
    return new PolarizedGaussianIE1(*this);
}
//----------------------------------------------------------------------------------------
//
//  Overlap type integrals
//
PolarizedGaussianIE1::SMat PolarizedGaussianIE1::MakeOverlap() const
{
   int N=size();
    SMat s(N);
    for (index_t ia=0;ia<N;ia++)
        for (index_t ib=ia;ib<N;ib++)
            s(ia+1,ib+1)=radials[ia]->Integrate(RadialFunction::Overlap2C,radials[ib],pols[ia],pols[ib],cache);

    Normalize(s);
    return s;
}

PolarizedGaussianIE1::SMat PolarizedGaussianIE1::MakeRepulsion() const
{
    int N=size();
    SMat s(N);
    for (index_t ia=0;ia<N;ia++)
        for (index_t ib=ia;ib<N;ib++)
            s(ia+1,ib+1)=radials[ia]->Integrate(RadialFunction::Repulsion2C,radials[ib],pols[ia],pols[ib],cache);
        
    Normalize(s);
    return s;
}

PolarizedGaussianIE1::SMat PolarizedGaussianIE1::MakeKinetic() const
{
    int N=size();
    SMat s(N);
    for (index_t ia=0;ia<N;ia++)
        for (index_t ib=ia;ib<N;ib++)
            s(ia+1,ib+1)=radials[ia]->Integrate(RadialFunction::Kinetic,radials[ib],pols[ia],pols[ib],cache);

    Normalize(s);
    return s;
}
//
PolarizedGaussianIE1::SMat PolarizedGaussianIE1::MakeNuclear(const Cluster& cl) const
{
    assert(&cl);
    int N=size();
    SMat s(N);
    for (index_t ia=0;ia<N;ia++)
        for (index_t ib=ia;ib<N;ib++)
            s(ia+1,ib+1)=radials[ia]->Integrate(RadialFunction::Nuclear,radials[ib],pols[ia],pols[ib],cache,&cl);

    Normalize(s);
    return s;
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
    SMat s(N);
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
    SMat s(N);
    for (index_t ia=0;ia<N;ia++)
        for (index_t ib=ia;ib<N;ib++)
            s(ia+1,ib+1)=rc->Integrate(RadialFunction::Repulsion3C,radials[ia],radials[ib],pols[ia],pols[ib],pc,cache);
        
    Normalize(s);
    return s;
}

//
PolarizedGaussianIE1::Mat PolarizedGaussianIE1::MakeRepulsion(const IE* ieb) const
{    
    const PolarizedGaussianIE1* pgb=dynamic_cast<const PolarizedGaussianIE1*>(ieb);
    assert(pgb);
    int Na=size(),Nb=pgb->size();
    Mat s(Na,Nb);
    for (index_t ia=0;ia<Na;ia++)
        for (index_t ib=0;ib<Nb;ib++)
            s(ia+1,ib+1)=radials[ia]->Integrate(RadialFunction::Repulsion2C,
                pgb->radials[ib],pols[ia],pgb->pols[ib],cache)*ns(ia+1)*pgb->ns(ib+1);
    assert(!isnan(s));
    return s;
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
    RVec c(size());
    int i=0;
    for (auto r:radials)
    {
        c(i+1)=r->GetCharge(pols[i])*ns(i+1); 
        i++;       
    }
    assert(!isnan(c));
    return c;
}

std::ostream& PolarizedGaussianIE1::Write(std::ostream& os) const
{
    return os;
}
std::istream& PolarizedGaussianIE1::Read (std::istream& is)
{
    return is;
}



