// File: SphericalGaussianIE1.C  Here is where all the integral get calculated.


#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianIE1.H"
#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianBS.H" //Just to get IEClient
#include "BasisSetImplementation/SphericalGaussian/GaussianIntegrals.H"
#include "BasisSetImplementation/SphericalGaussian/SlaterIntegrals.H"
#include "Cluster.H"
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include "Misc/ERI4.H"

double SphericalGaussianIE1::FourPi2=4*4*Pi*Pi;

//-----------------------------------------------------------------
//
//  Construction zone.
//
SphericalGaussianIE1::SphericalGaussianIE1(size_t _L, const RVec& alphas)
    : L(_L)
    , es(alphas) //exponents
    , ns(es.size())
{
    for (auto i:es.indices())  ns(i)=GaussianNorm(es(i),L);
};


//-----------------------------------------------------------------
//
//  Streamable Object stuff
//
AnalyticIE<double>* SphericalGaussianIE1::Clone() const
{
    return new SphericalGaussianIE1(*this);
}
//----------------------------------------------------------------------------------------
//
//  Overlap type integrals
//
SphericalGaussianIE1::SMat SphericalGaussianIE1::MakeOverlap() const
{
    size_t N=es.size();
    SMat s(N);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
            s(i,j)=GaussianIntegral(es(i)+es(j),2*L)*ns(i)*ns(j);

    return s;
}


SphericalGaussianIE1::SMat SphericalGaussianIE1::MakeOverlap(const bf_tuple& bf) const
{    
    size_t N=size();
    int Lo;
    double eo,no;
    std::tie(Lo,eo,no)=bf;
    SMat s(N);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
            s(i,j)=GaussianIntegral(es(i)+es(j)+eo,2*L+Lo)*ns(i)*ns(j)*no;
    return s;
}

SphericalGaussianIE1::ERI3 SphericalGaussianIE1::MakeOverlap3C(const IE* ie) const
{
    const SphericalGaussianIE1* other=dynamic_cast<const SphericalGaussianIE1*>(ie);;
    assert(other);

    ERI3 s3;
    for (auto i:other->es.indices()) s3.push_back(MakeOverlap((*other)(i)));
    return s3;
}

//----------------------------------------------------------------------------------------
//
//  Repulsion type integrals
//
SphericalGaussianIE1::SMat SphericalGaussianIE1::MakeRepulsion() const
{
    size_t N=size();
    SMat r(N,N);
    for (auto i:r.rows())
        for (auto j:r.cols(i))
            r(i,j)=GaussianRepulsionIntegral(es(i),es(j),L,L)*ns(i)*ns(j);

    return r;
}
//
SphericalGaussianIE1::Mat SphericalGaussianIE1::MakeRepulsion(const IE* iea,const IE* ieb) const
{
    const SphericalGaussianIE1* a=dynamic_cast<const SphericalGaussianIE1*>(iea);;
    assert(a);
    const SphericalGaussianIE1* b=dynamic_cast<const SphericalGaussianIE1*>(ieb);;
    assert(b);
    size_t Na=a->es.size(), Nb=b->es.size();
    Mat s(Na,Nb);
    for (auto i:s.rows())
        for (auto j:s.cols())
            s(i,j)=GaussianRepulsionIntegral(a->es(i),b->es(j),a->L,b->L)*a->ns(i)*b->ns(j);

    return s;
}

SphericalGaussianIE1::Mat SphericalGaussianIE1::MakeRepulsion(iec_t* iea,iec_t* ieb) const
{
    const SphericalGaussianIEClient* a=dynamic_cast<const SphericalGaussianIEClient*>(iea);;
    assert(a);
    const SphericalGaussianIEClient* b=dynamic_cast<const SphericalGaussianIEClient*>(ieb);;
    assert(b);
    size_t Na=a->es.size(), Nb=b->es.size();
    Mat s(Na,Nb);
    for (auto i:s.rows())
        for (auto j:s.cols())
            s(i,j)=GaussianRepulsionIntegral(a->es(i),b->es(j),a->Ls(i),b->Ls(j))*a->ns(i)*b->ns(j);

    return s;
}

//
SphericalGaussianIE1::SMat SphericalGaussianIE1::MakeRepulsion(const bf_tuple& bf) const
{    
    size_t N=es.size();
    int Lo;
    double eo,no;
    std::tie(Lo,eo,no)=bf;
//    int Lo=std::get<1>(bf);
//    double eo=std::get<2>(bf),no=std::get<3>(bf);
    SMat s(N,N);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
        {
            SlaterIntegrals R(es(i)+es(j),eo);
            s(i,j)=FourPi2*R(0,L,L,Lo,0)*ns(i)*ns(j)*no;
        }
    return s;
}


SphericalGaussianIE1::ERI3 SphericalGaussianIE1::MakeRepulsion3C(const IE* ie) const
{
    const SphericalGaussianIE1* other=dynamic_cast<const SphericalGaussianIE1*>(ie);;
    assert(other);

    ERI3 s3;
    for (auto i:other->es.indices()) s3.push_back(MakeRepulsion((*other)(i)));
    return s3;
}


SphericalGaussianIE1::SGparams::SGparams(const iev_t& iev)
{
    size_t N=0;
    for (auto ia: iev) N+=ia->size();
    Ls.SetLimits(N);
    es.SetLimits(N);
    ns.SetLimits(N);
    index_t i=1;
    for (auto ia: iev)
    {
        const SphericalGaussianIE1* sg=dynamic_cast<const SphericalGaussianIE1*>(ia);
        assert(sg);
        for (auto i1:sg->es.indices())
        {
            Ls(i)=sg->L;
            es(i)=sg->es(i1);
            ns(i)=sg->ns(i1);
            i++;
        }
    }
}

SphericalGaussianIE1::jk_t SphericalGaussianIE1::Make4C(const iev_t& iev) const
{
    SphericalGaussianIE1::SGparams sg(iev);
    size_t N=sg.size();
    ERI4 J(N,-1.0),K(N,-1.0);
    std::cout << N << " " << J.itsData.size() <<" " << K.itsData.size() << std::endl;

    for (index_t ia:sg.es.indices())
        for (index_t ib:sg.es.indices(ia))
            for (index_t ic:sg.es.indices())
                for (index_t id:sg.es.indices(ic))
                {
                    bool doJ = sg.Ls(ia)==sg.Ls(ib) && sg.Ls(ic)==sg.Ls(id) && J(ia,ib,ic,id)==-1.0;
                    bool doK = sg.Ls(ia)==sg.Ls(ic) && sg.Ls(ib)==sg.Ls(id) && K(ia,ib,ic,id)==-1.0;
                    if (doJ || doK)
                    {
                        double norm=sg.ns(ia)*sg.ns(ib)*sg.ns(ic)*sg.ns(id);
                        SlaterIntegrals R(sg.es(ia)+sg.es(ib),sg.es(ic)+sg.es(id));
                        if (doJ)
                            J(ia,ib,ic,id)=FourPi2*R(0,sg.Ls(ia),sg.Ls(ib),sg.Ls(ic),sg.Ls(id))*norm;
                        if (doK)
                            K(ia,ib,ic,id)=FourPi2*R.DoExchangeSum(sg.Ls(ia),sg.Ls(ib),sg.Ls(ic),sg.Ls(id))*norm;
                        else
                            K(ia,ib,ic,id)=0.0;
//                        std::cout << "L=(" << sg.Ls(ia) << "," << sg.Ls(ib) << "," << sg.Ls(ic) << "," << sg.Ls(id) 
//                        << ") abcd=(" << ia << "," << ib << "," << ic << "," << id << ")  J=" << J(ia,ib,ic,id) << std::endl;
                                
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
SphericalGaussianIE1::SMat SphericalGaussianIE1::MakeKinetic() const
{
    size_t N=size();
    SMatrix<double> Hk(N);
    for (auto i:Hk.rows())
        for (auto j:Hk.cols(i))
        {
            double t=es(i)+es(j);
            int L1=L+1;
            Hk(i,j)=0.5*ns(i)*ns(j)*
                   (
                       (L1*L1 + L*L1) * GaussianIntegral(t,2*L-2)
                       -2*L1 * t      * GaussianIntegral(t,2*L  )
                       +4*es(i)*es(j) * GaussianIntegral(t,2*L+2)
                   );
        }

    return Hk;
}
//
SphericalGaussianIE1::SMat SphericalGaussianIE1::MakeNuclear(const Cluster& cl) const
{
    size_t N=size();
    SMatrix<double> Hn(N);
    double Z=-cl.GetNuclearCharge();
    for (auto i:Hn.rows())
        for (auto j:Hn.cols(i))
            Hn(i,j)= Z*GaussianIntegral(es(i)+es(j),2*L-1)*ns(i)*ns(j);

    return Hn;
}

SphericalGaussianIE1::RVec SphericalGaussianIE1::MakeNormalization() const
{
    return ns;
}

SphericalGaussianIE1::RVec SphericalGaussianIE1::MakeCharge() const
{
    RVec c(es.size());
    for (auto i:es.indices())  c(i)=GaussianIntegral(es(i),L)*ns(i);
    return c;
}

std::ostream& SphericalGaussianIE1::Write(std::ostream& os) const
{
    return os << L << " " << es << ns;
}
std::istream& SphericalGaussianIE1::Read (std::istream& is)
{
    return is >> L >> es >> ns;
}



