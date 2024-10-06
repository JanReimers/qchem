// File: SphericalGaussianIE1.C  Here is where all the integral get calculated.


#include "BasisSetImplementation/SphericalGaussian/SphericalGaussianIE1.H"
#include "BasisSetImplementation/SphericalGaussian/GaussianIntegrals.H"
#include "BasisSetImplementation/SphericalGaussian/SlaterIntegrals.H"
#include "Cluster.H"
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include "Misc/ERIList.H"
#include "Misc/ERIProxy.H"
#include "Misc/MatrixList.H"

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
IntegralEngine1<double>* SphericalGaussianIE1::Clone() const
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

SphericalGaussianIE1::Mat SphericalGaussianIE1::MakeOverlap(const IE* ie) const
{
    assert(false);
    // No UT coverage.
    const SphericalGaussianIE1* other=dynamic_cast<const SphericalGaussianIE1*>(ie);;
    assert(other);
    size_t N=es.size(), No=other->es.size(), Lo=other->L;
    Mat s(N,No);
    for (auto i:s.rows())
        for (auto j:s.cols())
            s(i,j)=GaussianIntegral(es(i)+other->es(j),L+Lo)*ns(i)*other->ns(j);

    return s;
}

//
SphericalGaussianIE1::RVec SphericalGaussianIE1::MakeOverlap(const ScalarFunction<double>& f) const
{
    // No UT coverage.  Only used for numerical integrations.
    assert(false);
    return RVec();
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

void SphericalGaussianIE1::MakeOverlap3C(MList& mlist, const IE* ie) const
{
    const SphericalGaussianIE1* other=dynamic_cast<const SphericalGaussianIE1*>(ie);;
    assert(other);

    mlist.Empty();
    for (auto i:other->es.indices()) mlist.Add(MakeOverlap((*other)(i)));
    mlist.Clear();
}
//
////----------------------------------------------------------------------------------------
////
////  Repulsion type integrals
////
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
SphericalGaussianIE1::Mat SphericalGaussianIE1::MakeRepulsion(const IE* ie) const
{
    const SphericalGaussianIE1* other=dynamic_cast<const SphericalGaussianIE1*>(ie);;
    assert(other);
    size_t N=es.size(), No=other->es.size(), Lo=other->L;
    Mat s(N,No);
    for (auto i:s.rows())
        for (auto j:s.cols())
            s(i,j)=GaussianRepulsionIntegral(es(i),other->es(j),L,Lo)*ns(i)*other->ns(j);

    return s;
}
//
//
//
Vector<double> SphericalGaussianIE1::MakeRepulsion(const ScalarFunction<double>& f) const
{
    assert(false);
    return RVec();
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


void SphericalGaussianIE1::MakeRepulsion3C(MList& mlist, const IE* ie) const
{
    const SphericalGaussianIE1* other=dynamic_cast<const SphericalGaussianIE1*>(ie);;
    assert(other);

    mlist.Empty();
    for (auto i:other->es.indices()) mlist.Add(MakeRepulsion((*other)(i)));
    mlist.Clear();
}
//
////
////  This is where we do the big double loop over basis sets.
////
void SphericalGaussianIE1::MakeRepulsion4C(ERIList& Coulomb, ERIList& exchange, const iev_t& iev) const
{
    //No UT coverage
    assert(false);
    // TODO count integrals and % non zero.
    std::vector<const SphericalGaussianIE1*> ie1v;
    size_t N=0;
    for (auto ia: iev)
    {
        ie1v.push_back(dynamic_cast<const SphericalGaussianIE1*>(ia));
        N+=ie1v.back()->es.size();
    }
    Coulomb.SetSize(N);
    exchange.SetSize(N);
    
    std::cout << N << " " << Coulomb.itsData.size() <<" " << exchange.itsData.size() << std::endl;
    ERIList tracker_eris(Coulomb.GetSize()); //For debug purposes.  Check for double assigns
    int La,Lb,Lc,Ld;
    double ea,eb,ec,ed,na,nb,nc,nd;

    int start_a=1;
    for (auto ie_a=ie1v.begin();ie_a!=ie1v.end();ie_a++)
    {
        assert(*ie_a);
        int Na=(*ie_a)->es.size();
        int start_b=start_a;
        for (auto ie_b=ie_a;ie_b!=ie1v.end();ie_b++)
        {
            assert(*ie_b);
            int Nb=(*ie_b)->es.size();
            int start_c=start_a;
            for (auto ie_c=ie_a;ie_c!=ie1v.end();ie_c++)
            {
                assert(*ie_c);
                int Nc=(*ie_c)->es.size();
                int start_d=start_c;
                for (auto ie_d=ie_c;ie_d!=ie1v.end();ie_d++)
                {
                    assert(*ie_d);
                    int Nd=(*ie_d)->es.size();
//                    std::cout << start_a << " " << start_b << " " << start_c << " " << start_d << std::endl;

                    ERIProxy cp(Coulomb     ,start_a,start_b,start_c,start_d);
                    ERIProxy ep(exchange    ,start_a,start_b,start_c,start_d);
                    ERIProxy tp(tracker_eris,start_a,start_b,start_c,start_d);
                    for (index_t ia=1; ia<=Na; ia++)
                    {
                        std::tie(La,ea,na)=(**ie_a)(ia);
                        for (index_t ib=1; ib<=Nb; ib++)
                        {
                            std::tie(Lb,eb,nb)=(**ie_b)(ib);
                            for (index_t ic=1; ic<=Nc; ic++)
                            {
                                std::tie(Lc,ec,nc)=(**ie_c)(ic);
                                for (index_t id=1; id<=Nd; id++)
                                {
                                    std::tie(Ld,ed,nd)=(**ie_d)(id);
                                    double norm=na*nb*nc*nd;
                                    // Coulomb case
                                    double J=0.0,K=0.0;
                                    if (La==Lb && Lc==Ld) //TODO these ifs can go outside the abcd loops.
                                    {
                                        SlaterIntegrals R(ea+eb,ec+ed);
                                        J=FourPi2*R(0,La,Lb,Lc,Ld)*norm;
                                    }
                                    // Exchange case
                                    if (La==Lc && Lb==Ld) //TODO these ifs can go outside the abcd loops.
                                    {
                                        SlaterIntegrals R(ea+eb,ec+ed);
                                        K=FourPi2*R.DoExchangeSum(La,Lb,Lc,Ld)*norm;
                                    }

//                                    std::cout << start_a+ia-1 << " " << start_b+ib-1 << " " << start_c+ic-1 << " " << start_d+id-1 << std::endl;
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

                                } //for id
                            } //for ic
                        } //for ib
                    } //for ia

                    start_d+=Nd;
                } //for ie_d
                start_c+=Nc;
            }
            start_b+=Nb;
        }
        start_a+=Na;
    }
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
    ERIList1 J(N,-1.0),K(N,-1.0);
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



