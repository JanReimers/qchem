// File: AtomIE.C Common IE code for all atom basis sets.

#include "Imp/BasisSet/AtomIE.H"
#include "Imp/BasisSet/AtomIEClient.H"
#include <Cluster.H>
#include "oml/smatrix.h"

const AtomIEClient* AtomIE::dcast(iec_t* iea)
{
    const AtomIEClient* a=dynamic_cast<const AtomIEClient*>(iea);
    assert(a);
    return a;
}


AtomIE::SMat AtomIE::MakeOverlap(iec_t* iea ) const
{
    auto  a=dcast(iea);
    size_t N=a->size();
    SMat s(N);
    for (auto i:s.rows())
        for (auto j:s.cols(i))
            s(i,j)=Overlap(a->es(i),a->es(j),a->Ls(i))*a->ns(i)*a->ns(j);

    return s;
}

AtomIE::SMat AtomIE::MakeKinetic(iec_t* iea) const
{
    auto a=dcast(iea);;
    size_t N=a->size();
    SMatrix<double> Hk(N);
    for (auto i:Hk.rows())
        for (auto j:Hk.cols(i))
        {
            Hk(i,j)=Kinetic(a->es(i),a->es(j),a->Ls(i))*a->ns(i)*a->ns(j);
            assert(a->Ls(i)==a->Ls(j));
        }

    return Hk;
}

AtomIE::SMat AtomIE::MakeNuclear(iec_t* iea,const Cluster& cl) const
{
    auto a=dcast(iea);;
    size_t N=a->size(),L=a->Ls(1);
    SMatrix<double> Hn(N);
    double Z=-cl.GetNuclearCharge();
    for (auto i:Hn.rows())
        for (auto j:Hn.cols(i))
            Hn(i,j)= Z*Nuclear(a->es(i),a->es(j),L)*a->ns(i)*a->ns(j);

    return Hn;
}

AtomIE::RVec AtomIE::MakeCharge(iec_t* iea) const
{
    auto a=dcast(iea);;
    RVec c(a->size());
    for (auto i:a->es.indices())  c(i)=Charge(a->es(i),a->Ls(i))*a->ns(i);
    return c;
}

