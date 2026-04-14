// File PolarizedGaussian/Imp/BasisSet.C
module;
#include <memory>
#include <cmath>
#include <cassert>
#include <vector>

// namespace PolarizedGaussian{class Reader;} /* g++-15.2 BUG? not handling forward class decs as well as clang++ 20,21*/

module qchem.BasisSet.Molecule.PolarizedGaussian;
import qchem.BasisSet.Molecule.PolarizedGaussian.Reader;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.IEClient;
import qchem.Cluster;

namespace PolarizedGaussian
{


BasisSet::BasisSet( Reader* reader, const Cluster* cl)
{
    Insert(new Orbital_IBS(this,reader,cl));
}

void BasisSet::Insert(bs_t* bs)
{
    BS_Common::Insert(bs);
    auto oibs=dynamic_cast<const Orbital_HF_IBS<double>*>(bs);
    assert(oibs);
    Append(oibs);
}

ERI4 BasisSet::MakeDirect  (const Orbital_HF_IBS<double>* _a, const Orbital_HF_IBS<double>* _c) const
{
    const PGData* a=dynamic_cast<const PGData* >(_a);
    const PGData* c=dynamic_cast<const PGData* >(_c);
    assert(a);
    assert(c);
    size_t Na=a->size1(), Nc=c->size1();
    ERI4 J(Na,Nc);
    
    for (size_t ia:iv_t(0,Na))
        for (size_t ib:iv_t(ia,Na))
        {
            rsmat_t& Jab=J(ia,ib);
            for (size_t ic:iv_t(0,Nc))
                for (size_t id:iv_t(ic,Nc))
                {
                        //std::cout << "abcd=(" << ia << "," << ib << "," << ic << "," << id << ")" << std::endl;
                        double norm=a->ns1[ia]*a->ns1[ib]*c->ns1[ic]*c->ns1[id];
                        assert(c->radials1[id]);
                        Jab(ic,id)=norm * c->radials1[id]->Integrate(a->radials1[ia],a->radials1[ib],c->radials1[ic],a->pols1[ia],a->pols1[ib],c->pols1[ic],c->pols1[id],cache);
                }
        }
    return J;
}

ERI4 BasisSet::MakeExchange(const Orbital_HF_IBS<double>* _a, const Orbital_HF_IBS<double>* _b) const
{
    const PGData* a=dynamic_cast<const PGData* >(_a);
    const PGData* b=dynamic_cast<const PGData* >(_b);
    assert(a);
    assert(b);
    size_t Na=a->size1(), Nb=b->size1();
    ERI4 K(Na,Nb);
    for (size_t ia:iv_t(0,Na))
        for (size_t ib:iv_t(0,Nb))
           
            for (size_t ic:iv_t(ia,Na))
            {
                rsmat_t& Kac=K(ia,ic);
                for (size_t id:iv_t(0,Nb))
                {
                  //std::cout << "abcd=(" << ia << "," << ib << "," << ic << "," << id << ")" << std::endl;
                    double norm=a->ns1[ia]*b->ns1[ib]*a->ns1[ic]*b->ns1[id];
                    assert(b->radials1[id]);
                    if (ib==id)
                        Kac(ib,id)=norm * b->radials1[id]->Integrate(a->radials1[ia],b->radials1[ib],a->radials1[ic],a->pols1[ia],b->pols1[ib],a->pols1[ic],b->pols1[id],cache);
                    else if (ib<id)
                        Kac(ib,id)+=0.5*norm * b->radials1[id]->Integrate(a->radials1[ia],b->radials1[ib],a->radials1[ic],a->pols1[ia],b->pols1[ib],a->pols1[ic],b->pols1[id],cache);
                    else 
                        Kac(id,ib)+=0.5*norm * b->radials1[id]->Integrate(a->radials1[ia],b->radials1[ib],a->radials1[ic],a->pols1[ia],b->pols1[ib],a->pols1[ic],b->pols1[id],cache);
                }        
            }
    return K;
}

} //namespace
