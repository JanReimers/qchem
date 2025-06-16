// File PolarizedGaussian/BasisSet.C

#include "Molecule/PolarizedGaussian/BasisSet.H"
#include "Molecule/PolarizedGaussian/IrrepBasisSet.H"
#include <Cluster.H>

namespace PolarizedGaussian
{


BasisSet::BasisSet( Reader* reader, const Cluster* cl)
{
    Insert(new Orbital_IBS(this,reader,cl));
}

void BasisSet::Insert(bs_t* bs)
{
    BS_Common::Insert(bs);
    auto iec=dynamic_cast<const IrrepIEClient*>(bs);
    assert(iec);
    Append(iec);
}

ERI4 BasisSet::MakeDirect  (const ::IrrepIEClient* _a, const ::IrrepIEClient* _c) const
{
    const IrrepIEClient* a=dynamic_cast<const IrrepIEClient* >(_a);
    const IrrepIEClient* c=dynamic_cast<const IrrepIEClient* >(_c);
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    ERI4 J(Na,Nc);
    
    for (index_t ia:a->ns.indices())
        for (index_t ib:a->ns.indices(ia))
        {
            SMatrix<double>& Jab=J(ia,ib);
            for (index_t ic:c->ns.indices())
                for (index_t id:c->ns.indices(ic))
                {
                        //std::cout << "abcd=(" << ia << "," << ib << "," << ic << "," << id << ")" << std::endl;
                        double norm=a->ns(ia)*a->ns(ib)*c->ns(ic)*c->ns(id);
                        assert(c->radials[id-1]);
                        Jab(ic,id)=norm * c->radials[id-1]->Integrate(a->radials[ia-1],a->radials[ib-1],c->radials[ic-1],a->pols[ia-1],a->pols[ib-1],c->pols[ic-1],c->pols[id-1],cache);
                }
        }
    return J;
}

ERI4 BasisSet::MakeExchange(const ::IrrepIEClient* _a, const ::IrrepIEClient* _b) const
{
    const IrrepIEClient* a=dynamic_cast<const IrrepIEClient* >(_a);
    const IrrepIEClient* b=dynamic_cast<const IrrepIEClient* >(_b);
    assert(a);
    assert(b);
    size_t Na=a->size(), Nb=b->size();
    ERI4 K(Na,Nb);
    for (index_t ia:a->ns.indices())
        for (index_t ib:b->ns.indices())
           
            for (index_t ic:a->ns.indices(ia))
            {
                SMatrix<double>& Kac=K(ia,ic);
                for (index_t id:b->ns.indices())
                {
                  //std::cout << "abcd=(" << ia << "," << ib << "," << ic << "," << id << ")" << std::endl;
                    double norm=a->ns(ia)*b->ns(ib)*a->ns(ic)*b->ns(id);
                    assert(b->radials[id-1]);
                    if (ib==id)
                        Kac(ib,id)=norm * b->radials[id-1]->Integrate(a->radials[ia-1],b->radials[ib-1],a->radials[ic-1],a->pols[ia-1],b->pols[ib-1],a->pols[ic-1],b->pols[id-1],cache);
                    else if (ib<id)
                        Kac(ib,id)+=0.5*norm * b->radials[id-1]->Integrate(a->radials[ia-1],b->radials[ib-1],a->radials[ic-1],a->pols[ia-1],b->pols[ib-1],a->pols[ic-1],b->pols[id-1],cache);
                    else 
                        Kac(id,ib)+=0.5*norm * b->radials[id-1]->Integrate(a->radials[ia-1],b->radials[ib-1],a->radials[ic-1],a->pols[ia-1],b->pols[ib-1],a->pols[ic-1],b->pols[id-1],cache);
                }        
            }
    return K;
}

} //namespace
