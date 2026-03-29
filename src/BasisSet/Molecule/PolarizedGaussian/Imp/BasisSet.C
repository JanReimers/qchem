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
    const IrrepIEClient* a=dynamic_cast<const IrrepIEClient* >(_a);
    const IrrepIEClient* c=dynamic_cast<const IrrepIEClient* >(_c);
    assert(a);
    assert(c);
    size_t Na=a->size(), Nc=c->size();
    ERI4 J(Na,Nc);
    
    for (size_t ia:iv_t(0,Na))
        for (size_t ib:iv_t(ia,Na))
        {
            rsmat_t& Jab=J(ia,ib);
            for (size_t ic:iv_t(0,Nc))
                for (size_t id:iv_t(ic,Nc))
                {
                        //std::cout << "abcd=(" << ia << "," << ib << "," << ic << "," << id << ")" << std::endl;
                        double norm=a->ns[ia]*a->ns[ib]*c->ns[ic]*c->ns[id];
                        assert(c->radials[id]);
                        Jab(ic,id)=norm * c->radials[id]->Integrate(a->radials[ia],a->radials[ib],c->radials[ic],a->pols[ia],a->pols[ib],c->pols[ic],c->pols[id],cache);
                }
        }
    return J;
}

ERI4 BasisSet::MakeExchange(const Orbital_HF_IBS<double>* _a, const Orbital_HF_IBS<double>* _b) const
{
    const IrrepIEClient* a=dynamic_cast<const IrrepIEClient* >(_a);
    const IrrepIEClient* b=dynamic_cast<const IrrepIEClient* >(_b);
    assert(a);
    assert(b);
    size_t Na=a->size(), Nb=b->size();
    ERI4 K(Na,Nb);
    for (size_t ia:iv_t(0,Na))
        for (size_t ib:iv_t(0,Nb))
           
            for (size_t ic:iv_t(ia,Na))
            {
                rsmat_t& Kac=K(ia,ic);
                for (size_t id:iv_t(0,Nb))
                {
                  //std::cout << "abcd=(" << ia << "," << ib << "," << ic << "," << id << ")" << std::endl;
                    double norm=a->ns[ia]*b->ns[ib]*a->ns[ic]*b->ns[id];
                    assert(b->radials[id]);
                    if (ib==id)
                        Kac(ib,id)=norm * b->radials[id]->Integrate(a->radials[ia],b->radials[ib],a->radials[ic],a->pols[ia],b->pols[ib],a->pols[ic],b->pols[id],cache);
                    else if (ib<id)
                        Kac(ib,id)+=0.5*norm * b->radials[id]->Integrate(a->radials[ia],b->radials[ib],a->radials[ic],a->pols[ia],b->pols[ib],a->pols[ic],b->pols[id],cache);
                    else 
                        Kac(id,ib)+=0.5*norm * b->radials[id]->Integrate(a->radials[ia],b->radials[ib],a->radials[ic],a->pols[ia],b->pols[ib],a->pols[ic],b->pols[id],cache);
                }        
            }
    return K;
}

} //namespace
