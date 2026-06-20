// File: BasisSet/Molecule/PolarizedGaussian1/Imp/Orbital_IBS.C  Polarized Gaussian fit basis set, for MO calculations.
module;
#include <cassert>
#include <algorithm> //Need std::max
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <vector>

module qchem.BasisSet.Molecule.PolarizedGaussian1;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.GaussianRF;
import qchem.BasisSet.Molecule.PolarizedGaussian1.Internal.Readers.Gaussian94;
import qchem.BasisSet.Internal.IntegralEnums;
import qchem.BasisSet;
import qchem.Cluster;
import qchem.Symmetry.Unit;
import qchem.stl_io;
import qchem.Streamable;
import qchem.Math;
import qchem.Blaze;


namespace BasisSet::Molecule::PolarizedGaussian1
{

rsmat_t Orbital_IBS::Integrate(qchem::IType3C type , const GaussianRF* rc, const Polarization& pc) const
{
    auto ab=dynamic_cast<const PGData*>(this);
    int N=ab->size();
    rsmat_t s(N);
    for (size_t ia=0;ia<N;ia++)
        for (size_t ib=ia;ib<N;ib++)
            s(ia,ib)=rc->Integrate(type,*ab->radials[ia],*ab->radials[ib],ab->pols[ia],ab->pols[ib],pc)*ab->ns[ia]*ab->ns[ib];
        
    return s;    
}

rsmat_t MakeIntegrals(IType t2C,const PGData* ab, const Cluster* cl) 
{
    assert(ab);
    int N=ab->size();
    rsmat_t s(N);
    for (size_t ia=0;ia<N;ia++)
        for (size_t ib=ia;ib<N;ib++)
            s(ia,ib)=ab->radials[ia]->Integrate(t2C,*ab->radials[ib],ab->pols[ia],ab->pols[ib],cl)*ab->ns[ia]*ab->ns[ib];

    return s;
}

ERI3<double> Orbital_IBS::MakeOverlap3C(const Fit_IBS& _c) const
{
    auto c=dynamic_cast<const PGData*>(&_c);
    int Nc=c->size();
    ERI3<double> s3;
    for (size_t ic=0;ic<Nc;ic++)
    {
        rsmat_t s=Integrate(qchem::Overlap3C,c->radials[ic],c->pols[ic]);
        s*=c->ns[ic];
        s3.push_back(s);
    } 
    return s3;   
}
ERI3<double> Orbital_IBS::MakeRepulsion3C(const Fit_IBS& _c) const
{
    auto c=dynamic_cast<const PGData*>(&_c);
    int Nc=c->size();
    ERI3<double> s3;
    for (size_t ic=0;ic<Nc;ic++)
    {
        rsmat_t s=Integrate(qchem::Repulsion3C,c->radials[ic],c->pols[ic]);
        s*=c->ns[ic];
        s3.push_back(s);
    }    
    return s3;
}

ERI4 Orbital_IBS::MakeDirect  (const Orbital_HF_IBS<double>& _c) const
{
    const PGData* a=this;
    const PGData* c=dynamic_cast<const PGData* >(&_c);
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
                        Jab(ic,id)=norm * c->radials[id]->Integrate(*a->radials[ia],*a->radials[ib],*c->radials[ic],a->pols[ia],a->pols[ib],c->pols[ic],c->pols[id]);
                }
        }
    return J;
}

ERI4 Orbital_IBS::MakeExchange(const Orbital_HF_IBS<double>& _b) const
{
    const PGData* a=this;
    const PGData* b=dynamic_cast<const PGData* >(&_b);
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
                        Kac(ib,id)=norm * b->radials[id]->Integrate(*a->radials[ia],*b->radials[ib],*a->radials[ic],a->pols[ia],b->pols[ib],a->pols[ic],b->pols[id]);
                    else if (ib<id)
                        Kac(ib,id)+=0.5*norm * b->radials[id]->Integrate(*a->radials[ia],*b->radials[ib],*a->radials[ic],a->pols[ia],b->pols[ib],a->pols[ic],b->pols[id]);
                    else 
                        Kac(id,ib)+=0.5*norm * b->radials[id]->Integrate(*a->radials[ia],*b->radials[ib],*a->radials[ic],a->pols[ia],b->pols[ib],a->pols[ic],b->pols[id]);
                }        
            }
    return K;
}


} //namespace