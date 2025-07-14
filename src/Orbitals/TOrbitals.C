// File: OrbitalGroupImplementation.C  general orbital group implementation.



#include "TOrbitals.H"
#include "TOrbital.H"
#include <ChargeDensity/Factory.H>
#include <BasisSet/Irrep_BS.H>
#include <Symmetry/Symmetry.H>
#include "Common/stl_io.h"

import Common.Constants; //c_light

//-----------------------------------------------------------------
//
//  Construction zone
//
template <class T> TOrbitalsImp<T>::
TOrbitalsImp(const TOrbital_IBS<T>* bs, Spin ms)
    : itsBasisSet(bs)
    , itsQNs(ms,bs->GetSymmetry())
    , itsD     (bs->GetNumFunctions())
{
    assert(itsBasisSet->GetNumFunctions()>0);
    
    Fill(itsD,0.0);
  
};

template <class T> TOrbitalsImp<T>::~TOrbitalsImp()
{
 
}

//-----------------------------------------------------------------
//
//  Orbitals stuff.
//
template <class T> index_t TOrbitalsImp<T>::GetNumOrbitals() const
{
    return itsOrbitals.size();
}
template <class T> index_t TOrbitalsImp<T>::GetNumOccOrbitals() const
{
    index_t n=0;
    for (auto& o:*this)
        if (o->IsOccupied()) 
            n++;

    return n;
}
template <class T> double TOrbitalsImp<T>::GetEigenValueChange(const Orbitals& og) const
{
    // No UT coverage
    // TODO: OrbitalGroup should return a vector of energies.
    double del=0;
    auto b2=og.Iterate<Orbital>().begin();
    for (auto b1:Iterate<Orbital>())
    {
        const Orbital* o2=*b2;
        del+=Square(b1->GetEigenEnergy()-o2->GetEigenEnergy());
        ++b2;
    }
    return sqrt(del);
}

//
//  This is where the real SCF work gets done.
//
template <class T> void TOrbitalsImp<T>::UpdateOrbitals(const Mat& U, const Mat& UPrime, const RVec& e)
{
    itsOrbitals.clear();
    index_t n=e.size();
    //
    //  Strip out all the positron orbitals.
    //
    static const double e_positron=-c_light*c_light; //Max positron energy = -2mc^2.
    size_t index=1;
    for (index_t i=1; i<=n; i++)
    {
 //               std::cout << "o=" << o->GetEigenEnergy() << std::endl;
        if (e(i)<=e_positron) continue; //Strip out all the positron orbitals.
        size_t principle_QN=itsQNs.sym->GetPrincipleOffset() + index++;
        Orbital* o=new TOrbitalImp<T>(itsBasisSet,U.GetColumn(i), UPrime.GetColumn(i), e(i),Orbital_QNs(principle_QN,itsQNs));
        itsOrbitals.push_back(std::unique_ptr<Orbital>(o));

    }
}
template <class T> typename TOrbitalsImp<T>::ds_t TOrbitalsImp<T>::TakeElectrons(double ne)
{
    // Dump electrons into orbitals, starting from the lowest energy.
    for (auto o:this->Iterate<Orbital>())
    {
        ne=o->TakeElectrons(ne);
        if (ne<=0.0) break;
    }
    //
    //  Now the orbitals are accupied we can build the density matrix.
    //
    Fill(itsD,T(0.0));
    SMat DPrime(itsD.GetLimits());
    Fill(DPrime,T(0.0));
    for (auto o:Iterate<TOrbital<double>>()) o->AddDensityMatrix(itsD,DPrime);
    return std::make_tuple(ne,DPrime);
}

template <class T> DM_CD* TOrbitalsImp<T>::GetChargeDensity() const
{
    return IrrepCD_Factory(itsD,itsBasisSet,GetQNs());
}


template <class T>  Irrep_QNs TOrbitalsImp<T>::GetQNs() const
{
    return itsQNs;
}


//-----------------------------------------------------------------
//
//  VectorFunction stuff.
//
template <class T> typename TOrbitalsImp<T>::Vec TOrbitalsImp<T>::
operator()(const RVec3& r) const
{
    Vec ret(GetNumOrbitals());
    typename Vec::iterator i(ret.begin());
    // No UT coverage
    for (auto b:Iterate<TOrbital<double>>()) 
    {
        *i=(*b)(r);
        i++;
    }
    return ret;
}

template <class T> typename TOrbitalsImp<T>::Vec3Vec TOrbitalsImp<T>::
Gradient(const RVec3& r) const
{
    // No UT coverage
    Vec3Vec ret(GetNumOrbitals());
    typename Vec3Vec::iterator i(ret.begin());
    for (auto b:Iterate<TOrbital<double>>()) 
    {
        *i=b->Gradient(r);
        i++;
    } 
    return ret;
}

//-----------------------------------------------------------------
//
//  Streamable stuff.
//
template <class T> std::ostream& TOrbitalsImp<T>::Write(std::ostream& os) const
{
     if (!StreamableObject::Pretty())
    {
        os  << itsOrbitals;
        if (!StreamableObject::Binary()) os << std::endl;
    }
    else
    {
        os << "        Orbital group with " << GetNumOrbitals() << " " << itsBasisSet->GetSymmetry() << "orbitals:" << std::endl;
        os << "            Occupation      Energy      Eigenvector" << std::endl;
        os << itsOrbitals;
    }
    if (!StreamableObject::Pretty())
    {
        os  << itsBasisSet;
        if (StreamableObject::Ascii()) os << std::endl;
    }
    return os;
}


template class TOrbitalsImp<double>;
