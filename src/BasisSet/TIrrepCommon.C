// File: TBasisSetImplementation.C


#include "Imp/BasisSet/TIrrepCommon.H"
#include "Imp/Containers/ERI4.H"
#include <BasisFunction.H>
#include <QuantumNumber.H>
#include <AnalyticIE.H>
#include <LASolver.H>
#include <Hamiltonian.H>
#include "oml/vector.h"
#include <cassert>
#include <iostream>

//-----------------------------------------------------------------------------
//
//  Construction zone
//
template <class T> TIrrepBasisSetCommon<T>::TIrrepBasisSetCommon()
{};

template <class T> TIrrepBasisSetCommon<T>::TIrrepBasisSetCommon(const LAParams& lap)
    : itsLAParams      (lap)
{};

template <class T> TIrrepBasisSetCommon<T>::TIrrepBasisSetCommon(const TIrrepBasisSetCommon<T>& bs)
    : itsLAParams      (bs.itsLAParams)
{};

template <class T> TIrrepBasisSetCommon<T>::~TIrrepBasisSetCommon()
{
}

// template <class T> IntegralDataBase<T>* TIrrepBasisSetCommon<T>::GetDataBase() const
// {
//     assert(&*itsDataBase);
//     return itsDataBase;
// }

 

template <class T>  LASolver<double>* Orbital_IBS_Common<T>::CreateSolver() const
{
    StreamableObject::SetToPretty();
    // std::cout << "S_old=" << this->GetOverlap() << std::endl;
    // std::cout << "S_new=" << this->Integrals(qchem::Overlap1,this) << std::endl;
    LASolver<double>* las=LASolver<double>::Factory(TIrrepBasisSetCommon<T>::itsLAParams);
    //las->SetBasisOverlap(this->GetOverlap());
    las->SetBasisOverlap(this->Overlap());
    return las;
}


 using std::cout;
 using std::endl;
//-----------------------------------------------------------------------------
//
//  VectorFunction stuff.
//
template <class T> typename TIrrepBasisSetCommon<T>::Vec TIrrepBasisSetCommon<T>::
operator() (const RVec3& r) const
{
    Vec  ret(this->size());
    typename Vec::iterator i(ret.begin());
    for(auto b=this->beginT();b!=this->end();i++,b++) *i=(**b)(r);

    return ret;
}

template <class T> typename TIrrepBasisSetCommon<T>::Vec3Vec TIrrepBasisSetCommon<T>::
Gradient(const RVec3& r) const
{
    // No UT coverage
    Vec3Vec  ret(this->size());
    typename Vec3Vec::iterator i(ret.begin());
    for(auto b=this->beginT(); b!=this->end(); i++,b++) *i=b->Gradient(r);

    return ret;
}

template <class T> typename Orbital_DFT_IBS_Common<T>::Vec Orbital_DFT_IBS_Common<T>::
Overlap3C(const SMat& Dcd, const fbs_t* ff) const
{
    Vec ret(ff->size());
    const typename Integrals_DFT<T>::ERI3& S=this->Overlap3C(*ff);
    for(auto i:ret.indices())
        ret(i)=Dot(Dcd,S[i-1]);
    return ret;
}

template <class T> typename Orbital_DFT_IBS_Common<T>::Vec Orbital_DFT_IBS_Common<T>::
Repulsion3C(const SMat& Dcd, const fbs_t* ff) const
{
    Vec ret(ff->size());
    const typename Integrals_DFT<T>::ERI3& repulsion=this->Repulsion3C(*ff);
    for(auto i:ret.indices())
        ret(i)=Dot(Dcd,repulsion[i-1]);
    return ret;
}

template <class T> typename Orbital_HF_IBS_Common<T>::SMat Orbital_HF_IBS_Common<T>::
Direct(const SMat& Dcd, const obs_t* cd) const
{
    assert(!isnan(Dcd));
    assert(Max(fabs(Dcd))>0.0);  //Don't waste time!
    const TOrbital_HF_IBS<T>* ab=this;
    
    if (ab->GetID()<=cd->GetID())
    {
         ERI4 Jabcd=ab->Direct(*cd);
        return Jabcd*Dcd;
    }
    else
    {
        //ERI4 Jcdab=GetDataBase()->GetDirect__4C(*cd,*ab);
        const TOrbital_HF_IBS<T>* cdhf=dynamic_cast<const TOrbital_HF_IBS<T>*>(cd);
        assert(cdhf);
        ERI4 Jcdab=cdhf->Direct(*ab);
        return Dcd*Jcdab;        
    }
}

#include <iomanip>
template <class T> typename Orbital_HF_IBS_Common<T>::SMat Orbital_HF_IBS_Common<T>::
Exchange(const SMat& Dcd, const obs_t* cd) const
{
    assert(!isnan(Dcd));
    assert(Max(fabs(Dcd))>0.0);  //Don't waste time!
    const TOrbital_HF_IBS<T>* ab=this;

    if (ab->GetID()<=cd->GetID())
    {
        ERI4 Kabcd=ab->Exchange(*cd);
        return Kabcd*Dcd;
    }
    else
    {
        //ERI4 Kcdab=GetDataBase()->GetExchange4C(*cd,*ab);
        const TOrbital_HF_IBS<T>* cdhf=dynamic_cast<const TOrbital_HF_IBS<T>*>(cd);
        assert(cdhf);
        ERI4 Kcdab=cdhf->Exchange(*ab);
        return Dcd*Kcdab;        
    }

}

#include "Imp/Integrals/MeshIntegrator.H"

Fit_IBS_Common::Vec Fit_IBS_Common::MakeNorm   (const Mesh* m) const
{
    MeshIntegrator<double> mintegrator(m);
    return mintegrator.Normalize(*this);
}
Fit_IBS_Common::Vec Fit_IBS_Common::MakeCharge (const Mesh*  m) const
{
    assert(false);
    return *new Vec();
}
Fit_IBS_Common::Mat Fit_IBS_Common::MakeOverlap(const Mesh* m,const fbs_t& b) const
{
    assert(false);
    return *new Mat();
}
const Fit_IBS_Common::Vec Fit_IBS_Common::Overlap  (const Mesh* m,const Sf& f) const
{
    const Vec& n=Norm(m);
    MeshIntegrator<double> mintegrator(m);
    return DirectMultiply(mintegrator.Overlap(f,*this),n);
}  
const Fit_IBS_Common::Vec Fit_IBS_Common::Repulsion(const Mesh* m,const Sf& f) const
{
    const Vec& n=Norm(m);
    MeshIntegrator<double> mintegrator(m);
    return DirectMultiply(mintegrator.Repulsion(f,*this),n);
}

template class TIrrepBasisSetCommon<double>;
template class Orbital_IBS_Common<double>;
template class Orbital_DFT_IBS_Common<double>;
template class Orbital_HF_IBS_Common<double>;
//template class TBasisSetImplementation<std::complex<double> >;

