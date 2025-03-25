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

LAParams DefaultLAP({qchem::Lapack,qchem::SVD,1e-4,1e-12});
//-----------------------------------------------------------------------------
//
//  Construction zone
//
template <class T> TIrrepBasisSetCommon<T>::TIrrepBasisSetCommon()
    : itsLAParams      (DefaultLAP)
{};

template <class T> TIrrepBasisSetCommon<T>::TIrrepBasisSetCommon(const TIrrepBasisSetCommon<T>& bs)
    : itsLAParams      (bs.itsLAParams)
{};

template <class T> TIrrepBasisSetCommon<T>::~TIrrepBasisSetCommon()
{
}

template <class T> void TIrrepBasisSetCommon<T>::Set(const LAParams& lap)
{
    itsLAParams=lap;
} 

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
    for(auto b:IrrepBasisSet::Iterate<TBasisFunction<T> >()) 
    {
        *i=(*b)(r);
        i++;
    }

    return ret;
}

template <class T> typename TIrrepBasisSetCommon<T>::Vec3Vec TIrrepBasisSetCommon<T>::
Gradient(const RVec3& r) const
{
    // No UT coverage
    Vec3Vec  ret(this->size());
    typename Vec3Vec::iterator i(ret.begin());
    for(auto b:IrrepBasisSet::Iterate<TBasisFunction<T> >()) 
    {
        *i=b->Gradient(r);
        i++;
    }

    return ret;
}

template class TIrrepBasisSetCommon<double>;
template class Orbital_IBS_Common<double>;

