// File: TBasisSetImplementation.C


#include "Imp/BasisSet/TIBS_Common.H"
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
template <class T> TIBS_Common<T>::TIBS_Common()
    : itsLAParams      (DefaultLAP)
{};

template <class T> TIBS_Common<T>::TIBS_Common(const TIBS_Common<T>& bs)
    : itsLAParams      (bs.itsLAParams)
{};

template <class T> TIBS_Common<T>::~TIBS_Common()
{
}

template <class T> void TIBS_Common<T>::Set(const LAParams& lap)
{
    itsLAParams=lap;
} 

template <class T>  LASolver<double>* Orbital_IBS_Common<T>::CreateSolver() const
{
    StreamableObject::SetToPretty();
    // std::cout << "S_old=" << this->GetOverlap() << std::endl;
    // std::cout << "S_new=" << this->Integrals(qchem::Overlap1,this) << std::endl;
    LASolver<double>* las=LASolver<double>::Factory(TIBS_Common<T>::itsLAParams);
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
template <class T> typename TIBS_Common<T>::Vec TIBS_Common<T>::
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

template <class T> typename TIBS_Common<T>::Vec3Vec TIBS_Common<T>::
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

template class TIBS_Common<double>;
template class Orbital_IBS_Common<double>;

