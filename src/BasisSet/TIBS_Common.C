// File: TBasisSetImplementation.C


#include <tuple>
#include <iostream>
#include <vector>
#include <LASolver/LAParams.H>
#include <LASolver/LASolver.H>
#include "TIBS_Common.H"
import qchem.BasisFunction;
import oml;

LAParams DefaultLAP({qchem::Lapack,qchem::SVD,1e-10,1e-12});
//-----------------------------------------------------------------------------
//
//  Construction zone
//
template <class T> TIBS_Common<T>::TIBS_Common()
    : itsLAParams      (DefaultLAP) //gcc-15.0.1 segfault here
{
};

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
    LASolver<double>* las=LASolver<double>::Factory(TIBS_Common<T>::itsLAParams);
    las->SetBasisOverlap(this->Overlap());
    return las;
}


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

