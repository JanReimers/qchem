// File: OrbitalImplementation.C  General implementation of an orbital, the functional part.



#include "Imp/Orbitals/TOrbital.H"
#include <Irrep_BS.H>
#include "oml/smatrix.h"
#include <iostream>

//-----------------------------------------------------------------
//
//  Construction zone
//
template <class T> TOrbitalImp<T>::
TOrbitalImp(const TOrbital_IBS<T>* bs,const Vec& c,double e, const Orbital_QNs& qns)
    : OrbitalImp (e,qns)
    , itsCoeff   (c)
    , itsBasisSet(bs)
{
};

template <class T> void TOrbitalImp<T>::AddDensityMatrix(SMat& d) const
{
    
    if (IsOccupied()) 
        d+=SMat(OuterProduct(itsCoeff)*GetOccupation());
}



//----------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> T TOrbitalImp<T>::operator()(const RVec3& r) const
{
    return itsCoeff * (*itsBasisSet)(r);
}

//BUG
template <class T> typename TOrbitalImp<T>::Vec3 TOrbitalImp<T>::Gradient(const RVec3& r) const
{
    Vec3 ret(0,0,0);
    Vec3Vec grads=itsBasisSet->Gradient(r);
    auto c(itsCoeff.begin());
    for (auto b:grads) ret+=(*c++) * b;
    return ret;
}


//-----------------------------------------------------------------------
//
//  Streamabel stuff.
//
template <class T> std::ostream& TOrbitalImp<T>::Write(std::ostream& os) const
{
    OrbitalImp::Write(os);
    if (StreamableObject::Pretty())
        os << "              " << GetOccupation() << "/" << GetDegeneracy() << "       " << std::setw(12) << GetEigenEnergy() << "      ";
    else
        os << itsBasisSet;
    
    
        
    os << itsCoeff;
    if (StreamableObject::Pretty()) os << std::endl;
    return os;
}


template class TOrbitalImp<double>;
template class TOrbitalImp<std::complex<double> >;
