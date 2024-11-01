// File: OrbitalImplementation.C  General implementation of an orbital, the functional part.



#include "Imp/Orbitals/TOrbital.H"
#include <QuantumNumber.H>
#include "oml/vector.h"
#include "oml/vector3d.h"
#include "oml/smatrix.h"
#include "oml/matrix.h"
#include <iostream>
#include <stdlib.h>

//-----------------------------------------------------------------
//
//  Construction zone
//
template <class T> TOrbitalImp<T>::
TOrbitalImp(const TIrrepBasisSet<T>* bs,const Vec& c,double e, const Spin& S)
    : ElectronContainerImp(S,bs->GetQuantumNumber())
    , itsEigenEnergy(e)
    , itsCoeff      (c)
    , itsBasisSet   (bs)
{
};

template <class T> TOrbitalImp<T>::TOrbitalImp()
    : itsCoeff      ( )
{};

template <class T> double TOrbitalImp<T>::GetEigenEnergy() const
{
    return itsEigenEnergy;
}

template <class T> void TOrbitalImp<T>::AddDensityMatrix(SMat& d) const
{
    
    if (IsOccupied()) 
    {
//        std::cout << "Orbital occ=" << GetOccupation() << std::fixed << " E=" << GetEigenEnergy() 
//        << " QN=" << itsBasisSet->GetQuantumNumber() << std::endl;
        d+=SMat(OuterProduct(itsCoeff)*GetOccupation());
    }
}

template <class T> const QuantumNumber& TOrbitalImp<T>::GetQuantumNumber() const
{
    return itsBasisSet->GetQuantumNumber();
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
    ElectronContainerImp::Write(os);
    if (StreamableObject::Pretty())
        os << "              " << GetOccupation() << "/" << GetDegeneracy() << "       " << std::setw(12) << itsEigenEnergy << "      ";
    else
        os << itsBasisSet;
    
    if (StreamableObject::Binary())
        BinaryWrite(itsEigenEnergy,os);
    if (StreamableObject::Ascii ())
        os << itsEigenEnergy << " ";
        
    os << itsCoeff;
    if (StreamableObject::Pretty()) os << std::endl;
    return os;
}

template <class T> std::istream& TOrbitalImp<T>::Read (std::istream& is)
{
    ElectronContainerImp::Read(is);
    if (StreamableObject::Binary())
        BinaryRead(itsEigenEnergy,is);
    else
        is >> itsEigenEnergy >> std::ws;
        
    is >> itsCoeff;

    return is;
}


template class TOrbitalImp<double>;
template class TOrbitalImp<std::complex<double> >;
