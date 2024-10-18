// File: OrbitalImplementation.C  General implementation of an orbital, the functional part.



#include "OrbitalImplementation/TOrbitalImplementation.H"
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
template <class T> TOrbitalImplementation<T>::
TOrbitalImplementation(const IDRef<const IrrepBasisSet>& bs,
                       const Vec& c,
                       double e, const Spin& S)
    : OrbitalImplementation(bs,e,S)
    , itsCoeff      (c)
{
    if (!CastBasisSet())
    {
        std::cerr << "TOrbitalImplementation could not cast basis set." << std::endl;
        exit(-1);
    }
};

template <class T> TOrbitalImplementation<T>::TOrbitalImplementation()
    : itsCoeff      ( )
{};

//----------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> T TOrbitalImplementation<T>::operator()(const RVec3& r) const
{
    return itsCoeff * (*CastBasisSet())(r);
}

//BUG
template <class T> typename TOrbitalImplementation<T>::Vec3 TOrbitalImplementation<T>::Gradient(const RVec3& r) const
{
    Vec3 ret(0,0,0);
    Vec3Vec grads=CastBasisSet()->Gradient(r);
    typename Vec3Vec::const_iterator b(grads.begin());
    typename Vec    ::const_iterator c(itsCoeff.begin());
    for (; b!=grads.end(); b++,c++) ret+=(*c) * (*b);
    return ret;
}

template <class T> void TOrbitalImplementation<T>::AddDensityMatrix(SMat& d) const
{
    if (IsOccupied()) d+=SMat(OuterProduct(itsCoeff)*GetOccupation());
}

//-----------------------------------------------------------------------
//
//  Streamabel stuff.
//
template <class T> std::ostream& TOrbitalImplementation<T>::Write(std::ostream& os) const
{
    OrbitalImplementation::Write(os);
    os << itsCoeff;
    if (Pretty()) os << std::endl;
    return os;
}

template <class T> std::istream& TOrbitalImplementation<T>::Read (std::istream& is)
{
    OrbitalImplementation::Read(is);
    is >> itsCoeff;

    return is;
}


template class TOrbitalImplementation<double>;
template class TOrbitalImplementation<std::complex<double> >;
