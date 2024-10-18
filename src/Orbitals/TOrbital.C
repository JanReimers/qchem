// File: OrbitalImplementation.C  General implementation of an orbital, the functional part.



#include "Imp/Orbitals/TOrbital.H"
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
TOrbitalImp(const IDRef<const IrrepBasisSet>& bs,
                       const Vec& c,
                       double e, const Spin& S)
    : OrbitalImp(bs,e,S)
    , itsCoeff      (c)
{
    if (!CastBasisSet())
    {
        std::cerr << "TOrbitalImplementation could not cast basis set." << std::endl;
        exit(-1);
    }
};

template <class T> TOrbitalImp<T>::TOrbitalImp()
    : itsCoeff      ( )
{};

//----------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> T TOrbitalImp<T>::operator()(const RVec3& r) const
{
    return itsCoeff * (*CastBasisSet())(r);
}

//BUG
template <class T> typename TOrbitalImp<T>::Vec3 TOrbitalImp<T>::Gradient(const RVec3& r) const
{
    Vec3 ret(0,0,0);
    Vec3Vec grads=CastBasisSet()->Gradient(r);
    typename Vec3Vec::const_iterator b(grads.begin());
    typename Vec    ::const_iterator c(itsCoeff.begin());
    for (; b!=grads.end(); b++,c++) ret+=(*c) * (*b);
    return ret;
}

template <class T> void TOrbitalImp<T>::AddDensityMatrix(SMat& d) const
{
    if (IsOccupied()) d+=SMat(OuterProduct(itsCoeff)*GetOccupation());
}

//-----------------------------------------------------------------------
//
//  Streamabel stuff.
//
template <class T> std::ostream& TOrbitalImp<T>::Write(std::ostream& os) const
{
    OrbitalImp::Write(os);
    os << itsCoeff;
    if (Pretty()) os << std::endl;
    return os;
}

template <class T> std::istream& TOrbitalImp<T>::Read (std::istream& is)
{
    OrbitalImp::Read(is);
    is >> itsCoeff;

    return is;
}


template class TOrbitalImp<double>;
template class TOrbitalImp<std::complex<double> >;
