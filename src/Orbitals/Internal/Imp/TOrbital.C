// File: OrbitalImp.C  Implementation of an orbital.
module;
#include <iostream>
#include <iomanip>
#include <complex>
#include "blaze/Math.h" 
module qchem.Orbitals.Internal.OrbitalImp;
import qchem.IrrepBasisSet;
import qchem.Conversions;

//-----------------------------------------------------------------
//
//  Construction zone
//
template <class T> TOrbitalImp<T>::
TOrbitalImp(const Orbital_IBS<T>* bs,const vec_t<T>& _C,const vec_t<T>& _CPrime,double e, const Orbital_QNs& qns)
    : OrbitalImp   (e,qns)
    , itsCoeff     (_C)
    , itsCoeffPrime(_CPrime)
    , itsBasisSet  (bs)
{
};

template <class T> void TOrbitalImp<T>::AddDensityMatrix(smat_t<T>& D, smat_t<T>& DPrime) const
{
    
    if (IsOccupied()) 
    {
        smat_t<T> CCd=blaze::outer(itsCoeff,conj(itsCoeff))*GetOccupation();
        D+=CCd;
        smat_t<T> CCd_prime=blaze::outer(itsCoeffPrime,conj(itsCoeffPrime))*GetOccupation();
        DPrime+=CCd_prime;
        // D     +=SMatrix<T>(OuterProduct(itsCoeff     )*GetOccupation());
        // DPrime+=SMatrix<T>(OuterProduct(itsCoeffPrime)*GetOccupation());
    }
}



//----------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> T TOrbitalImp<T>::operator()(const RVec3& r) const
{
    return convert(itsCoeff) * (*itsBasisSet)(r);
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
    os << "              " << GetOccupation() << "/" << GetDegeneracy() << "       " << std::setw(12) << GetEigenEnergy() << "      ";
    os << itsCoeff << std::endl;
    return os;
}


template class TOrbitalImp<double>;
// template class TOrbitalImp<std::complex<double> >;
