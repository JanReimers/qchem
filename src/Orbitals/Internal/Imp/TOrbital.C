// File: OrbitalImp.C  Implementation of an orbital.
module;
#include <iostream>
#include <iomanip>
#include <complex>
module qchem.Orbitals.Internal.OrbitalImp;
import qchem.Irrep_BS;

//-----------------------------------------------------------------
//
//  Construction zone
//
template <class T> TOrbitalImp<T>::
TOrbitalImp(const TOrbital_IBS<T>* bs,const Vec& _C,const Vec& _CPrime,double e, const Orbital_QNs& qns)
    : OrbitalImp   (e,qns)
    , itsCoeff     (_C)
    , itsCoeffPrime(_CPrime)
    , itsBasisSet  (bs)
{
};

template <class T> void TOrbitalImp<T>::AddDensityMatrix(SMatrix<T>& D, SMatrix<T>& DPrime) const
{
    
    if (IsOccupied()) 
    {
        D     +=SMatrix<T>(OuterProduct(itsCoeff     )*GetOccupation());
        DPrime+=SMatrix<T>(OuterProduct(itsCoeffPrime)*GetOccupation());
    }
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
    os << "              " << GetOccupation() << "/" << GetDegeneracy() << "       " << std::setw(12) << GetEigenEnergy() << "      ";
    os << itsCoeff << std::endl;
    return os;
}


template class TOrbitalImp<double>;
template class TOrbitalImp<std::complex<double> >;
