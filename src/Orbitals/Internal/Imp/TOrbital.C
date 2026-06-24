// File: OrbitalImp.C  Implementation of an orbital.
module;
#include <iostream>
#include <iomanip>
#include <complex>
#include <cassert>
module qchem.Orbitals.Internal.OrbitalImp;
import qchem.Blaze;

namespace qchem::Orbitals
{

//-----------------------------------------------------------------
//
//  Construction zone
//
template <class T> TOrbitalImp<T>::
TOrbitalImp(const tobs_t<T>* bs,const vec_t<T>& _C,const vec_t<T>& _CPrime,double e, const Orbital_QNs& qns)
    : OrbitalImp   (e,qns)
    , itsCoeff     (_C)
    , itsCoeffPrime(_CPrime)
    , itsBasisSet  (bs)
{
};

template <class T> void TOrbitalImp<T>::AddDensityMatrix(hmat_t<T>& D, hmat_t<T>& DPrime) const
{
    // D += occ * |C><C|.  The outer product C_i conj(C_j) is HERMITIAN (= symmetric for real T), so it
    // fills a HermitianMatrix; the diagonal occ*|C_i|^2 is real.
    if (IsOccupied())
    {
        hmat_t<T> CCd=blazem::outer(itsCoeff,blazem::conj(itsCoeff))*GetOccupation();
        D+=CCd;
        hmat_t<T> CCd_prime=blazem::outer(itsCoeffPrime,blazem::conj(itsCoeffPrime))*GetOccupation();
        DPrime+=CCd_prime;
    }
}



//----------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> T TOrbitalImp<T>::operator()(const rvec3_t& r) const
{
    vec_t<T> gr=(*itsBasisSet)(r);   // was hardcoded rvec_t (basis values are complex for the PW case)
    assert(gr.size()==itsCoeff.size());
    return blazem::trans(itsCoeff) * gr;
}

//BUG
template <class T> vec3_t<T> TOrbitalImp<T>::Gradient(const rvec3_t& r) const
{
    vec3_t<T> ret(0,0,0);
    vec3vec_t<T> grads=itsBasisSet->Gradient(r);
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
template class TOrbitalImp<dcmplx>;

} //namespace