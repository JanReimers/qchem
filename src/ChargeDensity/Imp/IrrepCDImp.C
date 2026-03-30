// File: ExactIrrepCD.C  Exact implementation of the charged density.
module;
#include <cassert>
#include <complex>
#include <iostream>
#include <stdlib.h>
#include <vector>
#include "blaze/Math.h" 

module qchem.ChargeDensity.Imp.IrrepCD;
import qchem.Orbital_DFT_IBS;
import qchem.Symmetry;
import qchem.Molecule;
import qchem.Blaze;

typedef Vector3D<std::complex<double> > Vec3;

rvec3_t  GradientContraction(const vec_t<rvec3_t >&, const vec_t<double>&, const rsmat_t&);
// rvec3_t  GradientContraction(const Vector<Vec3 >&, const Vector<std::complex<double> >&, const smat_t<std::complex<double> >&);

//------------------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> IrrepCD<T>::IrrepCD()
{};

template <class T> IrrepCD<T>::IrrepCD
    (const DenSMat& D,
    const Orbital_IBS<T>* theBasisSet,
    Irrep_QNs qns)
    : itsDensityMatrix(D)
    , itsBasisSet(theBasisSet)
    , itsSpin(qns.ms)
    , itsIrrep(qns)
{
    assert(itsBasisSet);
};

template <> bool IrrepCD<double>::IsZero() const
{
    // return max(abs(itsDensityMatrix))==0.0;
    return isZero(itsDensityMatrix);
}


//-----------------------------------------------------------------------------
//
//  Total energy terms for a charge density.
//
template <> smat_t<double> IrrepCD<double>::GetRepulsion(const Orbital_HF_IBS<double>* bs_ab) const
{
    if (IsZero()) return zero<double>(bs_ab->GetNumFunctions());
    const Orbital_HF_IBS<double>* bs_cd=dynamic_cast<const Orbital_HF_IBS<double>*>(itsBasisSet);
    assert(bs_cd);
    return bs_ab->Direct(itsDensityMatrix,bs_cd);
}

template <> smat_t<double> IrrepCD<double>::GetExchange(const Orbital_HF_IBS<double>* bs_ab) const
{
    if (IsZero()) return zero<double>(bs_ab->GetNumFunctions());
    const Orbital_HF_IBS<double>* bs_cd=dynamic_cast<const Orbital_HF_IBS<double>*>(itsBasisSet);
    assert(bs_cd);
    return bs_ab->Exchange(itsDensityMatrix,bs_cd);
}

//------------------------------------------------------------------------------
//
//  Required by fitting routines.
//
template <class T> rvec_t IrrepCD<T>::GetRepulsion3C(const Fit_IBS* fbs) const
{
    if (IsZero()) return rvec_t(fbs->GetNumFunctions(),0.0);
    auto dftbs=dynamic_cast<const Orbital_DFT_IBS<T>*>(itsBasisSet);
    assert(dftbs);
    return dftbs->Repulsion3C(itsDensityMatrix,fbs);
}


template <class T> double IrrepCD<T>::DM_Contract(const Static_CC* v) const
{
    T ComplexE=sum(itsDensityMatrix % v->GetMatrix(itsBasisSet,itsSpin)); //% is the blaze op for the Shur (direct) product.
    assert(fabs(std::imag(ComplexE))<1e-8);
    return std::real(ComplexE);
}

template <class T> double IrrepCD<T>::DM_Contract(const Dynamic_CC* v,const DM_CD* cd) const
{
    T ComplexE=sum(itsDensityMatrix % v->GetMatrix(itsBasisSet,itsSpin,cd)); //% is the blaze op for the Shur (direct) product.
    assert(fabs(std::imag(ComplexE))<1e-8);
    return std::real(ComplexE);
}

template <class T> double IrrepCD<T>::GetTotalCharge() const
{
    return std::real(sum(itsDensityMatrix%itsBasisSet->Overlap())); //% is the blaze op for the Shur (direct) product.
}


//-------------------------------------------------------------------------
//
//  SCF convergence stuff.
//
template <class T> void IrrepCD<T>::ReScale(double factor)
{
    // No UT coverage
    itsDensityMatrix*=factor;
}

template <class T> void IrrepCD<T>::MixIn(const DM_CD& cd,double c)
{
    const IrrepCD<T>* eicd = dynamic_cast<const IrrepCD<T>*>(&cd);
    assert(eicd);
    assert(itsBasisSet->GetID() == eicd->itsBasisSet->GetID());
    itsDensityMatrix = itsDensityMatrix*(1-c) + eicd->itsDensityMatrix*c;
}

template <class T> double IrrepCD<T>::GetChangeFrom(const DM_CD& cd) const
{
    const IrrepCD<T>* eicd = dynamic_cast<const IrrepCD<T>*>(&cd);
    assert(eicd);
    assert(itsBasisSet->GetID() == eicd->itsBasisSet->GetID());
    return norm(itsDensityMatrix - eicd->itsDensityMatrix);
}

//-------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> double IrrepCD<T>::operator()(const rvec3_t& r) const
{
    vec_t<T> phir=(*itsBasisSet)(r);
    return std::real(trans(phir)*itsDensityMatrix*conj(phir));
}

template <class T> rvec3_t IrrepCD<T>::Gradient(const rvec3_t& r) const
{
    // No UT coverage
    vec_t<T> phir=(*itsBasisSet)(r);
    vec_t<rvec3_t > gphir=itsBasisSet->Gradient(r);
    return GradientContraction(gphir,phir,itsDensityMatrix);
}


//-----------------------------------------------------------------------
//
//  Streamable stuff.
//
template <class T> std::ostream& IrrepCD<T>::Write(std::ostream& os) const
{
    return os << itsDensityMatrix;
}

template class IrrepCD<double>;
//template class ExactIrrepCD<std::complex<double> >;




rvec3_t GradientContraction(const vec_t<rvec3_t>& g, const vec_t<double>& v, const rsmat_t& m)
{
    // No UT coverage
    assert(v.size      ()==m.columns());
    assert(g.size      ()==m.columns());

    rvec3_t ret(0,0,0);
    for (unsigned int i=0; i<v.size(); i++)
        for (unsigned int j=0; j<v.size(); j++)
            ret+=m(i,j)*(g[i]*v[j]+v[i]*g[j]);
    return ret;
}

// rvec3_t GradientContraction(const Vector<Vec3>& g, const Vector<std::complex<double> >& v, const smat_t<std::complex<double> >& m)
// {
//     // No UT coverage
//     assert(v.size      ()==m.columns());
//     assert(g.size      ()==m.columns());

//     Vec3 ret(0,0,0);
//     for (unsigned int i=1; i<=v.size(); i++)
//         for (unsigned int j=1; j<=v.size(); j++)
//             ret+=m(i-1,j-1)*(g(i)*conj(v(j))+v(i)*conj(g(j)));
//     return real(ret);
// }

