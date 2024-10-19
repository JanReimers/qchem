// File: ExactIrrepCD.C  Exact implementation of the charged density.


#include "Imp/Cluster/Molecule.H"
#include "Imp/ChargeDensity/IrrepCD.H"
#include "Imp/Hamiltonian/HamiltonianTerm.H"
#include "Imp/Fitting/FittedFunction.H"
#include <IntegralDataBase.H>
#include <QuantumNumber.H>
#include <ChargeDensity.H>
#include "oml/vector3d.h"
#include "oml/vector.h"
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include <cassert>
#include <complex>
#include <iostream>
#include <stdlib.h>

typedef Vector3D<std::complex<double> > Vec3;


RVec3  GradientContraction(const Vector<RVec3 >&, const Vector<double>&, const SMatrix<double>&);
RVec3  GradientContraction(const Vector<Vec3 >&, const Vector<std::complex<double> >&, const SMatrix<std::complex<double> >&);

//------------------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> IrrepCD<T>::IrrepCD()
{};

template <class T> IrrepCD<T>::IrrepCD(const DenSMat& theDensityMatrix,
                                                 const TIrrepBasisSet<T>* theBasisSet,
                                                 const Spin& s)
    : itsDensityMatrix(theDensityMatrix)
    , itsBasisSet(theBasisSet)
    , itsSpin(s)
{
    assert(itsBasisSet);
};

//-----------------------------------------------------------------------------
//
//  Total energy terms for a charge density.
//
template <> ChargeDensity::SMat IrrepCD<double>::GetRepulsion(const IrrepBasisSet* bs_ab) const
{
    const TIrrepBasisSet<double>* tbs_cd=dynamic_cast<const TIrrepBasisSet<double>*>(itsBasisSet);
    const TIrrepBasisSet<double>* tbs_ab=dynamic_cast<const TIrrepBasisSet<double>*>(bs_ab);
    return tbs_ab->GetRepulsion(itsDensityMatrix,tbs_cd);
}

template <> ChargeDensity::SMat IrrepCD<double>::GetExchange(const IrrepBasisSet* bs_ab) const
{
    const TIrrepBasisSet<double>* tbs_ab=dynamic_cast<const TIrrepBasisSet<double>*>(bs_ab);
    return tbs_ab->GetExchange(itsDensityMatrix,itsBasisSet);
}

//TODO: fix all complex fudges
// Fudge to get things to build.
template <> ChargeDensity::SMat IrrepCD<std::complex<double> >::GetRepulsion(const IrrepBasisSet* bs) const
{
    assert(itsBasisSet->GetID()==bs->GetID());
    assert(false);
    return SMat();
}

//template <class T> ChargeDensity::SMat ExactIrrepCD<T>::GetExchange(const BasisSet* bs) const
template <> ChargeDensity::SMat IrrepCD<std::complex<double> >::GetExchange(const IrrepBasisSet* bs) const
{
    assert(itsBasisSet->GetID()==bs->GetID());
 assert(false);
      return SMat();
}

template <class T> double IrrepCD<T>::GetEnergy(const HamiltonianTerm* v) const
{
    const HamiltonianTermImp* vi=dynamic_cast<const HamiltonianTermImp*>(v);
    assert(vi);
    T ComplexE=Sum(DirectMultiply(itsDensityMatrix,vi->GetCachedMatrix(itsBasisSet,itsSpin)));
    assert(fabs(imag(ComplexE))<1e-8);
    return real(ComplexE);
}


template <class T> double IrrepCD<T>::GetTotalCharge() const
{
    return real(Dot(itsDensityMatrix,itsBasisSet->GetOverlap()));
}
//------------------------------------------------------------------------------
//
//  Required by fitting routines.
//
template <class T> Vector<double> IrrepCD<T>::GetRepulsions(const IrrepBasisSet* fbs) const
{
    return itsBasisSet->GetRepulsion3C(itsDensityMatrix,fbs);
}


//-------------------------------------------------------------------------
//
//  SCF convergence stuff.
//
template <class T> void IrrepCD<T>::ReScale(double factor)
{
    itsDensityMatrix*=factor;
}

template <class T> void IrrepCD<T>::MixIn(const ChargeDensity& cd,double c)
{
    const IrrepCD<T>* eicd = dynamic_cast<const IrrepCD<T>*>(&cd);
    assert(eicd);
    assert(itsBasisSet->GetID() == eicd->itsBasisSet->GetID());
    itsDensityMatrix = itsDensityMatrix*(1-c) + eicd->itsDensityMatrix*c;
}

template <class T> double IrrepCD<T>::GetChangeFrom(const ChargeDensity& cd) const
{
    const IrrepCD<T>* eicd = dynamic_cast<const IrrepCD<T>*>(&cd);
    assert(eicd);
    assert(itsBasisSet->GetID() == eicd->itsBasisSet->GetID());
    return Max(fabs(itsDensityMatrix - eicd->itsDensityMatrix));
}

//-------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> void IrrepCD<T>::ShiftOrigin(const RVec3& newCenter)
{
    std::cerr << "ExactIrrepCD::ShiftOrigin this is an odd thing to do for an exact charge density" << std::endl;
    itsBasisSet=dynamic_cast<const TIrrepBasisSet<T>*>(itsBasisSet->Clone(newCenter));
}

template <class T> double IrrepCD<T>::operator()(const RVec3& r) const
{
    Vector<T> phir=(*itsBasisSet)(r);
    return real(phir*itsDensityMatrix*conj(phir));
}

template <class T> RVec3 IrrepCD<T>::Gradient(const RVec3& r) const
{
    Vector<T> phir=(*itsBasisSet)(r);
    Vector<RVec3 > gphir=itsBasisSet->Gradient(r);
    return GradientContraction(gphir,phir,itsDensityMatrix);
}

template <class T> void IrrepCD<T>::Eval(const Mesh& m, Vec& v) const
{
    index_t nm=v.size(), nb=itsBasisSet->GetNumFunctions();
    const Matrix<T>& ms((*itsBasisSet)(m));
    Vec::Subscriptor      vs(v);

    DenMat temp(nb,nm);
    Fill(temp,T(0));
    typename DenMat::Subscriptor ts(temp);

    for (index_t ib=1; ib<=nb; ib++)
        for (index_t jb=1; jb<=nb; jb++)
        {
            T dtemp=itsDensityMatrix(ib,jb);
            for (index_t iv=1; iv<=nm; iv++)
                ts(ib,iv)+=dtemp*conj(ms(jb,iv));
        }

    for (index_t ib=1; ib<=nb; ib++)
        for (index_t iv=1; iv<=nm; iv++)
            vs(iv)+=real(ms(ib,iv)*ts(ib,iv));
}

//-----------------------------------------------------------------------
//
//  Streamable stuff.
//
template <class T> std::ostream& IrrepCD<T>::Write(std::ostream& os) const
{
    return os << itsDensityMatrix;
}

template <class T> std::istream& IrrepCD<T>::Read(std::istream& is)
{
    is >> itsDensityMatrix;

//    TIrrepBasisSet<T>* tbs = dynamic_cast<TIrrepBasisSet<T>*>(IrrepBasisSet::Factory(is));
//    assert(tbs);
//    is >> *tbs;
//    itsBasisSet.reset(tbs);

    return is;
}

template class IrrepCD<double>;
//template class ExactIrrepCD<std::complex<double> >;




RVec3 GradientContraction(const Vector<RVec3>& g, const Vector<double>& v, const SMatrix<double>& m)
{
    assert(m.GetNumRows()==m.GetNumCols());
    assert(v.size      ()==m.GetNumCols());
    assert(g.size      ()==m.GetNumCols());

    RVec3 ret(0,0,0);
    for (unsigned int i=1; i<=v.size(); i++)
        for (unsigned int j=1; j<=v.size(); j++)
            ret+=m(i,j)*(g(i)*v(j)+v(i)*g(j));
    return ret;
}

RVec3 GradientContraction(const Vector<Vec3>& g, const Vector<std::complex<double> >& v, const SMatrix<std::complex<double> >& m)
{
    assert(m.GetNumRows()==m.GetNumCols());
    assert(v.size      ()==m.GetNumCols());
    assert(g.size      ()==m.GetNumCols());

    Vec3 ret(0,0,0);
    for (unsigned int i=1; i<=v.size(); i++)
        for (unsigned int j=1; j<=v.size(); j++)
            ret+=m(i,j)*(g(i)*conj(v(j))+v(i)*conj(g(j)));
    return real(ret);
}

