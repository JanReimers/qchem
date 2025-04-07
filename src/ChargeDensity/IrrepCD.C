// File: ExactIrrepCD.C  Exact implementation of the charged density.


#include "Imp/Cluster/Molecule.H"
#include "Imp/ChargeDensity/IrrepCD.H"
#include "Imp/Hamiltonian/HamiltonianTerm.H"
#include "Imp/Fitting/FittedFunction.H"
#include <HF_IBS.H>
#include <Fit_IBS.H>
#include <DFT_IBS.H>
#include <Symmetry.H>
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

template <class T> IrrepCD<T>::IrrepCD
    (const DenSMat& D,
    const TOrbital_IBS<T>* theBasisSet,
    Irrep_QNs qns)
    : itsDensityMatrix(D)
    , itsBasisSet(theBasisSet)
    , itsSpin(qns.ms)
    , itsQNs(qns)
{
    assert(itsBasisSet);
};

template <> bool IrrepCD<double>::IsZero() const
{
    return Max(fabs(itsDensityMatrix))==0.0;
}

template <> DM_CD::SMat IrrepCD<double>::ZeroM(size_t N) const
{
    SMat S(N);
    Fill(S,0.0);
    return S;
}

template <> IrrepCD<double>::RVec IrrepCD<double>::ZeroV(size_t N) const
{
    RVec V(N);
    Fill(V,0.0);
    return V;
}
//-----------------------------------------------------------------------------
//
//  Total energy terms for a charge density.
//
template <> DM_CD::SMat IrrepCD<double>::GetRepulsion(const TOrbital_HF_IBS<double>* bs_ab) const
{
    if (IsZero()) return ZeroM(bs_ab->size());
    const TOrbital_HF_IBS<double>* bs_cd=dynamic_cast<const TOrbital_HF_IBS<double>*>(itsBasisSet);
    assert(bs_cd);
    return bs_ab->Direct(itsDensityMatrix,bs_cd);
}

template <> DM_CD::SMat IrrepCD<double>::GetExchange(const TOrbital_HF_IBS<double>* bs_ab) const
{
    if (IsZero()) return ZeroM(bs_ab->size());
    const TOrbital_HF_IBS<double>* bs_cd=dynamic_cast<const TOrbital_HF_IBS<double>*>(itsBasisSet);
    assert(bs_cd);
    return bs_ab->Exchange(itsDensityMatrix,bs_cd);
}

//------------------------------------------------------------------------------
//
//  Required by fitting routines.
//
template <class T> Vector<double> IrrepCD<T>::GetRepulsion3C(const Fit_IBS* fbs) const
{
    if (IsZero()) return ZeroV(fbs->size());
    auto dftbs=dynamic_cast<const TOrbital_DFT_IBS<T>*>(itsBasisSet);
    assert(dftbs);
    return dftbs->Repulsion3C(itsDensityMatrix,fbs);
}


template <class T> double IrrepCD<T>::DM_Contract(const Static_HT* v) const
{
    T ComplexE=Sum(DirectMultiply(itsDensityMatrix,v->GetMatrix(itsBasisSet,itsSpin)));
    assert(fabs(imag(ComplexE))<1e-8);
    return real(ComplexE);
}

template <class T> double IrrepCD<T>::DM_Contract(const Dynamic_HT* v,const DM_CD* cd) const
{
    T ComplexE=Sum(DirectMultiply(itsDensityMatrix,v->GetMatrix(itsBasisSet,itsSpin,cd)));
    assert(fabs(imag(ComplexE))<1e-8);
    return real(ComplexE);
}

template <class T> double IrrepCD<T>::GetTotalCharge() const
{
    
    //std::cout << "D=" << itsDensityMatrix << " S=" << itsBasisSet->GetOverlap() << std::endl;
    // int N=itsDensityMatrix.GetNumRows();
    // assert(N%2==0);
    // int NL=N/2;
    // SMat S=itsBasisSet->GetOverlap();
    // SMat DLL=itsDensityMatrix.SubMatrix(MatLimits(1,NL,1,NL));
    // SMat DSS=itsDensityMatrix.SubMatrix(MatLimits(NL+1,N, NL+1,N));
    // SMat SLL=S.SubMatrix(MatLimits(1,NL,1,NL));
    // SMat SSS=S.SubMatrix(MatLimits(NL+1,N, NL+1,N));
    // std::cout.precision(10);
    // std::cout << "Charge LL=" << real(Dot(DLL,SLL)) << " SS=" << real(Dot(DSS,SSS)) << std::endl;
    return real(Dot(itsDensityMatrix,itsBasisSet->Overlap()));
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
    return Max(fabs(itsDensityMatrix - eicd->itsDensityMatrix));
}

//-------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> void IrrepCD<T>::ShiftOrigin(const RVec3& newCenter)
{
    // No UT coverage
    std::cerr << "ExactIrrepCD::ShiftOrigin this is an odd thing to do for an exact charge density" << std::endl;
    itsBasisSet=dynamic_cast<const TOrbital_IBS<T>*>(itsBasisSet->Clone(newCenter));
}

template <class T> double IrrepCD<T>::operator()(const RVec3& r) const
{
    Vector<T> phir=(*itsBasisSet)(r);
    return real(phir*itsDensityMatrix*conj(phir));
}

template <class T> RVec3 IrrepCD<T>::Gradient(const RVec3& r) const
{
    // No UT coverage
    Vector<T> phir=(*itsBasisSet)(r);
    Vector<RVec3 > gphir=itsBasisSet->Gradient(r);
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
    // No UT coverage
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
    // No UT coverage
    assert(m.GetNumRows()==m.GetNumCols());
    assert(v.size      ()==m.GetNumCols());
    assert(g.size      ()==m.GetNumCols());

    Vec3 ret(0,0,0);
    for (unsigned int i=1; i<=v.size(); i++)
        for (unsigned int j=1; j<=v.size(); j++)
            ret+=m(i,j)*(g(i)*conj(v(j))+v(i)*conj(g(j)));
    return real(ret);
}

