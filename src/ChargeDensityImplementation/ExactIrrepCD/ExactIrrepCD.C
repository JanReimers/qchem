// File: ExactIrrepCD.C  Exact implementation of the charged density.


#include "Cluster/Molecule.H"
#include "IntegralDataBase.H"
#include "QuantumNumber.H"
#include "ChargeDensity.H"
#include "ChargeDensityImplementation/ExactIrrepCD/ExactIrrepCD.H"
#include "HamiltonianImplementation/HamiltonianTermImplementation.H"
#include "FunctionsImp/FittedFunctionImplementation.H"
#include "oml/vector3d.h"
#include "oml/vector.h"
#include "oml/matrix.h"
#include "oml/smatrix.h"
#include <cassert>
#include <complex>
#include <iostream>
#include <stdlib.h>

typedef Vector3D<std::complex<double> > Vec3;

double FastContraction(const Vector<double>&, const SMatrix<double>&);
RVec3  FastContraction(const Vector<RVec3 >&, const Vector<double>&, const SMatrix<double>&);

double FastContraction(const Vector<std::complex<double> >&, const SMatrix<std::complex<double> >&);
RVec3  FastContraction(const Vector<Vec3 >&, const Vector<std::complex<double> >&, const SMatrix<std::complex<double> >&);

//------------------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> ExactIrrepCD<T>::ExactIrrepCD()
{};

template <class T> ExactIrrepCD<T>::ExactIrrepCD(const DenSMat& theDensityMatrix,
                                                 const rc_ptr<const IrrepBasisSet>& theBasisSet,
                                                 const Spin& s)
    : itsDensityMatrix(theDensityMatrix)
    , itsBasisSet(theBasisSet)
    , itsCastedBasisSet(dynamic_cast<const TIrrepBasisSet<T>*>(theBasisSet.get()))
    , itsSpin(s)
{
    assert(itsCastedBasisSet);
};

//-----------------------------------------------------------------------------
//
//  Totale energy terms for a charge density.
//
template <class T> ChargeDensity::SMat ExactIrrepCD<T>::GetOverlap  (const IrrepBasisSet*) const
{
    std::cerr << "ExactIrrepCD::GetOverlap 4 electron integrals not implementated yet" << std::endl;
    exit(-1);
    return SMat();
}

#include "Imp/BasisSet/SphericalGaussian/QuantumNumber.H"

//template <class T> ChargeDensity::SMat ExactIrrepCD<T>::GetRepulsion(const BasisSet* bs) const
template <> ChargeDensity::SMat ExactIrrepCD<double>::GetRepulsion(const IrrepBasisSet* bs_ab) const
{
//    assert(itsBasisSet->GetID()==bs->GetID()); basis sets get cloned, so this won't work.
//    std::cout << "   ExactIrrepCD GetRepulsion Lab=" << bs_ab->GetQuantumNumber() << ", Lcd=" << itsBasisSet->GetQuantumNumber() << std::endl;
    const TIrrepBasisSet<double>* tbs_cd=dynamic_cast<const TIrrepBasisSet<double>*>(itsBasisSet.get());
    const TIrrepBasisSet<double>* tbs_ab=dynamic_cast<const TIrrepBasisSet<double>*>(bs_ab);
    SMat Jab= tbs_ab->GetRepulsion(itsDensityMatrix,tbs_cd);
//    const SphericalSymmetryQN* qn_ab=dynamic_cast<const SphericalSymmetryQN*>(&bs_ab->GetQuantumNumber());
//    const SphericalSymmetryQN* qn_cd=dynamic_cast<const SphericalSymmetryQN*>(&itsBasisSet->GetQuantumNumber());
//    if (qn_ab->GetL()==1 && qn_cd->GetL()==1)
//    {
//    std::cout.setf(std::ios::fixed,std::ios::floatfield);
//
//    std::cout << "     Dcd=" << itsDensityMatrix << std::endl;
//    std::cout << "     Jab=" << Jab << std::endl;
//
//    }

    return Jab;
}

//template <class T> ChargeDensity::SMat ExactIrrepCD<T>::GetExchange(const BasisSet* bs) const
template <> ChargeDensity::SMat ExactIrrepCD<double>::GetExchange(const IrrepBasisSet* bs_ab) const
{
//    assert(itsBasisSet->GetID()==bs->GetID());
    const TIrrepBasisSet<double>* tbs_ab=dynamic_cast<const TIrrepBasisSet<double>*>(bs_ab);
    SMat Kab= tbs_ab->GetExchange(itsDensityMatrix,itsCastedBasisSet);
//    const SphericalSymmetryQN* qn_ab=dynamic_cast<const SphericalSymmetryQN*>(&bs_ab->GetQuantumNumber());
//    const SphericalSymmetryQN* qn_cd=dynamic_cast<const SphericalSymmetryQN*>(&itsBasisSet->GetQuantumNumber());
//    if (qn_ab->GetL() != qn_cd->GetL())
//    {
//        std::cout.precision(3);
//        std::cout.width(5);
//        std::cout << "     Lab=" << qn_ab->GetL() << " Lab=" << qn_cd->GetL() << std::endl;
//        std::cout << "     Dcd=" << itsDensityMatrix << std::endl;
//        std::cout << "     Kab=" << Kab << std::endl;
//
//    }
    return Kab;
//    return tbs_ab->GetExchange(itsDensityMatrix,itsCastedBasisSet);
}

//TODO: fix all complex fudges
// Fudge to get things to build.
template <> ChargeDensity::SMat ExactIrrepCD<std::complex<double> >::GetRepulsion(const IrrepBasisSet* bs) const
{
    assert(itsBasisSet->GetID()==bs->GetID());
    return SMat();
}

//template <class T> ChargeDensity::SMat ExactIrrepCD<T>::GetExchange(const BasisSet* bs) const
template <> ChargeDensity::SMat ExactIrrepCD<std::complex<double> >::GetExchange(const IrrepBasisSet* bs) const
{
    assert(itsBasisSet->GetID()==bs->GetID());
    return SMat();
}

template <class T> double ExactIrrepCD<T>::GetEnergy(const HamiltonianTerm* v) const
{
    const HamiltonianTermImplementation* vi=dynamic_cast<const HamiltonianTermImplementation*>(v);
    assert(vi);
    T ComplexE=Sum(DirectMultiply(itsDensityMatrix,vi->GetCachedMatrix(itsBasisSet.get(),itsSpin)));
    assert(fabs(imag(ComplexE))<1e-8);
    return real(ComplexE);
}

//template <class T> double ExactIrrepCD<T>::GetCoulombAtOrigin() const
//{
//    Molecule cl;
//    cl.Insert(new Atom(1,0,RVec3(0,0,0)));
//    return -real(Dot(itsDensityMatrix,itsCastedBasisSet->GetDataBase()->GetNuclear(cl)));
//}

template <class T> double ExactIrrepCD<T>::GetTotalCharge() const
{
    return real(Dot(itsDensityMatrix,itsBasisSet->GetOverlap()));
}
//------------------------------------------------------------------------------
//
//  Required by fitting routines.
//

template <class T> void ExactIrrepCD<T>::InjectOverlaps  (FittedFunction* ff, const IrrepBasisSet* fbs) const
{
    RVec delta_ff=itsBasisSet->GetOverlap3C(itsDensityMatrix,fbs);

    FittedFunctionImplementation<T>* ffi=dynamic_cast<FittedFunctionImplementation<T>*>(ff);
    assert(ffi);
    ffi->GetFitCoeff()+=delta_ff;
}

template <class T> void ExactIrrepCD<T>::InjectRepulsions(FittedFunction* ff, const IrrepBasisSet* fbs) const
{
    RVec delta_ff=itsBasisSet->GetRepulsion3C(itsDensityMatrix,fbs);

    FittedFunctionImplementation<T>* ffi=dynamic_cast<FittedFunctionImplementation<T>*>(ff);
    assert(ffi);
    ffi->GetFitCoeff()+=delta_ff;
}

//-------------------------------------------------------------------------
//
//  SCF convergence stuff.
//
template <class T> void ExactIrrepCD<T>::ReScale(double factor)
{
    itsDensityMatrix*=factor;
}

template <class T> void ExactIrrepCD<T>::MixIn(const ChargeDensity& cd,double c)
{
    const ExactIrrepCD<T>* eicd = dynamic_cast<const ExactIrrepCD<T>*>(&cd);
    assert(eicd);
    assert(itsBasisSet->GetID() == eicd->itsBasisSet->GetID());
    itsDensityMatrix = itsDensityMatrix*(1-c) + eicd->itsDensityMatrix*c;
}

template <class T> double ExactIrrepCD<T>::GetChangeFrom(const ChargeDensity& cd) const
{
    const ExactIrrepCD<T>* eicd = dynamic_cast<const ExactIrrepCD<T>*>(&cd);
    assert(eicd);
    assert(itsBasisSet->GetID() == eicd->itsBasisSet->GetID());
    return Max(fabs(itsDensityMatrix - eicd->itsDensityMatrix));
}

//-------------------------------------------------------------------------
//
//  Real space function stuff.
//
template <class T> void ExactIrrepCD<T>::ShiftOrigin(const RVec3& newCenter)
{
    std::cerr << "ExactIrrepCD::ShiftOrigin this is an odd thing to do for an exact charge density" << std::endl;
    itsBasisSet.reset(itsBasisSet->Clone(newCenter));
}

template <class T> double ExactIrrepCD<T>::operator()(const RVec3& r) const
{
    return FastContraction((*itsCastedBasisSet)(r),itsDensityMatrix);
}

template <class T> RVec3 ExactIrrepCD<T>::Gradient(const RVec3& r) const
{
    return FastContraction(itsCastedBasisSet->Gradient(r),(*itsCastedBasisSet)(r),itsDensityMatrix);
}

template <class T> void ExactIrrepCD<T>::Eval(const Mesh& m, Vec& v) const
{
    index_t nm=v.size(), nb=itsBasisSet->GetNumFunctions();
    const Matrix<T>& ms((*itsCastedBasisSet)(m));
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
template <class T> std::ostream& ExactIrrepCD<T>::Write(std::ostream& os) const
{
    return os << itsDensityMatrix << *itsBasisSet;
}

template <class T> std::istream& ExactIrrepCD<T>::Read(std::istream& is)
{
    is >> itsDensityMatrix;

    TIrrepBasisSet<T>* tbs = dynamic_cast<TIrrepBasisSet<T>*>(IrrepBasisSet::Factory(is));
    assert(tbs);
    is >> *tbs;
    itsBasisSet.reset(tbs);
    itsCastedBasisSet=tbs;

    return is;
}

template class ExactIrrepCD<double>;
//template class ExactIrrepCD<std::complex<double> >;

















double FastContraction(const Vector<double>& v, const SMatrix<double>& m)
{
    assert(m.GetNumRows()==m.GetNumCols());
    assert(v.size      ()==m.GetNumCols());
    double ret=0;
    for (unsigned int i=1; i<=v.size(); i++)
        for (unsigned int j=1; j<=v.size(); j++)
            ret+=v(i)*m(i,j)*v(j);
    return ret;
}

double FastContraction(const Vector<std::complex<double> >& v, const SMatrix<std::complex<double> >& m)
{
    assert(m.GetNumRows()==m.GetNumCols());
    assert(v.size      ()==m.GetNumCols());
    double ret=0;
    for (unsigned int i=1; i<=v.size(); i++)
        for (unsigned int j=1; j<=v.size(); j++)
            ret+=real(v(i)*m(i,j)*conj(v(j)));
    return ret;
}

RVec3 FastContraction(const Vector<RVec3>& g, const Vector<double>& v, const SMatrix<double>& m)
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

RVec3 FastContraction(const Vector<Vec3>& g, const Vector<std::complex<double> >& v, const SMatrix<std::complex<double> >& m)
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

