// File: TBasisSetImplementation.C


#include "BasisSetImplementation/TBasisSetImplementation.H"
#include "BasisSet.H"
#include "QuantumNumber.H"
#include "NumericalIE.H"
#include "AnalyticIE.H"
#include "IntegralDataBase.H"
#include "LASolver/LASolver.H"
#include "Hamiltonian.H"
#include "Mesh/Mesh.H"

#include "Misc/ERI4.H"
#include "OrbitalImplementation/TOrbitalGroupImplementation.H"
#include "FunctionsImp/FittedFunctionImplementation.H"
#include "ChargeDensityImplementation/ExactIrrepCD/ExactIrrepCD.H"
#include <cassert>
#include <iostream>

//-----------------------------------------------------------------------------
//
//  Construction zone
//
template <class T> TBasisSetImplementation<T>::TBasisSetImplementation()
    : VectorFunctionBuffer<T>(false,false) //don't pickle scalar or gradient.
    , itsNumericalIE(0)
    , itsDataBase      ( )
    , itsLASolver      (0)
{};

template <class T> TBasisSetImplementation<T>::TBasisSetImplementation(const LinearAlgebraParams& lap,IntegralDataBase<T>* theDataBase)
    : VectorFunctionBuffer<T>(false,false) //don't pickle scalar or gradient.
    , itsLAParams      (lap)
    , itsNumericalIE(0)
    , itsDataBase      (theDataBase)
    , itsLASolver      (0)
{};

template <class T> TBasisSetImplementation<T>::TBasisSetImplementation(const TBasisSetImplementation<T>& bs)
    : VectorFunctionBuffer<T>(false,false) //don't pickle scalar or gradient.
    , itsLAParams      (bs.itsLAParams)
    , itsNumericalIE   (bs.itsNumericalIE)
    , itsDataBase      (bs.itsDataBase)
    , itsLASolver      (0)
{};

template <class T> TBasisSetImplementation<T>::~TBasisSetImplementation()
{
    delete itsLASolver;
}

//-----------------------------------------------------------------------------
//
//  Post construction initializations called by dervied classes.
//
template <class T> void TBasisSetImplementation<T>::Insert(NumericalIE<T>* ie)
{
    assert(ie);
    itsNumericalIE.reset(ie);
    itsNumericalIE->Insert(this);
    itsDataBase->Insert(ie);
}

template <class T> void TBasisSetImplementation<T>::Insert(AnalyticIE<T>* ie)
{
    assert(ie);
    itsAnalyticIE.reset(ie);
    itsDataBase->Insert(ie);
    RVec ns=ie->MakeNormalization();
    RVec cs=ie->MakeCharge(this);
    int i=1;
    for (auto bf:*this) 
    {
        bf->Init(ns(i),cs(i));
        i++;
    }
}

template <class T> IntegralDataBase<T>* TBasisSetImplementation<T>::GetDataBase() const
{
    assert(&*itsDataBase);
    return itsDataBase;
}
template <class T> AnalyticIE<T>* TBasisSetImplementation<T>::GetAnalyticIE() const
{
    assert(&*itsAnalyticIE);
    return &*itsAnalyticIE;
}

 

//-----------------------------------------------------------------------------
//
//  Build up the hamiltonian matrix, and generate new orbitals from the eigen
//  vectors.  The hamiltonian controls which terms are included in H.
//
// TODO: Why do we need to pass in const rc_ptr<const BasisSet>& rc ????
//
template <class T> OrbitalGroup* TBasisSetImplementation<T>::
CreateOrbitals(const rc_ptr<const IrrepBasisSet>& rc,const Hamiltonian* ham, const Spin&S) const
{
    SMat H=ham->BuildHamiltonian(this,S);
    assert(!isnan(H));
    LASolver<T>* es=GetLASolver();
    assert(es);
    auto [U,e]=es->Solve(H);
    return new
           TOrbitalGroupImplementation<T>(rc,U,e,S);
}

template <class T> IrrepBasisSet::SMat TBasisSetImplementation<T>::
GetKinetic() const
{
    return GetDataBase()->GetKinetic(this);
}
template <class T> IrrepBasisSet::SMat TBasisSetImplementation<T>::
GetNuclear(const Cluster* cl) const
{
    return itsDataBase->GetNuclear(this,*cl);
}

template <class T> IrrepBasisSet::SMat TBasisSetImplementation<T>::
GetOverlap  (const FittedFunction* ff) const
{
    int n=this->GetNumFunctions();
    SMat J(n,n);
    Fill(J,0.0);
    const FittedFunctionImplementation<T>* ffi=dynamic_cast<const FittedFunctionImplementation<T>*>(ff);
    assert(ffi);
    assert(!isnan(ffi->itsFitCoeff));
    const ERI3& overlap=GetDataBase()->GetOverlap3C(this,ffi->CastBasisSet());
    typename Vector<T>::const_iterator f(ffi->itsFitCoeff.begin());
    for(index_t i=0; f!=ffi->itsFitCoeff.end(); f++,i++)
    {
        assert(!isnan(overlap[i]));
        SMat fo=(*f) * overlap[i];
        J+=fo;
    }
    assert(!isnan(J));
    return J;
}

template <class T> IrrepBasisSet::SMat TBasisSetImplementation<T>::
GetRepulsion(const FittedFunction* ff) const
{
    int n=this->GetNumFunctions();
    SMat J(n,n);
    Fill(J,0.0);
    const FittedFunctionImplementation<T>* ffi=dynamic_cast<const FittedFunctionImplementation<T>*>(ff);
    assert(ffi);
    const ERI3& repulsion=GetDataBase()->GetRepulsion3C(this,ffi->CastBasisSet());
    typename Vector<T>::const_iterator f(ffi->itsFitCoeff.begin());
    for(index_t i=0; f!=ffi->itsFitCoeff.end(); f++,i++) J+=SMat((*f) * repulsion[i]);
    assert(!isnan(J));
    return J;
}
#include "BasisSetImplementation/SphericalGaussian/SphericalSymmetryQN.H"
#include "BasisSetImplementation/BasisSetImplementation.H"
#include "Misc/DFTDefines.H"
#include "BasisSetImplementation/PolarizedGaussian/PolarizedGaussianBF.H"
#include "BasisSetImplementation/PolarizedGaussian/Gaussian/GaussianRF.H"


template <class T> IrrepBasisSet::SMat TBasisSetImplementation<T>::
GetRepulsion(const SMat& Dcd, const TIrrepBasisSet<T>* bs_cd) const
{
    assert(!isnan(Dcd));
    assert(itsBasisGroup);
//    std::cout << "    TBasisSetImplementation::GetRep Dcd=" << Dcd << std::endl;
//    const BasisSetImplementation* bsi=dynamic_cast<const BasisSetImplementation*>(this);
    BasisGroup::iecv_t iecs=itsBasisGroup->Flatten();
    const ERI4& Jfull=GetDataBase()->GetRepulsion4C(iecs);
    ERI4view J(Jfull,this->GetStartIndex(),bs_cd->GetStartIndex());
    int Nab=this->GetNumFunctions();
    int Ncd=bs_cd->GetNumFunctions();

    SMat Jab(Nab,Nab);
    for (int ia=1; ia<=Nab; ia++)
        for (int ib=ia; ib<=Nab; ib++)
        {
            T Jab_temp=0;
            for (int ic=1; ic<=Ncd; ic++)
                for (int id=1; id<=Ncd; id++) //Possible symmetric optimization here.  Need to be careful to handle complex D and ERIs.
                {
                    //std::cout << ia << " " << ib << " " << ic << " " << id << " " << J.GetIndex(ia,ib,ic,id) << " " << eris.GetIndex(ia,ib,ic,id) << " " << J(ia,ib,ic,id) << " " << eris(ia,ib,ic,id) << std::endl;
//                    assert(J(ia,ib,ic,id)==eris(ia,ib,ic,id));
                    Jab_temp+=J(ia,ib,ic,id)*Dcd(ic,id);
                }
            Jab(ia,ib)=Jab_temp;
        }
    assert(!isnan(Jab));

    return Jab;
}

template <class T> IrrepBasisSet::SMat TBasisSetImplementation<T>::
GetExchange(const SMat& Dcd, const TIrrepBasisSet<T>* bs_cd) const
{
    assert(!isnan(Dcd));
    BasisGroup::iecv_t iecs=itsBasisGroup->Flatten();
    const ERI4& Kfull=GetDataBase()->GetExchange4C(iecs);
    ERI4view K(Kfull,this->GetStartIndex(),bs_cd->GetStartIndex());
    int Nab=this->GetNumFunctions();
    int Ncd=bs_cd->GetNumFunctions();

    SMat Kab(Nab,Nab);
    for (int ia=1; ia<=Nab; ia++)
        for (int ib=ia; ib<=Nab; ib++)
        {
            T Kab_temp=0;
            for (int ic=1; ic<=Ncd; ic++)
                for (int id=1; id<=Ncd; id++) //Possible symmetric optimization here.  Need to be careful to handle complex D and ERIs.
                {
//                    std::cout << ia << " " << ib << " " << ic << " " << id << " " 
//                    << K.GetIndex(ia,ib,ic,id) << " " << eris.GetIndex(ia,ib,ic,id) 
//                    << " " << K(ia,ib,ic,id) << " " << eris(ia,ib,ic,id) << std::endl;
//                    assert(K(ia,ib,ic,id)==eris(ia,ib,ic,id));
//                    std::cout << ia << " " << ib << " " << ic << " " << id << " " 
//                    << K.GetIndex(ia,id,ic,ib) << " " << eris.GetIndex(ia,id,ic,ib) 
//                    << " " << K.Exchange(ia,id,ic,ib) << " " << eris.Exchange(ia,id,ic,ib) << std::endl;
//                    assert(K.Exchange(ia,id,ic,ib)==eris.Exchange(ia,id,ic,ib));
                    Kab_temp+=K.Exchange(ia,id,ic,ib)*Dcd(ic,id);
                }
            Kab(ia,ib)=Kab_temp;
        }
    assert(!isnan(Kab));

    return Kab;
}
//
//  Charge density repulsion calculations.
//
template <class T> double TBasisSetImplementation<T>::
GetCDRepulsion(const ChargeDensity* cd, const FittedFunction* ff) const
{
    assert(cd);
    assert(ff);
    const ExactIrrepCD<T>* icd=dynamic_cast<const ExactIrrepCD<T>*>(cd);
    assert(icd);
    assert(&*icd->itsBasisSet==static_cast<const IrrepBasisSet*>(this));
    const FittedFunctionImplementation<T>* ffi=dynamic_cast<const FittedFunctionImplementation<T>*>(ff);
    assert(ffi);
    double ret=0;
    typename Vector<T>::const_iterator c(ffi->itsFitCoeff.begin());
    const ERI3& repulsion=GetDataBase()->GetRepulsion3C(this,ffi->CastBasisSet());
    for(index_t i=0; c!=ffi->itsFitCoeff.end(); c++,i++)
        ret+=real((*c) * Dot(icd->itsDensityMatrix,repulsion[i]));
    return ret;
}

template <class T> double TBasisSetImplementation<T>::
GetCDOverlap  (const ChargeDensity* cd, const FittedFunction* ff) const
{
    assert(cd);
    assert(ff);
    const ExactIrrepCD<T>* icd=dynamic_cast<const ExactIrrepCD<T>*>(cd);
    assert(icd);
    assert(&*icd->itsBasisSet==static_cast<const IrrepBasisSet*>(this));
    const FittedFunctionImplementation<T>* ffi=dynamic_cast<const FittedFunctionImplementation<T>*>(ff);
    assert(ffi);
    double ret=0;
    typename Vector<T>::const_iterator c(ffi->itsFitCoeff.begin());
    const typename TIrrepBasisSet<T>::ERI3& overlap=GetDataBase()->GetOverlap3C(this,ffi->CastBasisSet());
    for(index_t i=0; c!=ffi->itsFitCoeff.end(); c++,i++)
        ret+=real((*c) * Dot(icd->itsDensityMatrix,overlap[i]));
    return ret;
}

//template <class T> double TBasisSetImplementation<T>::
//GetCDRepulsion(const ChargeDensity* cd) const
//{
//    assert(cd);
//    const ExactIrrepCD<T>* icd=dynamic_cast<const ExactIrrepCD<T>*>(cd);
//    assert(icd);
//    const TBasisSet<double>* tbs=dynamic_cast<const TBasisSet<double>*>(icd->itsBasisSet.get());
//    assert(tbs);
////    assert(&*icd->itsBasisSet==static_cast<const BasisSet*>(this));
//    SMatrix<T> Jab=GetRepulsion(icd->itsDensityMatrix,tbs);
//    return real(Dot(Jab,icd->itsDensityMatrix));
//}
//
//template <class T> double TBasisSetImplementation<T>::
//GetCDExchangeEnergy(const ChargeDensity* cd) const
//{
//    assert(cd);
//    const ExactIrrepCD<T>* icd=dynamic_cast<const ExactIrrepCD<T>*>(cd);
//    assert(icd);
//    assert(&*icd->itsBasisSet==static_cast<const BasisSet*>(this));
//    SMatrix<T> Kab=GetExchange(icd->itsDensityMatrix);
//    return real(Dot(Kab,icd->itsDensityMatrix));
//}

//
//  Load overlap (or repulsion) of this basis set with a scalar
//  funciton into a fitted function.
//
template <class T> void TBasisSetImplementation<T>::
SetFitOverlap  (FittedFunction* ff,const ScalarFunction<double>& sf) const
{
    assert(ff);
    FittedFunctionImplementation<T>* ffi=dynamic_cast<FittedFunctionImplementation<T>*>(ff);
    assert(ffi);
    assert(&*ffi->itsBasisSet==static_cast<const IrrepBasisSet*>(this));
    assert(!isnan(ffi->GetFitCoeff()));
    assert(!isnan(GetDataBase()->GetOverlap(sf)));
    ffi->GetFitCoeff()+=GetDataBase()->GetOverlap(sf);
}


template <class T> void TBasisSetImplementation<T>::
SetFitRepulsion(FittedFunction* ff,const ScalarFunction<double>& sf) const
{
    assert(ff);
    FittedFunctionImplementation<T>* ffi=dynamic_cast<FittedFunctionImplementation<T>*>(ff);
    assert(ffi);
    assert(&*ffi->itsBasisSet==static_cast<const IrrepBasisSet*>(this));
    assert(!isnan(ffi->GetFitCoeff()));
    assert(!isnan(GetDataBase()->GetRepulsion(sf)));
    ffi->GetFitCoeff()+=GetDataBase()->GetRepulsion(sf);
}



//-----------------------------------------------------------------------------
//
//  VectorFunction stuff.
//
template <class T> typename TBasisSetImplementation<T>::Vec TBasisSetImplementation<T>::
operator() (const RVec3& r) const
{
    Vec  ret(GetVectorSize());
    typename Vec::iterator i(ret.begin());
    for(auto b=this->beginT();b!=this->end();i++,b++) *i=(**b)(r);

    return ret;
}

template <class T> typename TBasisSetImplementation<T>::Vec3Vec TBasisSetImplementation<T>::
Gradient(const RVec3& r) const
{
    // No UT coverage
    Vec3Vec  ret(GetVectorSize());
    typename Vec3Vec::iterator i(ret.begin());
    for(auto b=this->beginT(); b!=this->end(); i++,b++) *i=b->Gradient(r);

    return ret;
}

template <class T> void TBasisSetImplementation<T>::Eval(const Mesh& mesh, Mat& mat) const
{
    // No UT coverage
    StreamableObject::SetToPretty();
    index_t i=1;
    for (auto b=this->beginT(); b!=this->end(); i++,b++)
    {
        mat.GetRow(i)=(**b)(mesh);
//		cout << "TBasisSetImplementation<T>::Eval (*b)(mesh)=" << (*b)(mesh) << std::endl;
//		cout << "TBasisSetImplementation<T>::Eval mat.GetRow(i)=" << mat.GetRow(i) << std::endl;
    }
}

template <class T> void TBasisSetImplementation<T>::EvalGrad(const Mesh& mesh, Vec3Mat& mat) const
{
    // No UT coverage
    index_t i=1;
    for (auto b=this->beginT(); b!=this->end(); i++,b++) mat.GetRow(i)=(**b).Gradient(mesh);
}

template <class T> LASolver<T>* TBasisSetImplementation<T>::GetLASolver() const
{
    if (!itsLASolver) 
    {
        SMat S=GetDataBase()->GetOverlap(this);
        itsLASolver=LASolver<T>::Factory(itsLAParams);
        itsLASolver->SetBasisOverlap(S);
    }
    assert(itsLASolver);
    return itsLASolver;
}

//-----------------------------------------------------------------------------
//
//  Streamable stuff.
//
template <class T> std::ostream& TBasisSetImplementation<T>::Write(std::ostream& os) const
{
    if(!StreamableObject::Pretty())
    {
        VectorFunctionBuffer<T>::Write(os);
        os << *itsAnalyticIE << itsDataBase;
    }
    return os;
}

template <class T> std::istream& TBasisSetImplementation<T>::Read(std::istream& is)
{
    VectorFunctionBuffer<T>::Read(is);
    AnalyticIE<T>* ie=AnalyticIE<T>::Factory(is);
    is >> *ie;
    is >> itsDataBase;
    Insert(ie); //Fix up lots of pointers.

    return is;
};

template class TBasisSetImplementation<double>;
//template class TBasisSetImplementation<std::complex<double> >;

