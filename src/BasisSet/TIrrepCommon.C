// File: TBasisSetImplementation.C


#include "Imp/BasisSet/TIrrepCommon.H"
#include <BasisSet.H>
#include <QuantumNumber.H>
#include <AnalyticIE.H>
#include <IntegralDataBase.H>
#include <LASolver/LASolver.H>
#include <Hamiltonian.H>

#include "Imp/Containers/ERI4.H"
#include "OrbitalImplementation/TOrbitalGroupImplementation.H"
#include <cassert>
#include <iostream>

//-----------------------------------------------------------------------------
//
//  Construction zone
//
template <class T> TIrrepBasisSetCommon<T>::TIrrepBasisSetCommon()
    : itsDataBase      ( )
    , itsLASolver      (0)
{};

template <class T> TIrrepBasisSetCommon<T>::TIrrepBasisSetCommon(const LinearAlgebraParams& lap,IntegralDataBase<T>* theDataBase)
    : itsLAParams      (lap)
    , itsDataBase      (theDataBase)
    , itsLASolver      (0)
{};

template <class T> TIrrepBasisSetCommon<T>::TIrrepBasisSetCommon(const TIrrepBasisSetCommon<T>& bs)
    : itsLAParams      (bs.itsLAParams)
    , itsDataBase      (bs.itsDataBase)
    , itsLASolver      (0)
{};

template <class T> TIrrepBasisSetCommon<T>::~TIrrepBasisSetCommon()
{
    delete itsLASolver;
}

//-----------------------------------------------------------------------------
//
//  Post construction initializations called by dervied classes.
//
template <class T> void TIrrepBasisSetCommon<T>::Insert(AnalyticIE<T>* ie)
{
    assert(ie);
    itsDataBase->Insert(ie);
}

template <class T> IntegralDataBase<T>* TIrrepBasisSetCommon<T>::GetDataBase() const
{
    assert(&*itsDataBase);
    return itsDataBase;
}

 

//-----------------------------------------------------------------------------
//
//  Build up the hamiltonian matrix, and generate new orbitals from the eigen
//  vectors.  The hamiltonian controls which terms are included in H.
//
// TODO: Why do we need to pass in const rc_ptr<const BasisSet>& rc ????
//
template <class T> OrbitalGroup* TIrrepBasisSetCommon<T>::
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

template <class T> typename TIrrepBasisSetCommon<T>::RVec TIrrepBasisSetCommon<T>::
GetCharge() const
{
    return GetDataBase()->GetCharge(this);
}

template <class T> typename TIrrepBasisSetCommon<T>::SMat TIrrepBasisSetCommon<T>::
GetOverlap() const
{
    return GetDataBase()->GetOverlap(this);
}

template <class T> typename TIrrepBasisSetCommon<T>::SMat TIrrepBasisSetCommon<T>::
GetInverseRepulsion() const
{
    return GetDataBase()->GetInverseRepulsion(this);
}

template <class T> typename TIrrepBasisSetCommon<T>::SMat TIrrepBasisSetCommon<T>::
GetInverseOverlap() const
{
    return GetDataBase()->GetInverseOverlap(this);
}

template <class T> IrrepBasisSet::SMat TIrrepBasisSetCommon<T>::
GetKinetic() const
{
    return GetDataBase()->GetKinetic(this);
}
template <class T> IrrepBasisSet::SMat TIrrepBasisSetCommon<T>::
GetNuclear(const Cluster* cl) const
{
    return itsDataBase->GetNuclear(this,*cl);
}

 template <class T> typename TIrrepBasisSetCommon<T>::Mat TIrrepBasisSetCommon<T>::
 GetRepulsion(const IrrepBasisSet* ff) const
 {
    const TIrrepBasisSet<T>* tff=dynamic_cast<const TIrrepBasisSet<T>*>(ff);
    assert(tff);
    return GetDataBase()->GetRepulsion(this,tff);
 }
 
 template <class T> typename TIrrepBasisSetCommon<T>::Mat TIrrepBasisSetCommon<T>::
 GetOverlap(const Mesh* m,const IrrepBasisSet* ff) const
 {
    const TIrrepBasisSet<T>* tff=dynamic_cast<const TIrrepBasisSet<T>*>(ff);
    assert(tff);
    return GetDataBase()->GetOverlap(m,*this,*tff);
 }

 template <class T> typename TIrrepBasisSetCommon<T>::RVec TIrrepBasisSetCommon<T>::
 GetOverlap(const Mesh* m,const ScalarFunction<double>* sf) const
 {
    return GetDataBase()->GetOverlap(m,*this,*sf);
 }
 template <class T> typename TIrrepBasisSetCommon<T>::RVec TIrrepBasisSetCommon<T>::
 GetRepulsion(const Mesh* m,const ScalarFunction<double>* sf) const
 {
    return GetDataBase()->GetRepulsion(m,*this,*sf);
 }
 
template <class T> IrrepBasisSet::SMat TIrrepBasisSetCommon<T>::
GetRepulsion(const SMat& Dcd, const TIrrepBasisSet<T>* bs_cd) const
{
    assert(!isnan(Dcd));
    ERI4view J=GetDataBase()->GetRepulsion4C(*this,*bs_cd);
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

template <class T> IrrepBasisSet::SMat TIrrepBasisSetCommon<T>::
GetExchange(const SMat& Dcd, const TIrrepBasisSet<T>* bs_cd) const
{
    assert(!isnan(Dcd));
    ERI4view K=GetDataBase()->GetExchange4C(*this,*bs_cd);
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
template <class T> typename TIrrepBasisSetCommon<T>::RVec TIrrepBasisSetCommon<T>::
GetOverlap3C(const SMat& Dcd, const IrrepBasisSet* ff) const
{
    RVec ret(ff->size());
    const ERI3& S=GetDataBase()->GetOverlap3C(this,ff);
    for(auto i:ret.indices())
        ret(i)=Dot(Dcd,S[i-1]);
    return ret;
}

template <class T> typename TIrrepBasisSetCommon<T>::RVec TIrrepBasisSetCommon<T>::
GetRepulsion3C(const SMat& Dcd, const IrrepBasisSet* ff) const
{
    RVec ret(ff->size());
    const ERI3& repulsion=GetDataBase()->GetRepulsion3C(this,ff);
    for(auto i:ret.indices())
        ret(i)=Dot(Dcd,repulsion[i-1]);
    return ret;
}

template <class T> const typename TIrrepBasisSetCommon<T>::ERI3& TIrrepBasisSetCommon<T>::
GetOverlap3C(const IrrepBasisSet* ff) const
{
    return GetDataBase()->GetOverlap3C(this,ff);
}

template <class T> const typename TIrrepBasisSetCommon<T>::ERI3& TIrrepBasisSetCommon<T>::
GetRepulsion3C(const IrrepBasisSet* ff) const
{
    return GetDataBase()->GetRepulsion3C(this,ff);
}

//-----------------------------------------------------------------------------
//
//  VectorFunction stuff.
//
template <class T> typename TIrrepBasisSetCommon<T>::Vec TIrrepBasisSetCommon<T>::
operator() (const RVec3& r) const
{
    Vec  ret(GetVectorSize());
    typename Vec::iterator i(ret.begin());
    for(auto b=this->beginT();b!=this->end();i++,b++) *i=(**b)(r);

    return ret;
}

template <class T> typename TIrrepBasisSetCommon<T>::Vec3Vec TIrrepBasisSetCommon<T>::
Gradient(const RVec3& r) const
{
    // No UT coverage
    Vec3Vec  ret(GetVectorSize());
    typename Vec3Vec::iterator i(ret.begin());
    for(auto b=this->beginT(); b!=this->end(); i++,b++) *i=b->Gradient(r);

    return ret;
}

template <class T> void TIrrepBasisSetCommon<T>::Eval(const Mesh& mesh, Mat& mat) const
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

template <class T> void TIrrepBasisSetCommon<T>::EvalGrad(const Mesh& mesh, Vec3Mat& mat) const
{
    // No UT coverage
    index_t i=1;
    for (auto b=this->beginT(); b!=this->end(); i++,b++) mat.GetRow(i)=(**b).Gradient(mesh);
}

template <class T> LASolver<T>* TIrrepBasisSetCommon<T>::GetLASolver() const
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
template <class T> std::ostream& TIrrepBasisSetCommon<T>::Write(std::ostream& os) const
{
    if(!StreamableObject::Pretty())
    {
        os <<  itsDataBase;
    }
    return os;
}

template <class T> std::istream& TIrrepBasisSetCommon<T>::Read(std::istream& is)
{
    is >> itsDataBase;

    return is;
};

template class TIrrepBasisSetCommon<double>;
//template class TBasisSetImplementation<std::complex<double> >;

