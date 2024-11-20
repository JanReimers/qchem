// File: TBasisSetImplementation.C


#include "Imp/BasisSet/TIrrepCommon.H"
#include "Imp/Containers/ERI4.H"
#include <BasisSet.H>
#include <QuantumNumber.H>
#include <AnalyticIE.H>
#include <IntegralDataBase.H>
#include <LASolver.H>
#include <Hamiltonian.H>
#include <cassert>
#include <iostream>

//-----------------------------------------------------------------------------
//
//  Construction zone
//
template <class T> TIrrepBasisSetCommon<T>::TIrrepBasisSetCommon()
    : itsDataBase      ( )
{};

template <class T> TIrrepBasisSetCommon<T>::TIrrepBasisSetCommon(const LAParams& lap,IntegralDataBase<T>* theDataBase)
    : itsLAParams      (lap)
    , itsDataBase      (theDataBase)
{};

template <class T> TIrrepBasisSetCommon<T>::TIrrepBasisSetCommon(const TIrrepBasisSetCommon<T>& bs)
    : itsLAParams      (bs.itsLAParams)
    , itsDataBase      (bs.itsDataBase)
{};

template <class T> TIrrepBasisSetCommon<T>::~TIrrepBasisSetCommon()
{
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

 

template <class T>  LASolver<double>* TIrrepBasisSetCommon<T>::CreateSolver() const
{
    LASolver<double>* las=LASolver<double>::Factory(itsLAParams);
    las->SetBasisOverlap(GetDataBase()->GetOverlap(this));
    return las;
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
GetInverseRepulsion(const LAParams& lap) const
{
    return GetDataBase()->GetInverseRepulsion(this,lap);
}

template <class T> typename TIrrepBasisSetCommon<T>::SMat TIrrepBasisSetCommon<T>::
GetInverseOverlap(const LAParams& lap) const
{
    return GetDataBase()->GetInverseOverlap(this,lap);
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

 using std::cout;
 using std::endl;
template <class T> IrrepBasisSet::SMat TIrrepBasisSetCommon<T>::
GetRepulsion(const SMat& Dcd, const TIrrepBasisSet<T>* cd) const
{
    assert(!isnan(Dcd));
    assert(Max(fabs(Dcd))>0.0);  //Don't waste time!
    const TIrrepBasisSet<T>* ab=this;
    if (ab->GetID()<=cd->GetID())
    {
        ERIJ1 Jabcd=GetDataBase()->GetRepulsion4C_new(*ab,*cd);
        return Jabcd*Dcd;
    }
    else
    {
        ERIJ1 Jcdab=GetDataBase()->GetRepulsion4C_new(*cd,*ab);
        return Dcd*Jcdab;        
    }
}

#include <iomanip>
template <class T> IrrepBasisSet::SMat TIrrepBasisSetCommon<T>::
GetExchange(const SMat& Dcd, const TIrrepBasisSet<T>* bs_cd) const
{
    assert(!isnan(Dcd));
    assert(Max(fabs(Dcd))>0.0);  //Don't waste time!
    ERIK  K1=GetDataBase()->GetExchange4C_new(*this,*bs_cd);
    int Nab=this->GetNumFunctions();
    int Ncd=bs_cd->GetNumFunctions();

    SMat Kab(Nab,Nab);
    for (int ia=1; ia<=Nab; ia++)
        for (int ib=ia; ib<=Nab; ib++)
        {
            T Kab_temp=0;
            for (int ic=1; ic<=Ncd; ic++)
            {
                for (int id=1; id<=Ncd; id++) //Possible symmetric optimization here.  Need to be careful to handle complex D and ERIs.
                {
                    assert(K1(ia,ib,ic,id)==K1(ia,ib,id,ic));
                    Kab_temp+= K1(ia,ib,ic,id)*Dcd(ic,id);
                }
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
    Vec  ret(this->size());
    typename Vec::iterator i(ret.begin());
    for(auto b=this->beginT();b!=this->end();i++,b++) *i=(**b)(r);

    return ret;
}

template <class T> typename TIrrepBasisSetCommon<T>::Vec3Vec TIrrepBasisSetCommon<T>::
Gradient(const RVec3& r) const
{
    // No UT coverage
    Vec3Vec  ret(this->size());
    typename Vec3Vec::iterator i(ret.begin());
    for(auto b=this->beginT(); b!=this->end(); i++,b++) *i=b->Gradient(r);

    return ret;
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

