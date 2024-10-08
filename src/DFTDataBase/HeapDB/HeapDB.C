// File: dbheap.cpp  Implement a heap storage integral data base.



#include "DFTDataBase/HeapDB/HeapDB.H"
#include "BasisSet.H"
#include "Cluster.H"
#include "NumericalIE.H"
#include "AnalyticIE.H"
#include "oml/vector.h"
#include "Misc/stl_io.h"
#include <iostream>
//#include <cassert>
#include <algorithm> //find

//------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> HeapDB<T>::HeapDB()
    :itsNumericalIE    (0)
    ,itsAnalyticIE   (0)
    ,itsBasisSet          (0)
    
{}



template <class T> void HeapDB<T>::WipeCleanAllData()
{
}

template <class T> void HeapDB<T>::Insert(bs_t* bs,const NumericalIE<T>* ie)
{
    assert(bs);
    assert(ie);
    itsBasisSet=bs;
    itsNumericalIE=ie;
}

template <class T> void HeapDB<T>::Insert(bs_t* bs,const AnalyticIE<T>* ie)
{
    assert(bs);
    assert(ie);
    itsBasisSet=bs;
    itsAnalyticIE=ie;
}

template <class T> void HeapDB<T>::Insert(const BasisGroup* bg)
{
    assert(bg);
    itsBasisGroup=bg;
}


template <class T> void HeapDB<T>::Insert(const ERI4& J, const ERI4& K)
{
    itsJTable=J;
    itsKTable=K;
}

template <class T> bool HeapDB<T>::operator==(const IntegralDataBase<T>& idb) const
{
    if (!itsBasisSet || itsBasisSet->GetID()==0) return false;
    const HeapDB* hdb=dynamic_cast<const HeapDB*>(&idb);
    if (!hdb) return false;
    if (!hdb->itsBasisSet || hdb->itsBasisSet->GetID()==0) return false;
    return itsBasisSet->GetID()==hdb->itsBasisSet->GetID();
}
//-------------------------------------------------------------------------
//
//  Assorted fluff.
//
template <class T> const NumericalIE<T>* HeapDB<T>::GetIntegralEngine() const
{
    assert(itsNumericalIE);
    return itsNumericalIE;
}

template <class T> IntegralDataBase<T>* HeapDB<T>::Clone() const
{
    return new HeapDB(*this);
}


//------------------------------------------------------------------
//
//  Pickle and un pickle the integral engine, all matricies, and
//  all lists of ID's.
//
template <class T> std::ostream& HeapDB<T>::Write(std::ostream& os) const
{
//
//  Write the IntegralEngine, which better be initialized.
//
    if(!Binary()) os << std::endl;
    UniqueID::Write(os);


    if (Binary())
    {
        
    }
    else
    {
        
    }
    return os;
}

template <class T> std::istream& HeapDB<T>::Read (std::istream& is)
{
    UniqueID::Read(is);

    if (Binary())
    {
       
    }
    else
    {
       
    }
    return is;
}





//---------------------------------------------------------------------------------
//
//  If possible return cached integral tables.  Calculate using integral engine only if required.
//  In some cases the numerical IE takes priority.
//
template <class T>  const typename HeapDB<T>::SMat& HeapDB<T>::GetOverlap(bs_t& a)
{
    assert(itsNumericalIE);
    id2c_t key=std::make_tuple(qchem::Overlap2C,a.GetID());
    if (auto i = its2C.find(key); i==its2C.end())
        return its2C[key] =itsNumericalIE->MakeOverlap(a);
    else
        return i->second;
}

template <class T>  const typename HeapDB<T>::SMat& HeapDB<T>::GetOverlap(iec_t* a)
{
    assert(itsAnalyticIE);
    id2c_t key=std::make_tuple(qchem::Overlap2C,a->GetID());
    if (auto i = its2C.find(key); i==its2C.end())
        return its2C[key] =itsAnalyticIE->MakeOverlap(a);
    else
        return i->second;
}

template <class T> const typename HeapDB<T>::Mat& HeapDB<T>::GetOverlap(bs_t& a,bs_t& b)
{
    // No UT coverage.
    assert(false);
    assert(itsNumericalIE);
    id2cx_t key=std::make_tuple(qchem::Overlap2C,a.GetID(),b.GetID());
    if (auto i = its2Cx.find(key); i==its2Cx.end())
        return its2Cx[key] = itsNumericalIE->MakeOverlap(a,b);
    else
        return i->second;
}


//
//  DO not try and cache these because the ScalarFunction f changes with iterations.
//
template <class T> const typename HeapDB<T>::Vec HeapDB<T>::GetOverlap(const ScalarFunction<double>& f)
{
    assert(itsNumericalIE);    
    return itsNumericalIE->MakeOverlap(f);;
}



template <class T> const typename HeapDB<T>::SMat& HeapDB<T>::GetRepulsion(iec_t* a )
{
    assert(itsAnalyticIE);
    id2c_t key=std::make_tuple(qchem::Repulsion2C,a->GetID());
    if (auto i = its2C.find(key); i==its2C.end())
        return its2C[key] =itsAnalyticIE->MakeRepulsion(a);
    else
        return i->second;
}

template <class T> const typename HeapDB<T>::Vec HeapDB<T>::GetRepulsion(const ScalarFunction<double>& f)
{
    assert(itsNumericalIE);
    return itsNumericalIE->MakeRepulsion(f);
}

template <class T> const typename HeapDB<T>::Mat& HeapDB<T>::GetRepulsion(bs_t& a,bs_t& b)
{
    assert(itsAnalyticIE || itsNumericalIE);
    id2cx_t key=std::make_tuple(qchem::Repulsion2C,a.GetID(),b.GetID());
    if (auto i = its2Cx.find(key); i==its2Cx.end())
        return its2Cx[key] = (itsNumericalIE 
            ? itsNumericalIE->MakeRepulsion(a,b) 
            :  itsAnalyticIE->MakeRepulsion(a.GetAnalyticIE(),b.GetAnalyticIE()));
    else
        return i->second;
}

template <class T> const typename HeapDB<T>::Mat& HeapDB<T>::GetRepulsion(iec_t* a,iec_t* b)
{
        assert(itsAnalyticIE);
    id2cx_t key=std::make_tuple(qchem::Repulsion2C,a->GetID(),b->GetID());
    if (auto i = its2Cx.find(key); i==its2Cx.end())
        return its2Cx[key] = itsAnalyticIE->MakeRepulsion(a,b);
    else
        return i->second;
}

template <class T> const typename HeapDB<T>::SMat& HeapDB<T>::GetKinetic(iec_t* a)
{
    assert(itsAnalyticIE);
    id2c_t key=std::make_tuple(qchem::Kinetic,a->GetID());
    if (auto i = its2C.find(key); i==its2C.end())
        return its2C[key] =itsAnalyticIE->MakeKinetic(a);
    else
        return i->second;
}

template <class T> const typename HeapDB<T>::RVec& HeapDB<T>::GetNormalization(bs_t& a)
{
    assert(itsNumericalIE);
    id2c_t key=std::make_tuple(qchem::Normalization,a.GetID());
    if (auto i = its1C.find(key); i==its1C.end())
        return its1C[key] =itsNumericalIE->MakeNormalization(a);
    else
        return i->second;
}

template <class T> const typename HeapDB<T>::RVec& HeapDB<T>::GetCharge(bs_t&  a)
{
    assert(itsNumericalIE || itsAnalyticIE);
    id2c_t key=std::make_tuple(qchem::Charge,a.GetID());
    if (auto i = its1C.find(key); i==its1C.end())
        return its1C[key] =itsNumericalIE->MakeCharge(a) ;
    else
        return i->second;
}

template <class T> const typename HeapDB<T>::RVec& HeapDB<T>::GetCharge(iec_t* a)
{
    assert(itsAnalyticIE);
    id2c_t key=std::make_tuple(qchem::Charge,a->GetID());
    if (auto i = its1C.find(key); i==its1C.end())
        return its1C[key] =itsAnalyticIE->MakeCharge(a);
    else
        return i->second;
}

template <class T> const typename HeapDB<T>::ERI3& HeapDB<T>::GetOverlap3C(iec_t* ab,iec_t* c )
{
    assert(itsAnalyticIE);
    id3c_t key=std::make_tuple(qchem::Overlap3C,ab->GetID(),c->GetID());
    if (auto i = its3C.find(key); i==its3C.end())
        return its3C[key] = itsAnalyticIE->MakeOverlap3C(ab,c);
    else
        return i->second;
}

template <class T> const typename HeapDB<T>::ERI3& HeapDB<T>::GetRepulsion3C(iec_t* ab,iec_t* c)
{
    assert(itsAnalyticIE);
    id3c_t key=std::make_tuple(qchem::Repulsion3C,ab->GetID(),c->GetID());
    if (auto i = its3C.find(key); i==its3C.end())
        return its3C[key] = itsAnalyticIE->MakeRepulsion3C(ab,c);
    else
        return i->second;
}

template <class T> void HeapDB<T>::BuildERIs()
{
    assert(itsAnalyticIE);
    assert(itsJTable.GetSize()==0);
    assert(itsKTable.GetSize()==0); //They should be synchronized.
    auto [J,K]=itsAnalyticIE->Make4C(itsBasisGroup->Flatten());
    itsBasisGroup->Insert(J,K); //Eeach basis set has its own HeapDB.
}

template <class T> ERI4view HeapDB<T>::GetRepulsion4C(bs_t* bs_cd)
{
    assert(itsBasisSet);
    assert(itsAnalyticIE);
   if (itsJTable.GetSize()==0) BuildERIs(); 
   assert(bs_cd);
   int ss=bs_cd->GetStartIndex();
   return ERI4view(itsJTable,itsBasisSet->GetStartIndex(),ss);
}

template <class T> ERI4view HeapDB<T>::GetExchange4C (bs_t* bs_cd)
{
   if (itsJTable.GetSize()==0) BuildERIs(); 
   if (itsKTable.GetSize()==0)
        return ERI4view(itsJTable,itsBasisSet->GetStartIndex(),bs_cd->GetStartIndex());
    else
        return ERI4view(itsKTable,itsBasisSet->GetStartIndex(),bs_cd->GetStartIndex());
}
//-------------------------------------------------------------------------
//
//  These guys have to check the ID lists to see if it integrals are
//  already known.
//

template <class T> const typename HeapDB<T>::SMat& HeapDB<T>::GetNuclear(iec_t* a,const Cluster& cl)
{
    assert(itsAnalyticIE);
    id2cx_t key=std::make_tuple(qchem::Nuclear,a->GetID(),cl.GetID());
    if (auto i = its2CNuc.find(key); i==its2CNuc.end())
        return its2CNuc[key] =itsAnalyticIE->MakeNuclear(a,cl);
    else
        return i->second;
}   



template <class T> const typename HeapDB<T>::SMat& HeapDB<T>::GetInverseOverlap(iec_t* a)
{
     //Only used for numerical IE.
    assert(itsNumericalIE);
    id2c_t key=std::make_tuple(qchem::InvOverlap,a->GetID());
    if (auto i = its2C.find(key); i==its2C.end())
        return its2C[key] =itsNumericalIE->MakeInverse(GetOverlap(a));
    else
        return i->second;
}

template <class T> const typename HeapDB<T>::SMat& HeapDB<T>::GetInverseRepulsion(iec_t* a)
{
    assert(itsAnalyticIE);
    assert(!itsNumericalIE); //Do we need to support this?
    id2c_t key=std::make_tuple(qchem::InvRepulsion,a->GetID());
    if (auto i = its2C.find(key); i==its2C.end())
        return its2C[key] =itsAnalyticIE->MakeInverse(GetRepulsion(a));
    else
        return i->second;
}

SMatrix<std::complex<double> > m; //dummy just to get the compiler to shut up.

template <> const HeapDB<std::complex<double> >::SMat& HeapDB<std::complex<double> >::
GetInverseOverlap(iec_t* a)
{
    std::cerr << "Sorry, inverse of complex matrix is not implemented" << std::endl;
    exit(-1);
    return m;
}

template <> const HeapDB<std::complex<double> >::SMat& HeapDB<std::complex<double> >::
GetInverseRepulsion(iec_t* a)
{
    std::cerr << "Sorry, inverse of complex matrix is not implemented" << std::endl;
    exit(-1);
    return m;
}

template class HeapDB<double>;
template class HeapDB<std::complex<double> >;
