// File: dbheap.cpp  Implement a heap storage integral data base.



#include "DFTDataBase/HeapDB/HeapDB.H"
#include "Mesh/MeshIntegrator.H"
#include "BasisSet.H"
#include "Cluster.H"
#include "AnalyticIE.H"
#include "oml/vector.h"
#include "Imp/Containers/stl_io.h"
#include <iostream>
//#include <cassert>
#include <algorithm> //find

//------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> HeapDB<T>::HeapDB()
    :itsAnalyticIE   (0)
{}

template <class T> HeapDB<T>::HeapDB(AnalyticIE<T>* ie, const IEClient* iec)
    :itsAnalyticIE(ie)   
    ,istIEClient  (iec)
{
    assert(itsAnalyticIE);
    assert(istIEClient);
}



template <class T> void HeapDB<T>::WipeCleanAllData()
{
}


template <class T> void HeapDB<T>::Insert(AnalyticIE<T>* ie)
{
    assert(ie);
    if (!itsAnalyticIE) 
        itsAnalyticIE=ie;
    else
        delete ie;
}


template <class T> bool HeapDB<T>::operator==(const IntegralDataBase<T>& idb) const
{
    const HeapDB* hdb=dynamic_cast<const HeapDB*>(&idb);
    return hdb;
    
}
//-------------------------------------------------------------------------
//
//  Assorted fluff.
//

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
template <class T> const typename HeapDB<T>::RVec& HeapDB<T>::GetNumericalNormalization(const Mesh* m,bs_t& a)
{
    MeshIntegrator<T> mintegrator(m);
    id2c_t key=std::make_tuple(qchem::Normalization,a.GetID());
    if (auto i = its1C.find(key); i==its1C.end())
        return its1C[key] =mintegrator.Normalize(a);
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


//
//  THis get called in code by FittedFunctionImplementation<T>::FitGet2CenterOverlap(const IrrepBasisSet* bs) const
//  But does not get used at run time.  The fit uses the repulsion version instead.
//
template <class T> const typename HeapDB<T>::Mat& HeapDB<T>::GetOverlap(const Mesh* m,bs_t& a,bs_t& b)
{
    // No UT coverage.
    assert(false);
    return *new Mat();
//    assert(itsNumericalIE);
//    id2cx_t key=std::make_tuple(qchem::Overlap2C,a.GetID(),b.GetID());
//    if (auto i = its2Cx.find(key); i==its2Cx.end())
//        return its2Cx[key] = itsNumericalIE->MakeOverlap(a,b);
//    else
//        return i->second;
}


//
//  DO not try and cache these because the ScalarFunction f changes with iterations.
//
template <class T> const typename HeapDB<T>::Vec HeapDB<T>::GetOverlap(const Mesh* m,bs_t& bs,Rf& f)
{
    const RVec& n=GetNumericalNormalization(m,bs);
    MeshIntegrator<T> mintegrator(m);
    return DirectMultiply(mintegrator.Overlap(f,bs),n);
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

template <class T> const typename HeapDB<T>::Vec HeapDB<T>::GetRepulsion(const Mesh* m,bs_t& bs,Rf& f)
{
    const RVec& n=GetNumericalNormalization(m,bs);
    MeshIntegrator<T> mintegrator(m);
    return DirectMultiply(mintegrator.Repulsion(f,bs),n);
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


//
//  THis get called in code by FittedFunctionImplementation<T>::FitGet2CenterOverlap(const IrrepBasisSet* bs) const
//  But does not get used at run time.  The fit uses the repulsion version instead.
//
template <class T> const typename HeapDB<T>::RVec& HeapDB<T>::GetCharge(const Mesh*,bs_t&  a)
{
    assert(false);
    return *new RVec();
//    assert(itsNumericalIE);
//    id2c_t key=std::make_tuple(qchem::Charge,a.GetID());
//    if (auto i = its1C.find(key); i==its1C.end())
//        return its1C[key] =itsNumericalIE->MakeCharge(a) ;
//    else
//        return i->second;
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
    assert(istIEClient);
    assert(itsJTable.GetSize()==0);
    assert(itsKTable.GetSize()==0); //They should be synchronized.
    itsAnalyticIE->Make4C(itsJTable,itsKTable,istIEClient);
}

template <class T> ERI4view  HeapDB<T>::GetRepulsion4C(bs_t& a,bs_t& b)
{
   if (itsJTable.GetSize()==0) BuildERIs(); 
   return ERI4view(itsJTable,a.GetStartIndex(),b.GetStartIndex());
}

template <class T> ERI4view  HeapDB<T>::GetExchange4C (bs_t& a,bs_t& b)
{
   if (itsJTable.GetSize()==0) BuildERIs(); 
   if (itsKTable.GetSize()==0)
        return ERI4view(itsJTable,a.GetStartIndex(),b.GetStartIndex());
    else
        return ERI4view(itsKTable,a.GetStartIndex(),b.GetStartIndex());
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
    assert(itsAnalyticIE);
    id2c_t key=std::make_tuple(qchem::InvOverlap,a->GetID());
    if (auto i = its2C.find(key); i==its2C.end())
        return its2C[key] =itsAnalyticIE->MakeInverse(GetOverlap(a));
    else
        return i->second;
}

template <class T> const typename HeapDB<T>::SMat& HeapDB<T>::GetInverseRepulsion(iec_t* a)
{
//     std::cout << GetOverlap(a) << std::endl;
    assert(itsAnalyticIE);
    id2c_t key=std::make_tuple(qchem::InvRepulsion,a->GetID());
    if (auto i = its2C.find(key); i==its2C.end())
        return its2C[key] =itsAnalyticIE->MakeInverse(GetRepulsion(a));
    else
        return i->second;
}

#ifndef UT_COVERAGE_ONLY

template <class T> const typename HeapDB<T>::Mat& HeapDB<T>::GetRepulsion(bs_t& a,bs_t& b)
{
    assert(itsNumericalIE);
    id2cx_t key=std::make_tuple(qchem::Repulsion2C,a.GetID(),b.GetID());
    if (auto i = its2Cx.find(key); i==its2Cx.end())
        return its2Cx[key] = itsNumericalIE->MakeRepulsion(a,b);
    else
        return i->second;
}

template <class T>  const typename HeapDB<T>::SMat& HeapDB<T>::GetOverlap(bs_t& a)
{
    assert(itsNumericalIE);
    id2c_t key=std::make_tuple(qchem::Overlap2C,a.GetID());
    if (auto i = its2C.find(key); i==its2C.end())
        return its2C[key] =itsNumericalIE->MakeOverlap(a);
    else
        return i->second;
}
#endif // UT_COVERAGE_ONLY



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
