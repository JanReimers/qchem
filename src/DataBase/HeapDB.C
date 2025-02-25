// File: HeapDB.H  Implement a heap storage integral data base.

#include "Imp/DataBase/HeapDB.H"
#include "Imp/Integrals/MeshIntegrator.H"
#include <BasisSet.H>
#include <Cluster.H>
#include <AnalyticIE.H>
#include "oml/vector.h"
#include <iomanip>

//------------------------------------------------------------------------
//
//  Construction zone.
//
template <class T> HeapDB<T>::HeapDB()
    :itsAnalyticIE   (0)
{}

template <class T> HeapDB<T>::HeapDB(AnalyticIE<T>* ie)
    :itsAnalyticIE(ie)   
{
    assert(itsAnalyticIE);
}

template <class T> HeapDB<T>::~HeapDB()
{
    Report(std::cout);
    delete itsAnalyticIE;
}

template <class T> size_t Size(const Vector <T>& m) {return m.size();}
template <class T> size_t Size(const Matrix <T>& m) {return m.size();}
template <class T> size_t Size(const SMatrix<T>& m) {return m.size();}
                   size_t Size(const ERI4      & m) {return m.size();}
template <class M> size_t Size(const std::vector<M> & v) 
{
    size_t N=0;
    for (auto i:v) 
        N+=Size(i);
    return N;
}

template <class K, class M> size_t Size(const std::map<K,M>& m)
{
    size_t N=0;
    for (auto i:m) 
        N+=Size(i.second);
    return N;
}

using std::setw;

template <class T> void HeapDB<T>::Report(std::ostream& os) const
{
    size_t N1=Size(its1C)+Size(its1Cx);
    size_t N2=Size(its2C)+Size(its2Cx)+Size(its2CNuc);
    size_t N3=Size(its3C);
    size_t N4=Size(Jac)+Size(Kab);
    
    os << "Heap DB storage report:" << std::endl;
    os << "    " << setw(10) << N1 << " 1 centre integrals." << std::endl;
    os << "    " << setw(10) << N2 << " 2 centre integrals." << std::endl;
    os << "    " << setw(10) << N3 << " 3 centre integrals." << std::endl;
    os << "    " << setw(10) << N4 << " 4 centre integrals." << std::endl;
    
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
template <class T> const typename HeapDB<T>::SMat& HeapDB<T>::GetRestMass(iec_t* a)
{
    assert(itsAnalyticIE);
    id2c_t key=std::make_tuple(qchem::RestMass,a->GetID());
    if (auto i = its2C.find(key); i==its2C.end())
        return its2C[key] =itsAnalyticIE->MakeRestMass(a);
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


using std::cout;
using std::endl;
template <class T> ERI4 HeapDB<T>::GetDirect__4C(bs_t& a,bs_t& c)
{
    assert(a.GetID()<=c.GetID());
    if (Jac.size()==0) 
        itsAnalyticIE->MakeDirect  (Jac);
    //cout << "GetRepulsion4C_new a,c=" << a.GetIndex() << " " << c.GetIndex() << endl;
    assert(Jac.find(a.GetID())!=Jac.end());
    assert(Jac[a.GetID()].find(c.GetID())!=Jac[a.GetID()].end());
    
    return Jac[a.GetID()][c.GetID()];
}
template <class T> ERI4 HeapDB<T>::GetExchange4C(bs_t& a,bs_t& b)
{
    assert(a.GetID()<=b.GetID());
    if (Kab.size()==0)
        itsAnalyticIE->MakeExchange(Kab); 
    //cout << "GetExchange4C_new a,b=" << a.GetIndex() << " " << b.GetIndex() << endl;
    assert(Kab.find(a.GetID())!=Kab.end());
    assert(Kab[a.GetID()].find(b.GetID())!=Kab[a.GetID()].end());
    
    return Kab[a.GetID()][b.GetID()];
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



template <class T> const typename HeapDB<T>::SMat& HeapDB<T>::GetInverseOverlap(iec_t* a, const LAParams& lap)
{
     //Only used for numerical IE.
    assert(itsAnalyticIE);
    id2c_t key=std::make_tuple(qchem::InvOverlap,a->GetID());
    if (auto i = its2C.find(key); i==its2C.end())
        return its2C[key] =itsAnalyticIE->MakeInverse(GetOverlap(a),lap);
    else
        return i->second;
}

template <class T> const typename HeapDB<T>::SMat& HeapDB<T>::GetInverseRepulsion(iec_t* a,const LAParams& lap)
{
//     std::cout << GetOverlap(a) << std::endl;
    assert(itsAnalyticIE);
    id2c_t key=std::make_tuple(qchem::InvRepulsion,a->GetID());
    if (auto i = its2C.find(key); i==its2C.end())
        return its2C[key] =itsAnalyticIE->MakeInverse(GetRepulsion(a),lap);
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
GetInverseOverlap(iec_t* a,const LAParams&)
{
    std::cerr << "Sorry, inverse of complex matrix is not implemented" << std::endl;
    exit(-1);
    return m;
}

template <> const HeapDB<std::complex<double> >::SMat& HeapDB<std::complex<double> >::
GetInverseRepulsion(iec_t* ,const LAParams&)
{
    std::cerr << "Sorry, inverse of complex matrix is not implemented" << std::endl;
    exit(-1);
    return m;
}

template class HeapDB<double>;
template class HeapDB<std::complex<double> >;
