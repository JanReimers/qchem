// File: HeapDB.H  Implement a heap storage integral data base.

#include "Imp/DataBase/HeapDB.H"
#include "Imp/Integrals/MeshIntegrator.H"
#include <Irrep_BS.H>
#include <Cluster.H>
#include <AnalyticIE.H>
#include "oml/vector.h"
#include <iomanip>

//------------------------------------------------------------------------
//
//  Construction zone.
//
// template <class T> HeapDB<T>::HeapDB()
//     :itsAnalyticIE   (0)
// {}

// template <class T> HeapDB<T>::HeapDB(AnalyticIE<T>* ie)
//     :itsAnalyticIE(ie)   
// {
//     assert(itsAnalyticIE);
// }

// template <class T> HeapDB<T>::~HeapDB()
// {
//     //Report(std::cout);
//     delete itsAnalyticIE;
// }

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

// template <class T> void HeapDB<T>::Report(std::ostream& os) const
// {
//     size_t N1=Size(its1C)+Size(its1Cx);
//     size_t N2=Size(its2C)+Size(its2Cx)+Size(its2CNuc);
//     size_t N3=Size(its3C);
//     //size_t N4=Size(Jac)+Size(Kab);
    
//     os << "Heap DB storage report:" << std::endl;
//     os << "    " << setw(10) << N1 << " 1 centre integrals." << std::endl;
//     os << "    " << setw(10) << N2 << " 2 centre integrals." << std::endl;
//     os << "    " << setw(10) << N3 << " 3 centre integrals." << std::endl;
//     //os << "    " << setw(10) << N4 << " 4 centre integrals." << std::endl;
    
// }
//---------------------------------------------------------------------------------

template <class T> typename Integrals_Base<T>::SMat_ref DB_Overlap<T>::Overlap() const
{
    auto cache(DB_Common<T>::itsCache);
    assert(cache);
    typename DB_cache<T>::id2c_t key=std::make_tuple(qchem::Overlap2C,this->GetID());
    if (auto i = cache->itsSMats.find(key); i==cache->itsSMats.end())
        return cache->itsSMats[key] = MakeOverlap();
    else
        return i->second;
}
template <class T> typename Integrals_Base<T>::SMat_ref DB_Kinetic<T>::Kinetic() const
{
    auto cache(DB_Common<T>::itsCache);
    assert(cache);
    typename DB_cache<T>::id2c_t key=std::make_tuple(qchem::Kinetic,this->GetID());
    if (auto i = cache->itsSMats.find(key); i==cache->itsSMats.end())
        return cache->itsSMats[key] = MakeKinetic();
    else
        return i->second;
}
template <class T> typename Integrals_Base<T>::SMat_ref DB_Nuclear<T>::Nuclear(const Cluster* cl) const
{
    auto cache(DB_Common<T>::itsCache);
    assert(cache);
    typename DB_cache<T>::id2c_t key=std::make_tuple(qchem::Nuclear,this->GetID());
    if (auto i = cache->itsSMats.find(key); i==cache->itsSMats.end())
        return cache->itsSMats[key] = MakeNuclear(cl);
    else
        return i->second;
}
template class DB_1E<double>;


template <class T> const typename DB_DFT<T>::ERI3& DB_DFT<T>::Overlap3C(const fbs_t& c) const
{ 
    auto cache(DB_Common<T>::itsCache);
    assert(cache);
    typename DB_cache<T>::id3c_t key=std::make_tuple(qchem::Overlap3C,this->GetID(),c.GetID());
    if (auto i = cache->itsERI3s.find(key); i==cache->itsERI3s.end())
    {
        return cache->itsERI3s[key] = MakeOverlap3C(c);
    }
    else
        return i->second;
}
template <class T> const typename DB_DFT<T>::ERI3& DB_DFT<T>::Repulsion3C(const fbs_t& c) const
{
    auto cache(DB_Common<T>::itsCache);
    assert(cache);
    typename DB_cache<T>::id3c_t key=std::make_tuple(qchem::Repulsion3C,this->GetID(),c.GetID());
    if (auto i = cache->itsERI3s.find(key); i==cache->itsERI3s.end())
    {
        return cache->itsERI3s[key] = MakeRepulsion3C(c);
    }
    else
        return i->second;
}
template class DB_DFT<double>;

template <class T> DB_2E<T>::DB_2E(const DB_BS_2E<T>* db) 
    : itsDB_BS_2E(db) 
    {
        assert(itsDB_BS_2E);
    };
template <class T> ERI4 DB_2E<T>::Direct(const obs_t& c) const
{
    assert(itsDB_BS_2E);
    return itsDB_BS_2E->Direct(this->GetID(),c.GetID());
}
template <class T> ERI4 DB_2E<T>::Exchange(const obs_t& b) const
{
    assert(itsDB_BS_2E);
    return itsDB_BS_2E->Exchange(this->GetID(),b.GetID()); 
}

template <class T> void DB_BS_2E<T>::Append(const IrrepIEClient* iec)
{
    itsIrreps.push_back(iec);
}
template <class T> ERI4 DB_BS_2E<T>::Direct(IDType a,IDType c) const
{
    assert(a<=c);
    if (Jac.size()==0) MakeDirect();
    //cout << "GetRepulsion4C_new a,c=" << a.GetIndex() << " " << c.GetIndex() << endl;
    assert(Jac.find(a)!=Jac.end());
    assert(Jac[a].find(c)!=Jac[a].end());
    
    return Jac[a][c];
}
template <class T> ERI4 DB_BS_2E<T>::Exchange(IDType a,IDType b) const
{
    assert(a<=b);
    if (Kab.size()==0) MakeExchange(); 
    //cout << "GetExchange4C_new a,b=" << a.GetIndex() << " " << b.GetIndex() << endl;
    assert(Kab.find(a)!=Kab.end());
    assert(Kab[a].find(b)!=Kab[a].end());
    
    return Kab[a][b];
}
template <class T> void DB_BS_2E<T>::MakeDirect() const
{
    Jac.clear();
    for (auto a: itsIrreps)
        for (auto c: itsIrreps) //TODO run from ia n
        {
            if (a->GetID()>c->GetID()) continue;
            Jac[a->GetID()][c->GetID()]=MakeDirect(a,c);
        }

}
template <class T> void DB_BS_2E<T>::MakeExchange() const
{
    Kab.clear();
    for (auto a: itsIrreps)
        for (auto b: itsIrreps) 
        {
            if (a->GetID()>b->GetID()) continue;
            Kab[a->GetID()][b->GetID()]=MakeExchange(a,b);            
        }
    
}
template class DB_2E<double>;
template class DB_BS_2E<double>;


template <class T> typename Integrals_Base<T>::SMat_ref  DB_RKB <T>::RestMass() const
{
    auto cache(DB_Overlap<T>::itsCache);
    assert(cache);
    typename DB_cache<T>::id2c_t key=std::make_tuple(qchem::RestMass,this->GetID());
    if (auto i = cache->itsSMats.find(key); i==cache->itsSMats.end())
        return cache->itsSMats[key] = MakeRestMass();
    else
        return i->second;
}
template <class T> typename Integrals_Base<T>::Mat_ref  DB_RKBL<T>::Kinetic(const Orbital_RKBS_IBS<T>* rkbs) const
{
    auto cache(DB_Overlap<T>::itsCache);
    assert(cache);
    typename DB_cache<T>::idx_t key=std::make_tuple(qchem::Kinetic,this->GetID(),rkbs->GetID());
    if (auto i = cache->itsMats.find(key); i==cache->itsMats.end())
        return cache->itsMats[key] = MakeKinetic(rkbs);
    else
        return i->second;
}
// template <class T> typename DB_RKBS<T>::SMat_ref DB_RKBS<T>::Overlap() const
// {
//     id2c_t key=std::make_tuple(qchem::Overlap2C,this->GetID());
//     if (auto i = itsBuffer.find(key); i==itsBuffer.end())
//     {
//         return itsBuffer[key] = MakeOverlap();
//     }
//     else
//         return i->second;
// }
// template <class T> typename DB_RKBS<T>::SMat_ref DB_RKBS<T>::Nuclear(const Cluster* cl) const
// {
//     assert(cl);
//     id2c_t key=std::make_tuple(qchem::Nuclear,this->GetID());
//     if (auto i = itsBuffer.find(key); i==itsBuffer.end())
//     {
//         return itsBuffer[key] = MakeNuclear(cl);
//     }
//     else
//         return i->second;
// }
// template <class T> typename DB_RKBS<T>::SMat_ref DB_RKBS<T>::RestMass() const
// {
//     auto cache(DB_Overlap<T>::itsCache);
//     assert(cache);
//     id2c_t key=std::make_tuple(qchem::RestMass,this->GetID());
//     if (auto i = cache->itsSMats.find(key); i==cache->itsSMats.end())
//         return cache->itsSMats[key] = MakeRestMass();
//     else
//         return i->second;
// }

#include "Imp/BasisSet/AtomIEClient.H"



template class DB_RKB<double>;
template class DB_RKBL<double>;
template class DB_RKBS<double>;


DB_Fit::Vec_ref DB_Fit::Charge() const
{
    assert(itsCache);
    DB_cache<double>::id2c_t key=std::make_tuple(qchem::Charge,this->GetID());
    if (auto i = itsCache->itsVecs.find(key); i==itsCache->itsVecs.end())
        return itsCache->itsVecs[key] = MakeCharge();
    else
        return i->second;
}
DB_Fit::SMat_ref DB_Fit::Repulsion() const
{
    DB_cache<double>::id2c_t key=std::make_tuple(qchem::Repulsion2C,this->GetID());
    if (auto i = itsCache->itsSMats.find(key); i==itsCache->itsSMats.end())
        return itsCache->itsSMats[key] = MakeRepulsion();
    else
        return i->second;
}
DB_Fit::Mat_ref DB_Fit::Repulsion(const fbs_t& b) const
{
    DB_cache<double>::idx_t key=std::make_tuple(qchem::Repulsion2C,this->GetID(),b.GetID());
    if (auto i = itsCache->itsMats.find(key); i==itsCache->itsMats.end())
        return itsCache->itsMats[key] = MakeRepulsion(b);
    else
        return i->second;
}
DB_Fit::SMat_ref DB_Fit::InvOverlap(const LAParams& lap) const
{
    DB_cache<double>::id2c_t key=std::make_tuple(qchem::InvOverlap,this->GetID());
    if (auto i = itsCache->itsSMats.find(key); i==itsCache->itsSMats.end())
        return itsCache->itsSMats[key] = MakeInverse(Overlap(),lap);
    else
        return i->second;
}
DB_Fit::SMat_ref DB_Fit::InvRepulsion(const LAParams& lap) const
{
    DB_cache<double>::id2c_t key=std::make_tuple(qchem::InvRepulsion,this->GetID());
    if (auto i = itsCache->itsSMats.find(key); i==itsCache->itsSMats.end())
        return itsCache->itsSMats[key] = MakeInverse(Repulsion(),lap);
    else
        return i->second;
}
DB_Fit::Vec_ref DB_Fit::Norm   (const Mesh* m        ) const
{
    DB_cache<double>::id2c_t key=std::make_tuple(qchem::NumNormalization,this->GetID());
    if (auto i = itsCache->itsVecs.find(key); i==itsCache->itsVecs.end())
        return itsCache->itsVecs[key] = MakeNorm(m);
    else
        return i->second;
}
DB_Fit::Vec_ref DB_Fit::Charge (const Mesh* m        ) const
{
    DB_cache<double>::id2c_t key=std::make_tuple(qchem::NumCharge,this->GetID());
    if (auto i = itsCache->itsVecs.find(key); i==itsCache->itsVecs.end())
        return itsCache->itsVecs[key] = MakeCharge(m);
    else
        return i->second;

}
DB_Fit::Mat_ref DB_Fit::Overlap(const Mesh* m,const fbs_t& b) const
{
    DB_cache<double>::idx_t key=std::make_tuple(qchem::NumOverlap,this->GetID(),b.GetID());
    if (auto i = itsCache->itsMats.find(key); i==itsCache->itsMats.end())
        return itsCache->itsMats[key] = MakeOverlap(m,b);
    else
        return i->second;
}


#include <LASolver.H>
#include <LAParams.H>
DB_Fit::SMat DB_Fit::MakeInverse(const SMat& S,const LAParams& lap) 
{
    LASolver<double>* las=LASolver<double>::Factory(lap);
    SMat Sinv=las->Inverse(S);
    delete las;
    return Sinv;
}
