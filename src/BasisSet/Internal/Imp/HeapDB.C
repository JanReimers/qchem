// File: HeapDBImp.C  Implement a heap storage integral data base.
module;
#include <iomanip>
#include <vector>
#include <memory>
#include <map>
#include <cassert>
#include "blaze/Math.h" 

module qchem.BasisSet.Internal.HeapDB;
import qchem.Fit_IBS;
import qchem.Orbital_DHF_IBS;
import qchem.IrrepBasisSet;
import qchem.Mesh.Integrator;
import qchem.Cluster;
import qchem.BasisSet.Internal.ERI4;

// template <class T> size_t Size(const Vector <T>& m) {return m.size();}
// template <class T> size_t Size(const Matrix <T>& m) {return m.size();}
// template <class T> size_t Size(const SMatrix<T>& m) {return m.size();}
//                    size_t Size(const ERI4      & m) {return m.size();}
// template <class M> size_t Size(const std::vector<M> & v) 
// {
//     size_t N=0;
//     for (auto i:v) 
//         N+=Size(i);
//     return N;
// }

// template <class K, class M> size_t Size(const std::map<K,M>& m)
// {
//     size_t N=0;
//     for (auto i:m) 
//         N+=Size(i.second);
//     return N;
// }

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

template <class T> const smat_t<T>& DB_Overlap <T>::Overlap() const
{
    auto cache(DB_Common<T>::itsCache);
    assert(cache);
    typename DB_cache<T>::id2c_t key=std::make_tuple(qchem::Overlap2C,this->GetID());
    if (auto i = cache->itsbSMats.find(key); i==cache->itsbSMats.end())
        return cache->itsbSMats[key] = MakeOverlap();
    else
        return i->second;
}
template <class T> const smat_t<T>& DB_Kinetic <T>::Kinetic() const
{
    auto cache(DB_Common<T>::itsCache);
    assert(cache);
    typename DB_cache<T>::id2c_t key=std::make_tuple(qchem::Grad2,this->GetID());
    if (auto i = cache->itsbSMats.find(key); i==cache->itsbSMats.end())
        return cache->itsbSMats[key] = MakeKinetic();
    else
        return i->second;
}
template <class T> const smat_t<T>& DB_Nuclear <T>::Nuclear(const Cluster* cl) const
{
    auto cache(DB_Common<T>::itsCache);
    assert(cache);
    typename DB_cache<T>::id2c_t key=std::make_tuple(qchem::Nuclear,this->GetID());
    if (auto i = cache->itsbSMats.find(key); i==cache->itsbSMats.end())
        return cache->itsbSMats[key] = MakeNuclear(cl);
    else
        return i->second;
}


template <class T> const mat_t<T>& DB_XKinetic<T>::Kinetic(const Orbital_RKBS_IBS<T>* rkbs) const
{
    auto cache(DB_Common<T>::itsCache);
    assert(cache);
    typename DB_cache<T>::idx_t key=std::make_tuple(qchem::Grad2,this->GetID(),rkbs->GetID());
    if (auto i = cache->itsbMats.find(key); i==cache->itsbMats.end())
        return cache->itsbMats[key] = MakeKinetic(rkbs);
    else
        return i->second;
}
template <class T> const smat_t<T>& DB_RestMass<T>::RestMass() const
{
    auto cache(DB_Common<T>::itsCache);
    assert(cache);
    typename DB_cache<T>::id2c_t key=std::make_tuple(qchem::RestMass,this->GetID());
    if (auto i = cache->itsbSMats.find(key); i==cache->itsbSMats.end())
        return cache->itsbSMats[key] = MakeRestMass();
    else
        return i->second;
}


template <class T> const ERI3<T>& DB_DFT<T>::Overlap3C(const Fit_IBS& c) const
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
template <class T> const ERI3<T>& DB_DFT<T>::Repulsion3C(const Fit_IBS& c) const
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


template <class T> DB_HF<T>::DB_HF(const Integrals_BS_HF<T>* db) 
    : itsDB_BS_HF(db) 
    {
        assert(itsDB_BS_HF);
    };
template <class T> ERI4 DB_HF<T>::Direct(const Orbital_HF_IBS<T>& c) const
{
    assert(itsDB_BS_HF);
    return itsDB_BS_HF->Direct(this->GetID(),c.GetID());
}
template <class T> ERI4 DB_HF<T>::Exchange(const Orbital_HF_IBS<T>& b) const
{
    assert(itsDB_BS_HF);
    return itsDB_BS_HF->Exchange(this->GetID(),b.GetID()); 
}

template class DB_HF<double>;
template class DB_RKB<double>;
template class DB_RKBL<double>;
template class DB_RKBS<double>;

const rvec_t& DB_Fit::Charge   () const
{
    assert(itsCache);
    DB_cache<double>::id2c_t key=std::make_tuple(qchem::Charge,this->GetID());
    if (auto i = itsCache->itsbVecs.find(key); i==itsCache->itsbVecs.end())
        return itsCache->itsbVecs[key] = MakeCharge();
    else
        return i->second;
}
const rsmat_t& DB_Fit::Repulsion() const
{
    DB_cache<double>::id2c_t key=std::make_tuple(qchem::Repulsion2C,this->GetID());
    if (auto i = itsCache->itsbSMats.find(key); i==itsCache->itsbSMats.end())
        return itsCache->itsbSMats[key] = MakeRepulsion();
    else
        return i->second;
}
const rmat_t& DB_Fit::Repulsion(const Fit_IBS& b) const
{
    DB_cache<double>::idx_t key=std::make_tuple(qchem::Repulsion2C,this->GetID(),b.GetID());
    if (auto i = itsCache->itsbMats.find(key); i==itsCache->itsbMats.end())
        return itsCache->itsbMats[key] = MakeRepulsion(b);
    else
        return i->second;
}
const rsmat_t& DB_Fit::InvOverlap  () const
{
    DB_cache<double>::id2c_t key=std::make_tuple(qchem::InvOverlap,this->GetID());
    if (auto i = itsCache->itsbSMats.find(key); i==itsCache->itsbSMats.end())
        return itsCache->itsbSMats[key] = MakeInverse(Overlap());
    else
        return i->second;
}
const rsmat_t& DB_Fit::InvRepulsion() const
{
    DB_cache<double>::id2c_t key=std::make_tuple(qchem::InvRepulsion,this->GetID());
    if (auto i = itsCache->itsbSMats.find(key); i==itsCache->itsbSMats.end())
        return itsCache->itsbSMats[key] = MakeInverse(Repulsion());
    else
        return i->second;
}
const rvec_t& DB_Fit::Norm   (const Mesh* m        ) const
{
    DB_cache<double>::id2c_t key=std::make_tuple(qchem::NumNormalization,this->GetID());
    if (auto i = itsCache->itsbVecs.find(key); i==itsCache->itsbVecs.end())
        return itsCache->itsbVecs[key] = MakeNorm(m);
    else
        return i->second;
}
// const rvec_t& DB_Fit::Charge (const Mesh* m        ) const
// {
//     DB_cache<double>::id2c_t key=std::make_tuple(qchem::NumCharge,this->GetID());
//     if (auto i = itsCache->itsbVecs.find(key); i==itsCache->itsbVecs.end())
//         return itsCache->itsbVecs[key] = MakeCharge(m);
//     else
//         return i->second;
// }

const rmat_t& DB_Fit::Overlap(const Mesh* m,const Fit_IBS& b) const
{
    DB_cache<double>::idx_t key=std::make_tuple(qchem::NumOverlap,this->GetID(),b.GetID());
    if (auto i = itsCache->itsbMats.find(key); i==itsCache->itsbMats.end())
        return itsCache->itsbMats[key] = MakeOverlap(m,b);
    else
        return i->second;
}


rsmat_t DB_Fit::MakeInverse(const rsmat_t& S) 
{
    return inv(S);
}
