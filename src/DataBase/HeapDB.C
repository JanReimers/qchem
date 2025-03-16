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

using std::cout;
using std::endl;

template <class T> typename DB_1E<T>::SMat_ref DB_1E<T>::Overlap() const
{
    assert(itsCache);
    id2c_t key=std::make_tuple(qchem::Overlap2C,this->GetID());
    if (auto i = itsCache->itsSMats.find(key); i==itsCache->itsSMats.end())
    {
        return itsCache->itsSMats[key] = MakeOverlap();
    }
    else
        return i->second;
}
template <class T> typename DB_1E<T>::SMat_ref DB_1E<T>::Kinetic() const
{
    id2c_t key=std::make_tuple(qchem::Kinetic,this->GetID());
    if (auto i = itsCache->itsSMats.find(key); i==itsCache->itsSMats.end())
    {
        return itsCache->itsSMats[key] = MakeKinetic();
    }
    else
        return i->second;
}
template <class T> typename DB_1E<T>::SMat_ref DB_1E<T>::Nuclear(const Cluster* cl) const
{
    assert(cl);
    id2c_t key=std::make_tuple(qchem::Nuclear,this->GetID());
    if (auto i = itsCache->itsSMats.find(key); i==itsCache->itsSMats.end())
    {
        return itsCache->itsSMats[key] = MakeNuclear(cl);
    }
    else
        return i->second;
}

template class DB_1E<double>;

DB_Fit::Vec_ref DB_Fit::Charge() const
{
    id2c_t key=std::make_tuple(qchem::Charge,this->GetID());
    if (auto i = itsCache->itsVecs.find(key); i==itsCache->itsVecs.end())
    {
        return itsCache->itsVecs[key] = MakeCharge();
    }
    else
        return i->second;
}

DB_Fit::SMat_ref DB_Fit::Overlap() const
{
    id2c_t key=std::make_tuple(qchem::Overlap2C,this->GetID());
    if (auto i = itsCache->itsSMats.find(key); i==itsCache->itsSMats.end())
    {
        return itsCache->itsSMats[key] = MakeOverlap();
    }
    else
        return i->second;
}

DB_Fit::SMat_ref DB_Fit::Repulsion() const
{
    id2c_t key=std::make_tuple(qchem::Repulsion2C,this->GetID());
    if (auto i = itsCache->itsSMats.find(key); i==itsCache->itsSMats.end())
    {
        return itsCache->itsSMats[key] = MakeRepulsion();
    }
    else
        return i->second;
}
DB_Fit::Mat_ref DB_Fit::Repulsion(const fbs_t& b) const
{
    idx_t key=std::make_tuple(qchem::Repulsion2C,this->GetID(),b.GetID());
    if (auto i = itsCache->itsMats.find(key); i==itsCache->itsMats.end())
    {
        return itsCache->itsMats[key] = MakeRepulsion(b);
    }
    else
        return i->second;
}
DB_Fit::SMat_ref DB_Fit::InvOverlap(const LAParams& lap) const
{
    id2c_t key=std::make_tuple(qchem::InvOverlap,this->GetID());
    if (auto i = itsCache->itsSMats.find(key); i==itsCache->itsSMats.end())
        return itsCache->itsSMats[key] = MakeInverse(Overlap(),lap);
    else
        return i->second;
}
DB_Fit::SMat_ref DB_Fit::InvRepulsion(const LAParams& lap) const
{
    id2c_t key=std::make_tuple(qchem::InvRepulsion,this->GetID());
    if (auto i = itsCache->itsSMats.find(key); i==itsCache->itsSMats.end())
        return itsCache->itsSMats[key] = MakeInverse(Repulsion(),lap);
    else
        return i->second;
}
DB_Fit::Vec_ref DB_Fit::Norm   (const Mesh* m        ) const
{
    id2c_t key=std::make_tuple(qchem::NumNormalization,this->GetID());
    if (auto i = itsCache->itsVecs.find(key); i==itsCache->itsVecs.end())
        return itsCache->itsVecs[key] = MakeNorm(m);
    else
        return i->second;
}
DB_Fit::Vec_ref DB_Fit::Charge (const Mesh* m        ) const
{
    id2c_t key=std::make_tuple(qchem::NumCharge,this->GetID());
    if (auto i = itsCache->itsVecs.find(key); i==itsCache->itsVecs.end())
        return itsCache->itsVecs[key] = MakeCharge(m);
    else
        return i->second;

}
DB_Fit::Mat_ref DB_Fit::Overlap(const Mesh* m,const fbs_t& b) const
{
    idx_t key=std::make_tuple(qchem::NumOverlap,this->GetID(),b.GetID());
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

template <class T> const typename DB_DFT<T>::ERI3& DB_DFT<T>::Overlap3C(const fbs_t& c) const
{
    id3c_t key=std::make_tuple(qchem::Overlap3C,GetID(),c.GetID());
    if (auto i = itsCache->itsERI3s.find(key); i==itsCache->itsERI3s.end())
    {
        return itsCache->itsERI3s[key] = MakeOverlap3C(c);
    }
    else
        return i->second;
}
template <class T> const typename DB_DFT<T>::ERI3& DB_DFT<T>::Repulsion3C(const fbs_t& c) const
{
    id3c_t key=std::make_tuple(qchem::Repulsion3C,GetID(),c.GetID());
    if (auto i = itsCache->itsERI3s.find(key); i==itsCache->itsERI3s.end())
    {
        return itsCache->itsERI3s[key] = MakeRepulsion3C(c);
    }
    else
        return i->second;
}
template class DB_DFT<double>;

template <class T> ERI4 DB_2E<T>::Direct(const obs_t& c) const
{
    assert(itsDB_BS_2E);
    return itsDB_BS_2E->Direct(GetID(),c.GetID());
}
template <class T> ERI4 DB_2E<T>::Exchange(const obs_t& b) const
{
    assert(itsDB_BS_2E);
    return itsDB_BS_2E->Exchange(GetID(),b.GetID()); 
}

template <class T> DB_2E<T>::DB_2E(const DB_BS_2E<T>* db) 
    : itsDB_BS_2E(db) 
    {
        assert(itsDB_BS_2E);
    };

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

template <class T> typename DB_RKB<T>::SMat_ref DB_RKB<T>::Kinetic() const
{
    id2c_t key=std::make_tuple(qchem::Kinetic,this->GetID());
    if (auto i = itsBuffer.find(key); i==itsBuffer.end())
    {
        return itsBuffer[key] = MakeKinetic();
    }
    else
        return i->second;
}
template <class T> typename DB_RKB<T>::SMat_ref DB_RKB<T>::RestMass() const
{
    id2c_t key=std::make_tuple(qchem::RestMass,this->GetID());
    if (auto i = itsBuffer.find(key); i==itsBuffer.end())
    {
        return itsBuffer[key] = MakeRestMass();
    }
    else
        return i->second;
}
template <class T> typename DB_RKBL<T>::SMat_ref DB_RKBL<T>::Overlap() const
{
    id2c_t key=std::make_tuple(qchem::Kinetic,this->GetID());
    if (auto i = itsBuffer.find(key); i==itsBuffer.end())
    {
        return itsBuffer[key] = MakeOverlap();
    }
    else
        return i->second;
}
template <class T> typename DB_RKBL<T>::Mat_ref DB_RKBL<T>::Kinetic(const Orbital_RKBS_IBS<T>* rkbs) const
{
    id2c_t key=std::make_tuple(qchem::Kinetic,this->GetID());
    if (auto i = itsBufferX.find(key); i==itsBufferX.end())
    {
        return itsBufferX[key] = MakeKinetic(rkbs);
    }
    else
        return i->second;
}
template <class T> typename DB_RKBL<T>::SMat_ref DB_RKBL<T>::Nuclear(const Cluster* cl)  const
{
    id2c_t key=std::make_tuple(qchem::Nuclear,this->GetID());
    if (auto i = itsBuffer.find(key); i==itsBuffer.end())
    {
        return itsBuffer[key] = MakeNuclear(cl);
    }
    else
        return i->second;
}
template <class T> typename DB_RKBS<T>::SMat_ref DB_RKBS<T>::Overlap() const
{
    id2c_t key=std::make_tuple(qchem::Overlap2C,this->GetID());
    if (auto i = itsBuffer.find(key); i==itsBuffer.end())
    {
        return itsBuffer[key] = MakeOverlap();
    }
    else
        return i->second;
}
template <class T> typename DB_RKBS<T>::SMat_ref DB_RKBS<T>::Nuclear(const Cluster* cl) const
{
    assert(cl);
    id2c_t key=std::make_tuple(qchem::Nuclear,this->GetID());
    if (auto i = itsBuffer.find(key); i==itsBuffer.end())
    {
        return itsBuffer[key] = MakeNuclear(cl);
    }
    else
        return i->second;
}
template <class T> typename DB_RKBS<T>::SMat_ref DB_RKBS<T>::RestMass() const
{
    id2c_t key=std::make_tuple(qchem::RestMass,this->GetID());
    if (auto i = itsBuffer.find(key); i==itsBuffer.end())
    {
        return itsBuffer[key] = MakeRestMass();
    }
    else
        return i->second;
}
template <class T> typename DB_RKBL<T>::SMat DB_RKBL<T>::MakeOverlap() const
{
    return this->MakeIntegrals(qchem::Overlap1);
}
template <class T> typename DB_RKBL<T>::SMat DB_RKBL<T>::MakeNuclear(const Cluster* cl) const
{
    assert(cl);
    int Z=cl->GetNuclearCharge();
    return -Z*this->MakeIntegrals(qchem::Nuclear1,cl);
}


#include "Imp/BasisSet/AtomIEClient.H"

template <class T> typename DB_RKBL<T>::SMat  DB_RKBL<T>::MakeIntegrals(qchem::IType it,const Cluster* cl)  const
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this);
    assert(a);
    size_t N=a->size(),l=a->l;
    SMatrix<double> H(N);
    for (auto i:H.rows())
        for (auto j:H.cols(i))
            H(i,j)= Integral(it,a->es(i),a->es(j),l)*a->ns(i)*a->ns(j);

    return H;
}
template <class T> typename DB_RKBS<T>::SMat DB_RKBS<T>::MakeOverlap() const
{
    return this->MakeIntegrals(qchem::Overlap1);
}
template <class T> typename DB_RKBL<T>::Mat  DB_RKBL<T>::MakeKinetic(const Orbital_RKBS_IBS<T>* rkbs) const
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this);
    const AtomIrrepIEClient* b=dynamic_cast<const AtomIrrepIEClient*>(rkbs);
    assert(a->l==b->l);
    size_t Na=a->size();
    size_t Nb=b->size();
    Matrix<double> Hk(Na,Nb);
    for (auto i:Hk.rows())
        for (auto j:Hk.cols())
            Hk(i,j)=Integral(qchem::Kinetic1,a->es(i),b->es(j),a->l)*a->ns(i)*b->ns(j);

    return Hk;
}
template <class T> typename DB_RKBS<T>::SMat DB_RKBS<T>::MakeNuclear(const Cluster* cl) const
{
    assert(cl);
    int Z=cl->GetNuclearCharge();
    return -Z*this->MakeIntegrals(qchem::Nuclear1,cl);
}
template <class T> typename DB_RKBS<T>::SMat DB_RKBS<T>::MakeRestMass() const
{
    return this->MakeIntegrals(qchem::RestMass1);
}
template <class T> typename DB_RKBS<T>::SMat  DB_RKBS<T>::MakeIntegrals(qchem::IType it,const Cluster* cl)  const
{
    const AtomIrrepIEClient* a=dynamic_cast<const AtomIrrepIEClient*>(this);
    assert(a);
    size_t N=a->size(),l=a->l;
    SMatrix<double> H(N);
    for (auto i:H.rows())
        for (auto j:H.cols(i))
            H(i,j)= Integral(it,a->es(i),a->es(j),l)*a->ns(i)*a->ns(j);

    return H;
}

template class DB_RKB<double>;
template class DB_RKBL<double>;
template class DB_RKBS<double>;