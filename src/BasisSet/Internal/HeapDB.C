// File: HeapDB.C  Implement an integral data base that uses heap for storage.
module;
#include <map>
#include <cassert>
#include <vector>
export module qchem.BasisSet.Internal.HeapDB;
import qchem.LAParams;
import qchem.BasisSet.Internal.IEClient;

import qchem.BasisSet.Internal.IntegralEnums;
import qchem.BasisSet.Internal.ERI4;

import qchem.Orbital_DHF_IBS;
import qchem.Orbital_HF_IBS;
import qchem.Fit_IBS;
import qchem.Orbital_DFT_IBS;

import Common.UniqueID;




export template  <class T> class DB_cache  
{
    typedef UniqueID::IDtype IDType;
    typedef SMatrix<T> SMat;    
    typedef Matrix<T> Mat;    
    typedef Vector<T> Vec;    
public:
    typedef std::map<IDType,std::map<IDType,ERI4> > erij_t;
    typedef std::tuple<qchem::IType2C,IDType> id2c_t;
    typedef std::tuple<qchem::IType2C,IDType,IDType> idx_t;
    typedef std::tuple<qchem::IType3C,IDType,IDType> id3c_t;
    
    mutable std::map<id2c_t ,SMat> itsSMats; 
    mutable std::map< idx_t , Mat> itsMats; 
    mutable std::map<id2c_t , Vec> itsVecs; 
    mutable std::map<id3c_t ,ERI3<T>> itsERI3s; 
    mutable erij_t Jac,Kab;
};
 
export template <class T> class DB_Common 
    : virtual public UniqueID
{
protected:
    DB_Common(const DB_cache<T>* db) : itsCache(db) {assert(itsCache);}
    const DB_cache<T>* itsCache;
};

// It would be nice use virtual inheretance from DB_Common<T>, but that opens a contructor hornets nest.
export template <class T> class DB_Overlap  : public DB_Common<T>, public virtual Integrals_Overlap<T>
{    
protected:
    DB_Overlap(const DB_cache<T>* db) : DB_Common<T>(db) {};
    virtual const SMatrix<T>& Overlap() const;
    virtual SMatrix<T> MakeOverlap() const=0;
};
export template <class T> class DB_Kinetic    : public DB_Common<T>, public virtual Integrals_Kinetic<T>
{    
protected:
    DB_Kinetic(const DB_cache<T>* db) : DB_Common<T>(db) {};
    virtual const SMatrix<T>& Kinetic() const;
    virtual SMatrix<T> MakeKinetic() const=0;
};
export template <class T> class DB_Nuclear  : public DB_Common<T>, public virtual Integrals_Nuclear<T>
{    
protected:
    DB_Nuclear(const DB_cache<T>* db) : DB_Common<T>(db) {};
    virtual const SMatrix<T>& Nuclear(const Cluster*) const;
    virtual SMatrix<T> MakeNuclear(const Cluster*) const=0;
};
export template <class T> class DB_XKinetic   : public DB_Common<T>, public virtual Integrals_XKinetic<T>
{    
protected:
    DB_XKinetic(const DB_cache<T>* db) : DB_Common<T>(db) {};
    virtual const Matrix<T>& Kinetic(const Orbital_RKBS_IBS<T>* rkbs) const;
    virtual Matrix<T> MakeKinetic(const Orbital_RKBS_IBS<T>* rkbs) const=0;
};
export template <class T> class DB_RestMass : public DB_Common<T>, public virtual Integrals_RestMass<T>
{    
protected:
    DB_RestMass(const DB_cache<T>* db) : DB_Common<T>(db) {};
    virtual const SMatrix<T>& RestMass() const;
    virtual SMatrix<T> MakeRestMass() const=0;
};

export template <class T> class DB_DFT 
    : virtual public Integrals_DFT<T>
    , public DB_Common<T>
{
protected:
    DB_DFT(const DB_cache<T>* db) : DB_Common<T>(db) {}
    virtual const ERI3<T>& Overlap3C  (const Fit_IBS& c) const; //<ab|c>
    virtual const ERI3<T>& Repulsion3C(const Fit_IBS& c) const; //<a(1)b(1)|1/r12|c(2)>
    virtual ERI3<T> MakeOverlap3C  (const Fit_IBS& c) const=0;
    virtual ERI3<T> MakeRepulsion3C(const Fit_IBS& c) const=0;
};

export class DB_Fit 
    : public virtual FitIntegrals
    , public DB_Common<double>
{
protected:
    DB_Fit(const DB_cache<double>* db) : DB_Common<double>(db) {};

    using Integrals_Overlap<double>::Overlap; 
    virtual const Vector<double>&  Charge   () const;   
    virtual const SMatrix<double>& Repulsion() const;
    virtual const  Matrix<double>& Repulsion(const Fit_IBS&) const;
    virtual const SMatrix<double>& InvOverlap(const LAParams&) const;
    virtual const SMatrix<double>& InvRepulsion(const LAParams&) const;
    virtual  const Vector<double>& Norm   (const Mesh*        ) const; //Numerical .
    virtual  const Vector<double>& Charge (const Mesh*        ) const; //Numerical .
    virtual  const  Matrix<double>& Overlap(const Mesh*,const Fit_IBS& b) const; //Numerical X overlap.


private:
    // One time calls to un-buffered integral calculations.
    using FitIntegrals::MakeCharge;
    virtual Vector<double>  MakeCharge() const=0;
    // virtual SMatrix<double> MakeOverlap() const=0;
    virtual SMatrix<double> MakeRepulsion() const=0;
    virtual  Matrix<double> MakeRepulsion(const Fit_IBS&) const=0;
    
    //! \brief Return the Penrose inverse of a symmetric matrix using SVD decomposition
    //! If \f$ S=UsV^{\dagger} \f$, then \f$ S^{-1}=V\frac{1}{s}U^{\dagger} \f$
    static  SMatrix<double> MakeInverse  (const SMatrix<double>&,const LAParams&); //Numerically stable algo required.

};

export template <class T> class DB_BS_2E : public virtual Integrals_BS_2E<T>, public DB_cache<T>
{
    typedef UniqueID::IDtype IDType;   
public: 
    virtual ERI4 Direct  (IDType a,IDType c) const;
    virtual ERI4 Exchange(IDType a,IDType b) const;
protected:
    void Append(const Orbital_HF_IBS<T>*);
    virtual ERI4 MakeDirect  (const IrrepIEClient* a, const IrrepIEClient* c) const=0;
    virtual ERI4 MakeExchange(const IrrepIEClient* a, const IrrepIEClient* b) const=0;
private:
    //! Internally called once to build direct and exchange supermatrix tables.
    virtual void MakeDirect  () const; 
    virtual void MakeExchange() const; 
    std::vector<const IrrepIEClient*> itsIrreps; //Used for 2-electron integrals.
    using DB_cache<T>::Jac;
    using DB_cache<T>::Kab;
};
export template <class T> class DB_2E 
    : virtual public Integrals_HF<T>
    , virtual public UniqueID
{
public:
    typedef typename Integrals_HF<T>::obs_t obs_t;
    virtual ERI4 Direct  (const obs_t& c) const;
    virtual ERI4 Exchange(const obs_t& b) const;
protected:
    DB_2E() : itsDB_BS_2E(0) {};
    DB_2E(const DB_BS_2E<T>* db);
private:
    const DB_BS_2E<T>* itsDB_BS_2E; //Database of all supermatrix tables.
};


export template <class T> class DB_RKB 
    : virtual public Integrals_RKB<T>
    , public DB_Overlap<T>
    , public DB_Kinetic<T>
    , public DB_Nuclear<T>
    , public DB_RestMass<T>
{    
protected:
    DB_RKB(const DB_cache<T>* db) :  DB_Overlap<T>(db), DB_Kinetic<T>(db), DB_Nuclear<T>(db), DB_RestMass<T>(db) {}; 
   
};
export template <class T> class DB_RKBL 
    : virtual public Integrals_RKBL<T>
    , public DB_Overlap<T>
    , public DB_XKinetic<T>
    , public DB_Nuclear<T>
{
protected:
    DB_RKBL(const DB_cache<T>* db) : DB_Overlap<T>(db), DB_XKinetic<T>(db), DB_Nuclear<T>(db) {};
};
export template <class T> class DB_RKBS 
    : virtual public Integrals_RKBS<T>
    , public DB_Kinetic<T>
    , public DB_Nuclear<T>
{
protected:
    DB_RKBS(const DB_cache<T>* db) :  DB_Kinetic<T>(db), DB_Nuclear<T>(db) {};  
   
};
