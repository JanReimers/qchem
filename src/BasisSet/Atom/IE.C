// File: AtomIE.C Common IE code for all atom basis sets.
export module qchem.BasisSet.Atom.IE;
export import qchem.BasisSet.Internal.HeapDB;
export import qchem.BasisSet.Internal.Cache4;
export import oml.Vector;
export import qchem.DHF_IBS;
export import qchem.Types;
import qchem.BasisSet.Atom.Internal.BFGrouper;
import qchem.BasisSet.Atom.IEClient;

import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.IEClient;
import qchem.Fit_IBS;
import qchem.DFT_IBS;

export
{

//  Generic
template <class T> class Primative_Overlap
{
public:
    virtual double Overlap(double ea ,double eb,size_t l_total) const=0;
};
template <class T> class Primative_Grad2
{
public:
    virtual double Grad2(double ea ,double eb,size_t la, size_t lb) const=0;
};
template <class T> class Primative_Inv_r1
{
public:
    virtual double Inv_r1(double ea ,double eb,size_t l_total) const=0;
};
template <class T> class Primative_Inv_r2
{
public:
    virtual double Inv_r2(double ea ,double eb,size_t l_total) const=0;
};
template <class T> class Primative_Repulsion
{
public:
    virtual double Repulsion(double ea ,double ec,size_t la, size_t lc) const=0;
};
template <class T> class Primative_Charge
{
public:
    virtual double Charge   (double ea, size_t l) const=0;
};

template <class T> class AtomIE_Overlap
: public virtual Primative_Overlap<T>
, public DB_Overlap<T>
{
protected:
    using Primative_Overlap<T>::Overlap;
    virtual SMatrix<T> MakeOverlap() const;
    AtomIE_Overlap(const DB_cache<T>* db) : DB_Overlap<T>(db) {};
};
template <class T> class AtomIE_Kinetic
: public virtual Primative_Grad2<T>
, public virtual Primative_Inv_r2<T> //for centrifugal term.
, public DB_Kinetic<T>
{
protected:
    using Primative_Grad2 <T>::Grad2;
    using Primative_Inv_r2<T>::Inv_r2;
    virtual SMatrix<T> MakeKinetic() const;
    AtomIE_Kinetic(const DB_cache<T>* db) : DB_Kinetic<T>(db) {};
};
template <class T> class AtomIE_Nuclear
: public virtual Primative_Inv_r1<T>
, public DB_Nuclear<T>
{
protected:
    using Primative_Inv_r1<T>::Inv_r1;
    virtual SMatrix<T> MakeNuclear(const Cluster* cl) const;
    AtomIE_Nuclear(const DB_cache<T>* db) : DB_Nuclear<T>(db) {};
};
template <class T> class AtomIE_XKinetic
: public virtual Primative_Grad2<T>
, public virtual Primative_Inv_r2<T>
, public DB_XKinetic<T>
{
protected:
    using Primative_Grad2<T>::Grad2;
    using Primative_Inv_r2<T>::Inv_r2;
    virtual Matrix<T> MakeKinetic(const Orbital_RKBS_IBS<T>* rkbs) const;
    AtomIE_XKinetic(const DB_cache<T>* db) : DB_XKinetic<T>(db) {};
};

// HF
class AtomIE_BS_2E_Angular
{
public:
    typedef AtomIrrepIEClient iec_t;
    virtual RVec Coulomb_AngularIntegrals(const iec_t* a,const iec_t* c) const=0;
    virtual RVec ExchangeAngularIntegrals(const iec_t* a,const iec_t* c) const=0;
};

template <class T> class AtomIE_BS_2E 
    : public virtual AtomIE_BS_2E_Angular
    , public virtual Cache4
    , public DB_BS_2E<T>
    , public BFGrouper
{
public:
    virtual ERI4 MakeDirect  (const IrrepIEClient* a, const IrrepIEClient* c) const;
    virtual ERI4 MakeExchange(const IrrepIEClient* a, const IrrepIEClient* c) const;

    // Cach4 functions
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc) const=0;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc) const=0;
protected:
    virtual void Append(const IrrepIEClient*);
};

// DFT
template <class T> class AtomIE_DFT 
: public virtual Primative_Overlap<T>
, public virtual Primative_Repulsion<T>
, public DB_DFT<T>
{
protected:
    AtomIE_DFT(const DB_cache<T>* db) : DB_DFT<T>(db) {};
    
    virtual ERI3<T> MakeOverlap3C  (const Fit_IBS& c) const;
    virtual ERI3<T> MakeRepulsion3C(const Fit_IBS& c) const;
private:
    typedef AtomIrrepIEClient::bf_tuple bf_tuple;
    SMatrix<T> MakeOverlap  (const bf_tuple& c) const; //ab loops
    SMatrix<T> MakeRepulsion(const bf_tuple& c) const; //ab loops
};
// DHF
template <class T> class AtomIE_RKBL 
    : public AtomIE_Overlap<T>
    , public AtomIE_XKinetic<T>
    , public AtomIE_Nuclear<T>
{
protected:
    AtomIE_RKBL(const DB_cache<T>* db) : AtomIE_Overlap<T>(db),AtomIE_XKinetic<T>(db),AtomIE_Nuclear<T>(db) {};

};
template <class T> class AtomIE_RKBS 
: public AtomIE_Kinetic<T>
, public AtomIE_Nuclear<T>
{
protected:
    AtomIE_RKBS(const DB_cache<T>* db) : AtomIE_Kinetic<T>(db), AtomIE_Nuclear<T>(db) {};
};
// Fit
class AtomIE_Fit 
: public virtual Primative_Repulsion<double>
, public virtual Primative_Charge<double>
, public DB_Fit
{
    protected:
    AtomIE_Fit(const DB_cache<double>* db) : DB_Fit(db) {};

    virtual  Vector<double> MakeCharge() const;
    virtual SMatrix<double> MakeRepulsion() const;
    virtual  Matrix<double> MakeRepulsion(const Fit_IBS&) const;
private:
    // Derived classes must provide the actual integral calculations.
    using DB_Fit::Charge; //un hide
    using DB_Fit::Repulsion; //un hide
    using Primative_Repulsion<double>::Repulsion;
    using Primative_Charge   <double>::Charge;
};

} // export block