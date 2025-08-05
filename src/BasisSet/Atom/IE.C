// File: AtomIE.C Common IE code for all atom basis sets.
module;
#include <memory>
export module qchem.BasisSet.Atom.IE;
export import qchem.BasisSet.Internal.HeapDB;
export import qchem.BasisSet.Internal.Cache4;
export import oml.Vector;
export import qchem.Orbital_DHF_IBS;
export import qchem.Types;
export import qchem.BasisSet.Atom.IE_Primatives;
import qchem.BasisSet.Atom.Internal.BFGrouper;
import qchem.BasisSet.Atom.IEClient;

import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet.Internal.IEClient;
import qchem.Fit_IBS;
import qchem.Orbital_DFT_IBS;

export
{

template <class T> class AtomIE_Overlap
: public DB_Overlap<T>
{
protected:
    AtomIE_Overlap(const DB_cache<T>* db,const IE_Primatives* _pie) : DB_Overlap<T>(db), pie(_pie) {};
    virtual SMatrix<T> MakeOverlap() const;
private:
    const IE_Primatives* pie;
};
template <class T> class AtomIE_Kinetic
: public DB_Kinetic<T>
{
protected:
    AtomIE_Kinetic(const DB_cache<T>* db,const IE_Primatives* _pie) : DB_Kinetic<T>(db), pie(_pie) {};
    virtual SMatrix<T> MakeKinetic() const;
private:
    const IE_Primatives* pie;
};
template <class T> class AtomIE_Nuclear
: public DB_Nuclear<T>
{
protected:
    virtual SMatrix<T> MakeNuclear(const Cluster* cl) const;
    AtomIE_Nuclear(const DB_cache<T>* db,const IE_Primatives* _pie) : DB_Nuclear<T>(db), pie(_pie) {};
private:
    const IE_Primatives* pie;
};
template <class T> class AtomIE_XKinetic
: public DB_XKinetic<T>
{
protected:
    virtual Matrix<T> MakeKinetic(const Orbital_RKBS_IBS<T>* rkbs) const;
    AtomIE_XKinetic(const DB_cache<T>* db,const IE_Primatives* _pie) : DB_XKinetic<T>(db), pie(_pie) {};
private:
    const IE_Primatives* pie;
};

// HF
class AtomIE_BS_2E_Angular
{
public:
    virtual ~AtomIE_BS_2E_Angular() {};
    typedef AtomIrrepIEClient iec_t;
    virtual RVec Coulomb_AngularIntegrals(const iec_t* a,const iec_t* c) const=0;
    virtual RVec ExchangeAngularIntegrals(const iec_t* a,const iec_t* c) const=0;
};

template <class T> class AtomIE_BS_2E 
    : public virtual Cache4
    , public DB_BS_2E<T>
    , public BFGrouper
{
public:
    AtomIE_BS_2E(AtomIE_BS_2E_Angular* a) : itsAngular(a) {};
    virtual ERI4 MakeDirect  (const IrrepIEClient* a, const IrrepIEClient* c) const;
    virtual ERI4 MakeExchange(const IrrepIEClient* a, const IrrepIEClient* c) const;

    // Cach4 functions
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc) const=0;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc) const=0;
protected:
    virtual void Append(const IrrepIEClient*);
private: 
    std::unique_ptr<AtomIE_BS_2E_Angular> itsAngular;
};

// DFT
template <class T> class AtomIE_DFT 
: public DB_DFT<T>
{
protected:
    AtomIE_DFT(const DB_cache<T>* db,const IE_Primatives* _pie) : DB_DFT<T>(db), pie(_pie) {};
    
    virtual ERI3<T> MakeOverlap3C  (const Fit_IBS& c) const;
    virtual ERI3<T> MakeRepulsion3C(const Fit_IBS& c) const;
private:
    typedef AtomIrrepIEClient::bf_tuple bf_tuple;
    SMatrix<T> MakeOverlap  (const bf_tuple& c) const; //ab loops
    SMatrix<T> MakeRepulsion(const bf_tuple& c) const; //ab loops
private:
    const IE_Primatives* pie;
};
// DHF
template <class T> class AtomIE_RKBL 
    : public AtomIE_Overlap<T>
    , public AtomIE_XKinetic<T>
    , public AtomIE_Nuclear<T>
{
protected:
    AtomIE_RKBL(const DB_cache<T>* db,const IE_Primatives* pie) 
    : AtomIE_Overlap<T>(db,pie)
    , AtomIE_XKinetic<T>(db,pie)
    , AtomIE_Nuclear<T>(db,pie) 
    {};

};
template <class T> class AtomIE_RKBS 
: public AtomIE_Kinetic<T>
, public AtomIE_Nuclear<T>
{
protected:
    AtomIE_RKBS(const DB_cache<T>* db,const ::IE_Primatives* pie) : AtomIE_Kinetic<T>(db,pie), AtomIE_Nuclear<T>(db,pie) {};
};
// Fit
class AtomIE_Fit 
: public DB_Fit
{
    protected:
    AtomIE_Fit(const DB_cache<double>* db,const IE_Primatives* _pie) : DB_Fit(db), pie(_pie) {};

    virtual  Vector<double> MakeCharge() const;
    virtual SMatrix<double> MakeRepulsion() const;
    virtual  Matrix<double> MakeRepulsion(const Fit_IBS&) const;
private:
    using DB_Fit::Charge; //un hide
    using DB_Fit::Repulsion; //un hide
private:
    const IE_Primatives* pie;

};

} // export block