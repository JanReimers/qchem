// File: BSpline/IE.C Common IE code for BSpline basis sets.
module;
#include <bspline/Core.h>
#include "radial/BSpline/BFGrouper.H"

export module qchem.Basisset.Atom.radial.BSpline.IE;
import qchem.BasisSet.Atom.BFGrouper;
import qchem.Basisset.Atom.radial.BSpline.IEC;
import qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Imp.HeapDB;
import qchem.BasisSet.Imp.Cache4;
import qchem.BasisSet.Imp.HeapDB;
import qchem.BasisSet.Imp.IEClient;
import qchem.BasisSet.Integrals;
import qchem.BasisSet.ERI4;
import qchem.BasisSet;
import qchem.HF_IBS;
import qchem.DHF_IBS;
import qchem.Fit_IBS;

export namespace BSpline
{
// K is the spline order
template <class T, size_t K> class Primative_Overlap
{
    typedef bspline::Spline<T, K> spline_t;
public:
    virtual double Overlap(const spline_t& a , const spline_t& b,size_t l_total) const=0;
};
template <class T, size_t K> class Primative_Grad2
{
    typedef bspline::Spline<T, K> spline_t;
public:
    virtual double Grad2(const spline_t& a , const spline_t& b,size_t la, size_t lb) const=0;
};
template <class T, size_t K> class Primative_Inv_r1
{
    typedef bspline::Spline<T, K> spline_t;
public:
    virtual double Inv_r1(const spline_t& a , const spline_t& b,size_t l_total) const=0;
};
template <class T, size_t K> class Primative_Inv_r2
{
    typedef bspline::Spline<T, K> spline_t;
public:
    virtual double Inv_r2(const spline_t& a , const spline_t& b,size_t l_total) const=0;
};
template <class T, size_t K> class Primative_Repulsion
{
    typedef bspline::Spline<T, K> spline_t;
public:
    virtual double Repulsion(const spline_t& a , const spline_t& b,size_t la, size_t lc) const=0;
};
template <class T, size_t K> class Primative_Charge
{
    typedef bspline::Spline<T, K> spline_t;
public:
    virtual double Charge   (const spline_t& a, size_t l) const=0;
};

template <class T, size_t K> class IE_Overlap
: public virtual Primative_Overlap<T,K>
, public DB_Overlap<T>
{
protected:
    typedef bspline::Spline<T, K> spline_t;
    using Primative_Overlap<T,K>::Overlap;
    virtual typename Integrals_Base<T>::SMat MakeOverlap() const;
    IE_Overlap(const DB_cache<T>* db) : DB_Overlap<T>(db) {};
};
template <class T, size_t K> class IE_Kinetic
: public virtual Primative_Grad2<T,K>
, public virtual Primative_Inv_r2<T,K>
, public DB_Kinetic<T>
{
protected:
    using Primative_Grad2 <T,K>::Grad2;
    using Primative_Inv_r2<T,K>::Inv_r2;
    virtual typename Integrals_Base<T>::SMat MakeKinetic() const;
    IE_Kinetic(const DB_cache<T>* db) : DB_Kinetic<T>(db) {};
};
template <class T, size_t K> class IE_Inv_r1
: public virtual Primative_Inv_r1<T,K>
, public DB_Nuclear<T>
{
protected:
    using Primative_Inv_r1<T,K>::Inv_r1;
    virtual typename Integrals_Base<T>::SMat MakeNuclear(const Cluster* cl) const;
    IE_Inv_r1(const DB_cache<T>* db) : DB_Nuclear<T>(db) {};
};
template <class T, size_t K> class IE_XGrad2
: public virtual Primative_Grad2<T,K>
, public DB_XKinetic<T>
{
protected:
    using Primative_Grad2<T,K>::Grad2;
    virtual typename Integrals_Base<T>::Mat MakeKinetic(const Orbital_RKBS_IBS<T>* rkbs) const;
    IE_XGrad2(const DB_cache<T>* db) : DB_XKinetic<T>(db) {};
};

template <class T, size_t K> class IE_BS_2E 
    : public virtual ::AtomIE_BS_2E_Angular
    , public virtual Cache4
    , public DB_BS_2E<T>
    , public BFGrouper<K>
{
    typedef typename ::AtomIE_BS_2E_Angular::RVec RVec;
public:
    virtual ERI4 MakeDirect  (const ::IrrepIEClient* a, const ::IrrepIEClient* c) const;
    virtual ERI4 MakeExchange(const ::IrrepIEClient* a, const ::IrrepIEClient* c) const;

    // Cach4 functions
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc) const=0;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc) const=0;
protected:
    virtual void Append(const ::IrrepIEClient*);
};

template <class T, size_t K> class IE_DFT 
: public virtual Primative_Overlap<T,K>
, public virtual Primative_Repulsion<T,K>
, public DB_DFT<T>
{
    typedef Integrals_Base<T> Base;
    typedef typename Base::SMat SMat;
    typedef typename Base::ERI3 ERI3;
    typedef bspline::Spline<T, K> spline_t;
protected:
    IE_DFT(const DB_cache<T>* db) : DB_DFT<T>(db) {};
    
    virtual ERI3 MakeOverlap3C  (const Fit_IBS& c) const;
    virtual ERI3 MakeRepulsion3C(const Fit_IBS& c) const;
private:
    typedef typename BSpline::IrrepIEClient<K>::bf_tuple bf_tuple;
    SMat MakeOverlap  (const bf_tuple& c) const; //ab loops
    SMat MakeRepulsion(const bf_tuple& c) const; //ab loops
};
template <class T, size_t K> class IE_RKBL 
    : public IE_Overlap<T,K>
    , public IE_XGrad2 <T,K>
    , public IE_Inv_r1<T,K>
{
protected:
    IE_RKBL(const DB_cache<T>* db) : IE_Overlap<T,K>(db),IE_XGrad2<T,K>(db),IE_Inv_r1<T,K>(db) {};

};
template <class T, size_t K> class IE_RKBS 
: public IE_Kinetic  <T,K>
, public IE_Inv_r1<T,K>
{
protected:
    IE_RKBS(const DB_cache<T>* db) : IE_Kinetic<T,K>(db), IE_Inv_r1<T,K>(db) {};
};
template <size_t K> class IE_Fit 
: public virtual Primative_Repulsion<double,K>
, public virtual Primative_Charge<double,K>
, public DB_Fit
{
    protected:
    IE_Fit(const DB_cache<double>* db) : DB_Fit(db) {};

    virtual Vec  MakeCharge() const;
    virtual SMat MakeRepulsion() const;
    virtual  Mat MakeRepulsion(const Fit_IBS&) const;
private:
    // Derived classes must provide the actual integral calculations.
    using DB_Fit::Charge; //un hide
    using DB_Fit::Repulsion; //un hide
    using Primative_Repulsion<double,K>::Repulsion;
    using Primative_Charge   <double,K>::Charge;
};

} //namespace
