// File: BSpline/IE.C Common IE code for BSpline basis sets.
module;
#include <bspline/Core.h>
#include <memory>
export module qchem.Basisset.Atom.radial.BSpline.IE;
import qchem.Basisset.Atom.radial.BSpline.BFGrouper;
import qchem.Basisset.Atom.radial.BSpline.IEC;

import qchem.BasisSet.Atom.Internal.BFGrouper;
import qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Internal.Cache4;
import qchem.BasisSet.Internal.HeapDB;
import qchem.BasisSet.Internal.IEClient;

import qchem.BasisSet.Internal.ERI4;
import qchem.BasisSet;
import qchem.Orbital_HF_IBS;
import qchem.Orbital_DHF_IBS;
import qchem.Orbital_DFT_IBS;
import qchem.Fit_IBS;

export namespace BSpline
{

template <size_t K> class IE_Primatives
 {
    typedef bspline::Spline<double, K> spline_t;
public:
    virtual double Overlap  (const spline_t& a , const spline_t& b,size_t l_total     ) const=0;
    virtual double Grad2    (const spline_t& a , const spline_t& b,size_t la,size_t lb) const=0;
    virtual double Inv_r1   (const spline_t& a , const spline_t& b,size_t l_total     ) const=0; //! <a|1/r^1|b>
    virtual double Inv_r2   (const spline_t& a , const spline_t& b,size_t l_total     ) const=0; //! <a|1/r^2|b>
    virtual double Repulsion(const spline_t& a , const spline_t& b,size_t la,size_t lc) const=0;
    virtual double Charge   (const spline_t& a ,                   size_t l           ) const=0;
};

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
    virtual SMatrix<T> MakeOverlap() const;
    IE_Overlap(const DB_cache<T>* db,const ::BSpline::IE_Primatives<K>* _pie) : DB_Overlap<T>(db), pie(_pie) {};
private:
    const IE_Primatives<K>* pie;
};
template <class T, size_t K> class IE_Kinetic
: public virtual Primative_Grad2<T,K>
, public virtual Primative_Inv_r2<T,K>
, public DB_Kinetic<T>
{
protected:
    using Primative_Grad2 <T,K>::Grad2;
    using Primative_Inv_r2<T,K>::Inv_r2;
    virtual SMatrix<T> MakeKinetic() const;
    IE_Kinetic(const DB_cache<T>* db,const ::BSpline::IE_Primatives<K>* _pie) : DB_Kinetic<T>(db), pie(_pie)  {};
private:
    const IE_Primatives<K>* pie;
};
template <class T, size_t K> class IE_Inv_r1
: public virtual Primative_Inv_r1<T,K>
, public DB_Nuclear<T>
{
protected:
    using Primative_Inv_r1<T,K>::Inv_r1;
    virtual SMatrix<T> MakeNuclear(const Cluster* cl) const;
    IE_Inv_r1(const DB_cache<T>* db,const ::BSpline::IE_Primatives<K>* _pie) : DB_Nuclear<T>(db), pie(_pie)  {};
private:
    const IE_Primatives<K>* pie;
};

template <class T, size_t K> class IE_BS_2E 
    : public virtual Cache4
    , public DB_BS_2E<T>
    , public BFGrouper<K>
{
public:
    IE_BS_2E(AtomIE_BS_2E_Angular* a) : itsAngular(a) {};
    virtual ERI4 MakeDirect  (const ::IrrepIEClient* a, const ::IrrepIEClient* c) const;
    virtual ERI4 MakeExchange(const ::IrrepIEClient* a, const ::IrrepIEClient* c) const;

    // Cach4 functions
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc) const=0;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc) const=0;
protected:
    virtual void Append(const ::IrrepIEClient*);
private: 
    std::unique_ptr<AtomIE_BS_2E_Angular> itsAngular;
};

template <class T, size_t K> class IE_DFT 
: public virtual Primative_Overlap<T,K>
, public virtual Primative_Repulsion<T,K>
, public DB_DFT<T>
{
    typedef bspline::Spline<T, K> spline_t;
protected:
    IE_DFT(const DB_cache<T>* db,const ::BSpline::IE_Primatives<K>* _pie) : DB_DFT<T>(db), pie(_pie)  {};
    
    virtual ERI3<T> MakeOverlap3C  (const Fit_IBS& c) const;
    virtual ERI3<T> MakeRepulsion3C(const Fit_IBS& c) const;
private:
    typedef typename BSpline::IrrepIEClient<K>::bf_tuple bf_tuple;
    SMatrix<T> MakeOverlap  (const bf_tuple& c) const; //ab loops
    SMatrix<T> MakeRepulsion(const bf_tuple& c) const; //ab loops
private:
    const IE_Primatives<K>* pie;
};

template <size_t K> class IE_Fit 
: public virtual Primative_Repulsion<double,K>
, public virtual Primative_Charge<double,K>
, public DB_Fit
{
    protected:
    IE_Fit(const DB_cache<double>* db,const ::BSpline::IE_Primatives<K>* _pie) : DB_Fit(db), pie(_pie)  {};

    virtual  Vector<double> MakeCharge   () const;
    virtual SMatrix<double> MakeRepulsion() const;
    virtual  Matrix<double> MakeRepulsion(const Fit_IBS&) const;
private:
    // Derived classes must provide the actual integral calculations.
    using DB_Fit::Charge; //un hide
    using DB_Fit::Repulsion; //un hide
    using Primative_Repulsion<double,K>::Repulsion;
    using Primative_Charge   <double,K>::Charge;
private:
    const IE_Primatives<K>* pie;
};

} //namespace
