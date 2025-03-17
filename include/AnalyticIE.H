// File: AnalyticIE.H  Abtract integral engine interface.
#ifndef _AnalyticIE_H_
#define _AnalyticIE_H_


#include "types.H"
#include "Imp/Misc/UniqueID/UniqueID.H"
#include <vector>
#include <map>
#include <iosfwd>

class IrrepIEClient;
//--------------------------------------------------------------------------------
//
//!  \brief This class is abstract interface used for calculating integrals over the basis functions.
//!
//!  The method of integral evaluation is of course strongly dependant on the
//!  precise details of basis function or basis set.  Hence for each concrete type
//!  of basis function there is a corresponding basis set and integral engine.
//!  The integral engine does not see the basis function/set objects directly. Instead it
//!  works with the IrrepIEClient interface which simply supply list of
//!  exponents, polarizations and normalization constants for the particular irrep basis set.
//!  All functions except MakeNormalization return normalized integrals.
//!
//
// Interface for non-relativistic 1 electron integrals.  
// The calls return matrix refrences which implies they are buffered behind the scenes.
//
template <class T> class Integrals_Base
{
public:
    typedef TOrbital_IBS<T> obs_t;
    typedef TOrbital_HF_IBS<T> hf_obs_t;
    typedef TOrbital_DFT_IBS<T> dft_obs_t;
    typedef Fit_IBS fbs_t;

    typedef  Vector<T>  Vec;
    typedef  Matrix<T>  Mat;
    typedef SMatrix<T> SMat;
    typedef std::vector<SMat> ERI3;
    typedef const  Vec&  Vec_ref;
    typedef const  Mat&  Mat_ref;
    typedef const SMat& SMat_ref;
};

template <class T> class Integrals_Overlap : public virtual Integrals_Base<T>
{
public:
    //! Single basis set Overlap \f$ \left\langle a\left|1\right|b\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right) \f$ 
    virtual typename Integrals_Base<T>::SMat_ref Overlap() const=0;
};

template <class T> class Integrals_Kinetic : public virtual Integrals_Base<T>
{
public:
    //! Kinetic energy \f$ \left\langle a\left|-\frac{1}{2}\nabla^{2}\right|b\right\rangle =-\frac{1}{2}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\nabla^{2}g_{b}\left(\vec{r}\right)\f$
    virtual typename Integrals_Base<T>::SMat_ref Kinetic() const=0;
};

template <class T> class Integrals_Nuclear : public virtual Integrals_Base<T>
{
public:
    //! Nuclear attraction \f$ \sum_{i}\left\langle a\left|\frac{-Z_{i}}{\left|\vec{r}-\vec{R}_{c}\right|}\right|b\right\rangle =-\sum_{i}Z_{i}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\frac{1}{\left|\vec{r}-\vec{R}_{c}\right|}g_{b}\left(\vec{r}\right)\f$
    virtual typename Integrals_Base<T>::SMat_ref Nuclear(const Cluster*) const=0;   
};



template <class T> class Integrals_1E 
    : public virtual Integrals_Overlap<T>
    , public virtual Integrals_Kinetic<T>
    , public virtual Integrals_Nuclear<T>
{
};

template <class T> class Integrals_DFT : public virtual Integrals_Base<T>
{
    typedef Integrals_Base<T> Base;
    typedef typename Base::fbs_t fbs_t;
    typedef typename Base::SMat SMat;
    typedef typename Base::ERI3 ERI3;
public:
    //! 3 centre overlap used for DFT \f$ \left\langle ab\left|1\right|c\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right)f_{c}\left(\vec{r}\right) \f$
    virtual const ERI3& Overlap3C  (const fbs_t& c) const=0; 
    //! 3 centre repulsion used for DFT \f$\left\langle a\left(1\right)b\left(1\right)\left|\frac{1}{r_{12}}\right|c\left(2\right)\right\rangle =\int d^{3}\vec{r}_{1}\:d^{3}\vec{r}_{2}\:g_{a}\left(\vec{r}_{1}\right)g_{b}\left(\vec{r}_{1}\right)\frac{1}{r_{12}}f_{c}\left(\vec{r}_{2}\right) \f$
    virtual const ERI3& Repulsion3C(const fbs_t& c) const=0; 

};

template <class T> class Integrals_2E : public virtual Integrals_Base<T>
{
    typedef Integrals_Base<T> Base;
    typedef typename Base::obs_t obs_t;
public:
    virtual ERI4 Direct  (const obs_t& c) const=0; //! <ab|1/r_12|cd>, this=a.
    virtual ERI4 Exchange(const obs_t& b) const=0; //! <ac|1/r_12|bd>, this=a.

};

template <class T> class Integrals_BS_2E
{
public:
    typedef UniqueID::IDtype IDType;   
    virtual ERI4 Direct  (IDType a,IDType c) const=0;
    virtual ERI4 Exchange(IDType a,IDType b) const=0;
};

template <class T> class Integrals_RKB
: public virtual Integrals_1E<T>
{
public:
    virtual  typename Integrals_Base<T>::SMat_ref RestMass() const=0;
};

template <class T> class Integrals_RKBL 
{
public:
    typedef SMatrix<T> SMat;
    typedef  Matrix<T>  Mat;
    typedef const SMat& SMat_ref;
    typedef const  Mat&  Mat_ref;
    virtual SMat_ref Overlap() const=0;
    //! L/S cross Kinetic energy \f$ \left\langle a\left|-\frac{1}{2}\nabla^{2}\right|b\right\rangle =-\frac{1}{2}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\nabla^{2}g_{b}\left(\vec{r}\right)\f$
    virtual  Mat_ref Kinetic(const Orbital_RKBS_IBS<T>* rkbs) const=0;
    virtual SMat_ref Nuclear(const Cluster*) const=0;   
};

template <class T> class Integrals_RKBS 
{
public:
    typedef SMatrix<T> SMat;
    typedef const SMat& SMat_ref;
    virtual SMat_ref Overlap() const=0;
    virtual SMat_ref Nuclear(const Cluster*) const=0;   
    //! Rest mass \f$ \left\langle a\left|\left(\beta-\alpha\right)c^{2}\right|b\right\rangle =\left(\beta-\alpha\right)c^{2}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right) \f$
    virtual SMat_ref RestMass() const=0;
};

class FitIntegrals_1E  
    : public virtual Integrals_Base<double>
    , public virtual Integrals_Overlap<double>
{
public:
    //! Single basis set Overlap \f$ \left\langle a\left|1\right|b\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right) \f$ 
    using Integrals_Overlap<double>::Overlap;
    virtual Vec_ref  Charge   () const=0;   
    virtual SMat_ref Repulsion() const=0;
    virtual  Mat_ref Repulsion(const fbs_t&) const=0;
    virtual SMat_ref InvOverlap(const LAParams&) const=0;
    virtual SMat_ref InvRepulsion(const LAParams&) const=0;
    // Pure numerial versions
    virtual  Vec_ref Norm   (const Mesh*        ) const=0; //Numerical .
    virtual  Vec_ref Charge (const Mesh*        ) const=0; //Numerical .
    virtual  Mat_ref Overlap(const Mesh*,const fbs_t& b) const=0; //Numerical X overlap.
//
//  These are used for charge and Vxc fitting.  They change with iterations
//  So they MUST not be cached.
//
    typedef ScalarFunction<double> Sf;
    virtual const Vec Overlap    (const Mesh*,const Sf&) const=0; //Numerical  
    virtual const Vec Repulsion  (const Mesh*,const Sf&) const=0; //Numerical 

    protected:
    virtual  Vec MakeNorm   (const Mesh*        ) const=0; //Numerical .
    virtual  Vec MakeCharge (const Mesh*        ) const=0; //Numerical .
    virtual  Mat MakeOverlap(const Mesh*,const fbs_t& b) const=0; //Numerical X overlap.
};

#endif //_AnalyticIE_H_
