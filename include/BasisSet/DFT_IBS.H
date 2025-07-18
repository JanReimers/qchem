// File: DFT_IBS.H  Interface for a Density Functional Theory (DFT) Orbital Irrep Basis Set.
#ifndef _DFT_IBS_H_
#define _DFT_IBS_H_

#include <BasisSet/Irrep_BS.H>
#include <BasisSet/Integrals.H>

//! \brief Interface for 3-center integrals used in DFT calculations.
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


template <class T> class TOrbital_DFT_IBS
    : public virtual TOrbital_IBS<T>
    , public virtual Integrals_DFT<T> //DFT integrals
    
{
    typedef typename Integrals_Base<T>::SMat SMat;
    typedef typename Integrals_Base<T>::fbs_t fbs_t;
    typedef Vector<T> Vec;
public:
    virtual fbs_t*    CreateCDFitBasisSet(const BasisSet*,const Cluster*) const=0;
    virtual fbs_t*    CreateVxcFitBasisSet(const BasisSet*,const Cluster*) const=0;
    using Integrals_DFT<T>::Overlap3C; //Unhide
    using Integrals_DFT<T>::Repulsion3C; //Unhide
    virtual Vec Overlap3C  (const SMat& Dcd, const fbs_t* ff) const=0;
    virtual Vec Repulsion3C(const SMat& Dcd, const fbs_t* ff) const=0;
};

#endif //_DFT_IBS_H_
