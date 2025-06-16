// File: DFT_IE.H Interface for 3-center integrals used in DFT calculations.
#ifndef _DFT_IE_H_
#define _DFT_IE_H_

#include <BasisSet/IntegralEngine.H>

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

#endif //_DFT_IE_H_
