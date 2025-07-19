// File: IntegralEngine.H  Abtract integral engine interfaces.
#ifndef _IntegralEngine_H_
#define _IntegralEngine_H_


#include <vector>
#include <BasisSet/fwd.H>

import Common.UniqueID; //Need IDType
import qchem.Cluster;
import oml;

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


//! \brief Define all commonly used typedefs in one place.
template <class T> class Integrals_Base
{
public:
    // typedef TOrbital_IBS<T> obs_t;
    // typedef TOrbital_HF_IBS<T> hf_obs_t;
    typedef TOrbital_DFT_IBS<T> dft_obs_t;
    // typedef Fit_IBS fbs_t;

    typedef  Vector<T>  Vec;
    typedef  Matrix<T>  Mat;
    typedef SMatrix<T> SMat;
    typedef std::vector<SMat> ERI3;
    typedef const  Vec&  Vec_ref;
    typedef const  Mat&  Mat_ref;
    typedef const SMat& SMat_ref;
};

//! \brief Interface for overlap integrals.
template <class T> class Integrals_Overlap : public virtual Integrals_Base<T>
{
public:
    //! Single basis set Overlap \f$ \left\langle a\left|1\right|b\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right) \f$ 
    virtual typename Integrals_Base<T>::SMat_ref Overlap() const=0;
};
//! \brief Interface for Laplacian integrals using in kinetic energy calculations.
template <class T> class Integrals_Kinetic : public virtual Integrals_Base<T>
{
public:
    //! Grad^2 \f$ \left\langle a\left|-\frac{1}{2}\nabla^{2}\right|b\right\rangle =-\frac{1}{2}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\nabla^{2}g_{b}\left(\vec{r}\right)\f$
    virtual typename Integrals_Base<T>::SMat_ref Kinetic() const=0;
};
//! \brief Interface for electron-nucleus attraction integrals.
template <class T> class Integrals_Nuclear : public virtual Integrals_Base<T>
{
public:
    //! Nuclear attraction \f$ \sum_{i}\left\langle a\left|\frac{-Z_{i}}{\left|\vec{r}-\vec{R}_{c}\right|}\right|b\right\rangle =-\sum_{i}Z_{i}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\frac{1}{\left|\vec{r}-\vec{R}_{c}\right|}g_{b}\left(\vec{r}\right)\f$
    virtual typename Integrals_Base<T>::SMat_ref Nuclear(const Cluster*) const=0;   
};
//! \brief Interface for rest mass matrix used to shift relativistic energies.
template <class T> class Integrals_RestMass : public virtual Integrals_Base<T>
{
public:
    //! Rest mass \f$ \left\langle a\left|\left(\beta-\alpha\right)c^{2}\right|b\right\rangle =\left(\beta-\alpha\right)c^{2}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right) \f$
    virtual typename Integrals_Base<T>::SMat_ref RestMass() const=0;   
};
//! \brief Interface for L-S cross kinetic matrix used in relativistic calculations.
template <class T> class Integrals_XKinetic : public virtual Integrals_Base<T>
{
public:
    //! L/S cross Grad^2 \f$ \left\langle a\left|-\frac{1}{2}\nabla^{2}\right|b\right\rangle =-\frac{1}{2}\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)\nabla^{2}g_{b}\left(\vec{r}\right)\f$
    virtual typename Integrals_Base<T>::Mat_ref Kinetic(const Orbital_RKBS_IBS<T>* rkbs) const=0;
};


#endif //_IntegralEngine_H_
