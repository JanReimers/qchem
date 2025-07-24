// File: DFT_IBS.H  Interface for a Density Functional Theory (DFT) Orbital Irrep Basis Set.
module;
#include <vector>
export module qchem.DFT_IBS;
export import qchem.BasisSet;
export import qchem.Irrep_BS;
import qchem.BasisSet.Internal.Integrals;

export template <class T> using ERI3=std::vector<SMatrix<T>>;

//! \brief Interface for 3-center integrals used in DFT calculations.
export template <class T> class Integrals_DFT 
{
public:
    //! 3 centre overlap used for DFT \f$ \left\langle ab\left|1\right|c\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right)f_{c}\left(\vec{r}\right) \f$
    virtual const ERI3<T>& Overlap3C  (const Fit_IBS& c) const=0; 
    //! 3 centre repulsion used for DFT \f$\left\langle a\left(1\right)b\left(1\right)\left|\frac{1}{r_{12}}\right|c\left(2\right)\right\rangle =\int d^{3}\vec{r}_{1}\:d^{3}\vec{r}_{2}\:g_{a}\left(\vec{r}_{1}\right)g_{b}\left(\vec{r}_{1}\right)\frac{1}{r_{12}}f_{c}\left(\vec{r}_{2}\right) \f$
    virtual const ERI3<T>& Repulsion3C(const Fit_IBS& c) const=0; 

};


export template <class T> class TOrbital_DFT_IBS
    : public virtual TOrbital_IBS<T>
    , public virtual Integrals_DFT<T> //DFT integrals
    
{
public:
    virtual Fit_IBS*    CreateCDFitBasisSet(const BasisSet*,const Cluster*) const=0;
    virtual Fit_IBS*    CreateVxcFitBasisSet(const BasisSet*,const Cluster*) const=0;
    using Integrals_DFT<T>::Overlap3C; //Unhide
    using Integrals_DFT<T>::Repulsion3C; //Unhide
    virtual Vector<T> Overlap3C  (const SMatrix<T>& Dcd, const Fit_IBS* ff) const=0;
    virtual Vector<T> Repulsion3C(const SMatrix<T>& Dcd, const Fit_IBS* ff) const=0;
};

