// File: BasisSet/Orbital_DFT_IBS.C  Interface for a Density Functional Theory (DFT) Orbital Irrep Basis Set.
module;
export module qchem.BasisSet.Orbital_DFT_IBS;
export import qchem.BasisSet.IrrepBasisSet;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Fit_IBS;
export import qchem.BasisSet.Internal.ERI3;

export namespace BasisSet
{

template <class T> class Orbital_DFT_IBS
    : public virtual Orbital_1E_IBS<T>
    , public virtual IrrepBasisSet_IDs //avoid using statements for RadialID,AngularID
{
public:
    //! 3 centre overlap used for DFT \f$ \left\langle ab\left|1\right|c\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right)f_{c}\left(\vec{r}\right) \f$
    virtual const ERI3<T>& Overlap3C  (const Fit_IBS& c) const;
    //! 3 centre repulsion used for DFT \f$\left\langle a\left(1\right)b\left(1\right)\left|\frac{1}{r_{12}}\right|c\left(2\right)\right\rangle =\int d^{3}\vec{r}_{1}\:d^{3}\vec{r}_{2}\:g_{a}\left(\vec{r}_{1}\right)g_{b}\left(\vec{r}_{1}\right)\frac{1}{r_{12}}f_{c}\left(\vec{r}_{2}\right) \f$
    virtual const ERI3<T>& Repulsion3C(const Fit_IBS& c) const;

    virtual Fit_IBS* CreateCDFitBasisSet (const Cluster*) const=0;
    virtual Fit_IBS* CreateVxcFitBasisSet(const Cluster*) const=0;
    virtual vec_t<T> Overlap3C  (const smat_t<T>& Dcd, const Fit_IBS* c) const;
    virtual vec_t<T> Repulsion3C(const smat_t<T>& Dcd, const Fit_IBS* c) const;
protected:
    virtual ERI3<T> MakeOverlap3C  (const Fit_IBS& c) const=0;
    virtual ERI3<T> MakeRepulsion3C(const Fit_IBS& c) const=0;
};

typedef Orbital_DFT_IBS<double> Real_DFT_OIBS;
} //namespace