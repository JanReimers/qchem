// File: BasisSet/Orbital_DFT_IBS.C  Interface for a Density Functional Theory (DFT) Orbital Irrep Basis Set.
module;
export module qchem.BasisSet.Orbital_DFT_IBS;
export import qchem.BasisSet.IrrepBasisSet;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Fit_IBS;
export import qchem.BasisSet.Internal.ERI3;

export namespace qchem::BasisSet
{

template <class T> class Orbital_DFT_IBS
    : public virtual Orbital_1E_IBS<T>
    , public virtual IrrepBasisSet_IDs   // Name / BasisSetID identity face
{
public:
    //! 3 centre overlap used for DFT \f$ \left\langle ab\left|1\right|c\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right)f_{c}\left(\vec{r}\right) \f$.  The fit \a c is an overlap-metric (scalar-function) aux basis.
    virtual const ERI3<T>& Overlap3C  (const rFIT_SF_ABS& c) const;
    //! 3 centre repulsion used for DFT \f$\left\langle a\left(1\right)b\left(1\right)\left|\frac{1}{r_{12}}\right|c\left(2\right)\right\rangle =\int d^{3}\vec{r}_{1}\:d^{3}\vec{r}_{2}\:g_{a}\left(\vec{r}_{1}\right)g_{b}\left(\vec{r}_{1}\right)\frac{1}{r_{12}}f_{c}\left(\vec{r}_{2}\right) \f$.  The fit \a c is a Coulomb-metric (charge-density) aux basis.
    virtual const ERI3<T>& Repulsion3C(const rFIT_CD_ABS& c) const;

    virtual FIT_CD_ABS<T>* CreateCDFitBasisSet (const Structure*, const qcMesh::MeshParams&) const=0;
    virtual FIT_SF_ABS<T>* CreateVxcFitBasisSet(const Structure*, const qcMesh::MeshParams&) const=0;
    // NB: there is deliberately NO Repulsion3C(D,c)/Overlap3C(D,c) here -- the density-matrix contraction
    // <rho|c> = Sum_ab D_ab <ab|c> is a DENSITY operation (it lives in the charge density, which owns D),
    // NOT a basis one.  The basis exposes only the D-free integral tensors above; D never enters qcBasisSet.
protected:
    virtual ERI3<T> MakeOverlap3C  (const rFIT_SF_ABS& c) const=0;
    virtual ERI3<T> MakeRepulsion3C(const rFIT_CD_ABS& c) const=0;
};

typedef Orbital_DFT_IBS<double> Real_DFT_OIBS;
} //namespace