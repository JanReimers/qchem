// File: BasisSet/Orbital_DFT_IBS.C  Interface for a Density Functional Theory (DFT) Orbital Irrep Basis Set.
module;
export module qchem.BasisSet.Orbital_DFT_IBS;
export import qchem.BasisSet.IrrepBasisSet;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Fit_IBS;
export import qchem.BasisSet.Internal.ERI3;

export namespace qchem::BasisSet
{

//! \brief What an IBS must supply to support DFT -- on TWO independent axes (V1.1).
//!
//! \tparam T    the ORBITAL representation: real (molecular / a TRIM Bloch block) or dcmplx (general k).
//! \tparam TFit the FIT-BASIS representation, defaulting to \a T.
//!
//! They are independent because the fit basis is a property of the RUN, not of a block: every call site
//! builds it once from the whole \c tBasisSet at Hamiltonian construction, and all irrep/k blocks share
//! that one object.  So under doc/RealComplexPlan.md, where \c BlochQN::IsReal() turns TRIM k-blocks REAL
//! while the rest stay complex, those real blocks still share the run's COMPLEX G-space fit basis.  That
//! case -- \c Orbital_DFT_IBS<double,dcmplx> -- had no spelling at all before this parameter.
//!
//! Band_FT_IBS's own comment already knew ("never assuming orbital==fit"); it simply had no second
//! parameter to say it with, which is why it exists as a parallel dcmplx-hardwired class instead of an
//! INSTANTIATION of this one.  Before this the class even contradicted itself: \c CreateCDFitBasisSet
//! returned \c FIT_CD_ABS<T>* while \c Repulsion3C accepted only \c FIT_CD_ABS<double>, so
//! \c Orbital_DFT_IBS<dcmplx> could not consume the fit basis it created.
//!
//! Three of the four combinations are meaningful (user 2026-08-16): \c <double,double> molecular,
//! \c <dcmplx,dcmplx> GPW, \c <double,dcmplx> a real block on a periodic run.  \c <dcmplx,double> is not
//! used -- complex orbitals do not get a real fit basis -- and is simply never instantiated.
template <class T, class TFit = T> class Orbital_DFT_IBS
    : public virtual Orbital_1E_IBS<T>
    , public virtual IrrepBasisSet_IDs   // Name / BasisSetID identity face
{
public:
    //! 3 centre overlap used for DFT \f$ \left\langle ab\left|1\right|c\right\rangle \f$.  The fit \a c is an
    //! overlap-metric (scalar-function) aux basis, in the FIT representation.
    virtual const ERI3<TFit>& Overlap3C  (const FIT_SF_ABS<TFit>& c) const;
    //! 3 centre repulsion used for DFT.  The fit \a c is a Coulomb-metric (charge-density) aux basis.
    virtual const ERI3<TFit>& Repulsion3C(const FIT_CD_ABS<TFit>& c) const;

    virtual FIT_CD_ABS<TFit>* CreateCDFitBasisSet (const Structure*, const qcMesh::MeshParams&) const=0;
    virtual FIT_SF_ABS<TFit>* CreateVxcFitBasisSet(const Structure*, const qcMesh::MeshParams&) const=0;
    // NB: there is deliberately NO Repulsion3C(D,c)/Overlap3C(D,c) here -- the density-matrix contraction
    // <rho|c> = Sum_ab D_ab <ab|c> is a DENSITY operation (it lives in the charge density, which owns D),
    // NOT a basis one.  The basis exposes only the D-free integral tensors above; D never enters qcBasisSet.
    //
    // V1.1(ii), the resolvable half: the tensor element type follows TFit, NOT T.  A no-op today (every
    // instantiation has T==TFit), but it DOCUMENTS which axis decides -- which was the ambiguous part.  The
    // integrand is g_a g_b f_c, so complex FIT functions make <ab|c> complex even with real orbitals: TFit
    // is right for all three meaningful combinations and wrong only for <dcmplx,double>, which is never
    // instantiated.
    //
    // STILL OPEN, and it is the rest of (ii): the periodic path returns G_ERI3, not ERI3.  Reconciling the
    // two is what decides whether Band_FT_IBS can BE Orbital_DFT_IBS<dcmplx,dcmplx> or must stay parallel --
    // the head of the type flow, so it belongs with doc/RealComplexPlan.md rather than here.
protected:
    virtual ERI3<TFit> MakeOverlap3C  (const FIT_SF_ABS<TFit>& c) const=0;
    virtual ERI3<TFit> MakeRepulsion3C(const FIT_CD_ABS<TFit>& c) const=0;
};

typedef Orbital_DFT_IBS<double> Real_DFT_OIBS;
} //namespace