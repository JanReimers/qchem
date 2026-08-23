// File: BasisSet.C Quantum Chemistry basis set expressed as a sequence of Irrep basis sets.
module;
#include <vector>
#include <memory>

export module qchem.BasisSet;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Fit_IBS;
export import qchem.Structure;
export import qchem.Symmetry;
export import qchem.ElectronConfiguration;
export import qchem.Matrix3D;   // Matrix3D -- the crystal point-group ops a periodic basis exposes (IBZ symmetrization)
export import qchem.Symmetry.Lattice_3D.SpaceGroup;   // ReciprocalOp {U|τ} -- the density-symmetrization ops (glide phase)

import qchem.Iterators;
export import qchem.Streamable;

export namespace qchem::BasisSet
{
typedef std::vector<Irrep> irrepv_t;

//----------------------------------------------------------------------------
//
//  Interface for a BasisSet which is assumed to a list of Irrep Basis Sets.
//
template <class T> class tBasisSet
    : public virtual Streamable
{
public:
    typedef Orbital_1E_IBS<T> obs_t;

    virtual ~tBasisSet() {};
    virtual size_t   GetNumFunctions() const=0;
    virtual irrepv_t GetIrreps(const Spin& ms) const=0;

    virtual FIT_CD_ABS<T>* CreateCDFitBasisSet(const Structure* cl, const qcMesh::MeshParams&) const;
    //! \brief The \f$v_{xc}\f$ fit basis for this run: \a fit picks the REPRESENTATION and \a mp the
    //! POINTS, and ONE object comes back carrying both -- because a fit basis is a family of weight vectors
    //! over shared points, so its mesh is constitutive of it, not a sibling return value.
    //!
    //! THIS is the level that chooses between representations (a Bloch BLOCK's own factory builds only the
    //! lineage's fitted basis): \c VxcFit::Delta returns the \f$\delta\f$ basis over the run's XC
    //! quadrature (\c CreateXCQuadrature -- mesh + fold + Shubnikov tags), anything else the lineage's
    //! fitted one.  A molecular caller passes \c Auto and gets its Gaussian auxiliary basis, as always.
    //! \param partition OPTIONAL out: the real-space quadrature the returned basis was built over, when it
    //! has one (\c VxcFit::Delta; null otherwise).  INJECTION, not a getter (user, 2026-08-23): this
    //! factory CREATES the mesh, so handing it to a second collaborator is what a creator does -- whereas
    //! a \c Mesh() accessor on the fit basis would hand out its own private state, which is exactly the
    //! escape the 2026-08-22 arc spent itself closing.  The point is that an atom-partitioned quadrature
    //! is a GENERALLY useful object with nothing to do with \f$v_{xc}\f$ or \f$\rho\f$ fitting, and
    //! today this is the only place in a run where one is built; a consumer of atomic observables gets the
    //! SAME object the basis integrates on rather than rebuilding it.
    virtual FIT_SF_ABS<T>* CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams&,
                                                VxcFit fit=VxcFit::Auto,
                                                std::shared_ptr<const qcMesh::Mesh>* partition=nullptr) const;
    //! The DELTA-fit sibling of \c CreateVxcFitBasisSet: the finished real-space XC quadrature (mesh +
    //! symmetry fold), assembled by the basis -- which owns the cell and the §3-imposed ops -- so the
    //! Hamiltonian does no mesh work.  Default/plain path: the Structure's integration mesh, no fold.
    virtual FitQuadrature   CreateXCQuadrature (const Structure* cl, const qcMesh::MeshParams&) const;

    //! The crystal RECIPROCAL point group as \f${U|\tau}\f$ ops for IBZ density symmetrization.  Default {} =
    //! trivial {E} = no-op -- molecules / Γ / unfolded bases.  A periodic GPW basis returns the ops WHEN it
    //! folds the mesh (imposeSymmetry), so the composite density ctor-injects them and the reduced density is star-
    //! averaged with the glide phase \f$e^{+2\pi i(Um)\cdot\tau}\f$ (doc/GPWPlan1.md items 3 + 5).  It is a basis
    //! property (the basis computes the space group), so the density/WF read it here rather than via a setter.
    virtual std::vector<Symmetry::Lattice_3D::ReciprocalOp> GetReciprocalPointOps() const {return {};}

    //! The FULL DETECTED crystal reciprocal point group \f${U|\tau}\f$ -- populated on every periodic run,
    //! IMPOSED or not (unlike \c GetReciprocalPointOps, which is the ops-to-impose set and stays empty on a
    //! free run).  This is what the §3 order-parameter diagnostic (doc/SymmetryUpgradePlan.md) measures the
    //! converged density against: a free run reports the symmetry it actually found, per op
    //! (\c SymmetryDefects).  Default {} = no lattice / nothing detected.
    virtual std::vector<Symmetry::Lattice_3D::ReciprocalOp> GetDetectedReciprocalOps() const {return {};}

    // Iterate() with no type argument yields the base obs_t* directly (no cast);
    // Iterate<D>() dynamic_cast's each IBS to the requested derived type D.
    // Built on the two primitives below, so storage stays private to the
    // concrete BasisSet.
    auto Iterate() const {return IndexProxy<const tBasisSet>(this, GetNumIBS());}
    template <class D> auto Iterate() const
    {
        return D_IndexProxy<const D, const tBasisSet>(this, GetNumIBS());
    }
    const obs_t* operator[](size_t i) const {return GetIBS(i);}

    virtual size_t GetNumIBS() const=0; //!< Number of Irrep basis sets (public since Step 3c-2: the
                                        //!< mixed-aware per-index walk needs the count beside GetRealIBS).
    //! \brief The CROSS-SCALAR block view (doc/RealComplexPlan.md Step 3c-2): non-null iff block \a i is a
    //! REAL block inside a complex-faced set (a TRIM block after the Step-3c-3 factory decision).  The
    //! same-scalar face stays \c GetIBS / \c Iterate; a mixed-aware consumer (MakeIrrepWFs, the GPW
    //! preflight) checks this FIRST and falls back to \c GetIBS.  Default null: a homogeneous set -- and
    //! every real-faced set, whose blocks are already real via \c GetIBS -- has no cross-scalar children.
    virtual const Orbital_1E_IBS<double>* GetRealIBS(size_t) const {return nullptr;}

protected:
    virtual const obs_t* GetIBS(size_t) const=0; //The only storage-specific primitive.
};

typedef tBasisSet<double>    Real_BS;
typedef tBasisSet<dcmplx> Complex_BS;

}//namespace