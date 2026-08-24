// File: BasisSet/Fit_IBS.C  Interfaces for a fitting (auxiliary) Basis Set.
module;
#include <memory>
#include <string>
export module qchem.BasisSet.Fit_IBS;
export import qchem.BasisSet.IrrepBasisSet;
export import qchem.ScalarFunction;
export import qchem.Mesh;            // qcMesh::Mesh / MeshParams -- the fit quadrature mesh + knobs
import qchem.Structure;               // Structure (the ctor builds the quadrature mesh from it)
export import qchem.BasisSet.Orbital_1E_IBS;  // Orbital_1E_IBS<U> -- the block Integrals_Overlap3C takes
export import qchem.BasisSet.Internal.Projector3;  // Projector3<U> -- the house CONTRACTIBLE 3-centre object

export namespace qchem::BasisSet
{

//! \brief Common base interface for the two 3-centre overlap integrals.
template <class U> class Integrals_Overlap3C
{
public:
    virtual ~Integrals_Overlap3C() = default;
    virtual const Projector3<U>& Overlap3C(const Orbital_1E_IBS<U>& orb) const=0;
};

//! Forward declaration of the two types the FIT_SF_ABS::Overlap3C default needs (defined below /
//! in the implementation unit).
template <class T> class FIT_SF_ABS;
//! \brief \c FIT_SF_ABS::Overlap3C's default body, as a free template.
//!
//! DEFINED in the IMPLEMENTATION unit (Imp/Fit_IBS.C) and explicitly instantiated there, because it
//! needs \c Orbital_DFT_IBS -- whose interface imports THIS one, so importing it back here would close
//! a module cycle.  An implementation unit may import it; the interface may not.
template <class T> const Projector3<T>& OrbitalOverlap3C(const Orbital_1E_IBS<T>& orb, const FIT_SF_ABS<T>& fit);


//! \brief Abstract interface for an orthogonal (not necessarily normalized) charge density auxilliary fit basis set.
//! No behaviour beyond IrrepBasisSet is defined at this level.
template <class T> class FIT_CD_ABS
    : public virtual IrrepBasisSet<T>
{
public:
    virtual bool isOrtho() const=0; //!< Is the metric/self-overlap DIAGONAL?, \f$\langle f_a|f_b\rangle=0\f$ for a!=b.
};
using rFIT_CD_ABS = FIT_CD_ABS<double>;  //!< real (Gaussian/Slater/BSpline) density-fit basis
using cFIT_CD_ABS = FIT_CD_ABS<dcmplx>;  //!< complex (plane-wave, G-space) density-fit basis

//! \brief Extension of FIT_CD_ABS for non-orthogonal (Gaussian/Slater/BSpline) density-fit basis sets.
//! For molecules the charge density fit error is minimized by using the Dunlap algorithm (B. I. Dunlap, J. W. D. Connolly, and J. R. Sabin, The Journal of Chemical Physics 71, 3396 (1979)). 
//! Density fitting solves \f$\min_c \|\rho-\sum_c c_c f_c\|_V\f$ in the
//! Coulomb norm under a charge constraint, so this face serves the Coulomb metric \c Repulsion (the
//! \f$\langle f_a|1/r_{12}|f_b\rangle\f$ system matrix), its inverse, the cross-repulsion against another CD
//! fit basis (self-energy), and the per-function \c Charge.  SOLE consumer: the non-ortho \c ConstrainedFF
//! density fitter. So far only the real version is required.
class FIT_CD_NonOrtho
    : public virtual rFIT_CD_ABS
{
public:
    
    virtual rvec_t Charge() const=0; //!< \copydoc BasisSet::FIT_SF_ABS::Charge 
    virtual const rsmat_t& Repulsion   () const=0;  //!< Coulomb metric <f_a|1/r12|f_b>, cached.
    virtual const  rmat_t& Repulsion   (const rFIT_CD_ABS& g) const=0; //!< cross Coulomb <f_a|1/r12|g_b>, cached.
    virtual const rsmat_t& InvRepulsion() const=0;  //!< inverse of the Coulomb metric, cached.
};

//! \brief Abstract interface for an orthogonal (not necessarily normalized) auxilliary fit basis set, for fitting  scalar functions.
//! Typical scalar functions that require fitting are Vxc, Vc, Vex potentials.  Charge density is also a scalar function,
//! but we almost never fit the CD using op(r), insterad we take advantage of density matrix representation and use analytic
//! integrals to compute the projection onto the fit basis.  
template <class T> class FIT_SF_ABS
    : public virtual IrrepBasisSet<T>
    , public virtual Integrals_Overlap3C<T>   // Support overlap between this fit,f and and orbital basis block, <chi_i|f|chi_j>
{
public:
    virtual bool isOrtho() const=0; //!< \copydoc BasisSet::FIT_CD_ABS::isOrtho 
    virtual vec_t<T> Overlap(const ScalarFunction<double>& f) const=0; //!< \brief \f$\langle f_a|f\rangle\f$ -- the projection of a scalar function onto this basis, NOT cached.
    //! \brief \f$\langle\chi_i|f_a|\chi_j\rangle\f$ -- the 3-CENTRE overlap between an orbital block's
    //! functions and this fit, used by Hamiltonian terms to form \f$H_{ij}\f$ and \f$E\f$ 
    const Projector3<T>& Overlap3C(const Orbital_1E_IBS<T>& orb) const override
        {return OrbitalOverlap3C(orb, *this);}

    
    virtual vec_t<T> Charge() const=0; //!< \copydoc BasisSet::FIT_SF_ABS::Charge 
    

    //! \brief The DIAGONAL of the fit basis self overlap metric, \f$\langle f_a|f_a\rangle\f$ 
    virtual vec_t<T> OverlapDiagonal() const=0;

    //! \brief This function does not belong in the fit IBS interface.
    //! STAR-AVERAGE an EXPANSION OVER ME, in place, over the crystal point group (the IBZ density
    //! symmetrization): the argument is a coefficient vector, one entry per FUNCTION, and comes back
    //! projected onto the group-invariant subspace.  REAL-space, so it PRESERVES ρ≥0 -- XC stays on the
    //! non-negative ρ_DM samples, never routed onto ρ̃ (doc/GPWPlan1.md item 3).
    //!
    //! ONE operation, two mechanisms, because the basis owns its own geometry: a raster-backed basis
    //! permutes voxels (g→W·g) and applies the glide τ by the FFT shift theorem; a δ basis applies its
    //! mesh's orbit-mean projector.  Default NO-OP -- molecules, unfolded crystals, any fit basis with no
    //! symmetry structure -- so a caller never asks whether symmetry was imposed, it just symmetrizes.
    //! (Was \c SymmetrizeRaster: named for the PW mechanism, which is why the δ route grew a duplicate
    //! declaration of its own before this was noticed -- 2026-08-22.)
    //! \note The two implementations permute FUNCTIONS the group maps onto each other, and for both of
    //! them that IS a permutation of points -- because both are representations whose functions are keyed
    //! by position.  A Gaussian auxiliary basis would symmetrize by permuting shells instead; it takes the
    //! default no-op because a molecular run imposes nothing.
    virtual void Symmetrize(rvec_t&) const {}

    //! \brief This function does not belong in the fit IBS interface.
    //! The MAGNETIC sibling: project the \f$(\rho,m)\f$ PAIR, which is what diagonalizes
    //! \f$\sigma\f$ -- \f$\rho\f$ EVEN under the orbit mean, \f$m\f$ ODD under the \f$\chi\f$-signed
    //! one with the flip-fixed functions zeroed first (Shubnikov S3, doc/SymmetryUpgradePlan.md §7).
    //!
    //! HERE, not on a \f$\delta\f$-only face, since 2026-08-23 (user: why would this differ by representation?).
    //! It does not.  The DEFAULT below is the grey/free semantics -- average each channel independently --
    //! and is BIT-IDENTICAL to the branch the \f$\delta\f$ basis used to run for itself whenever the run
    //! carried no \f$\sigma\f$ tags.  A representation that knows nothing of magnetic symmetry therefore
    //! gets the right answer for free (a Gaussian basis: two no-ops; a raster: two star-averages), and the
    //! ONE override that remains is \f$\delta\f$'s genuinely different Shubnikov projection.
    //!
    //! The old justification for keeping it δ-only -- "δ is the only representation a polarized run can
    //! use" -- was about which representation gets SELECTED, not about what the operation means; a fact
    //! of the Hamiltonian's \c VxcFit::Auto policy has no business shaping a basis face.
    virtual void SymmetrizeSpin(rvec_t& rho, rvec_t& m) const {Symmetrize(rho); Symmetrize(m);}
};
using rFIT_SF_ABS = FIT_SF_ABS<double>;  //!< real (Gaussian/Slater/BSpline) potential-fit basis
using cFIT_SF_ABS = FIT_SF_ABS<dcmplx>;  //!< complex (plane-wave, G-space) potential-fit basis

//! \brief Extension of FIT_SF_ABS for non-orthogonal (Gaussian/Slater/BSpline) density-fit basis sets.
class FIT_SF_NonOrtho
    : public virtual rFIT_SF_ABS
    , public virtual Integrals_Overlap<double> 
{
public:
    using Integrals_Overlap<double>::Overlap;       // the metric <f_a|f_b> (un-hidden past Overlap(Sf))
    using Integrals_Overlap<double>::MakeOverlap;
    using rFIT_SF_ABS::Overlap;                     // ...and the projection <f_a|f>, now inherited
    typedef ScalarFunction<double> Sf;
    virtual const  rvec_t& Norm   ()            const=0; //!< 1/sqrt(<f_a|f_a>), cached
    virtual const rsmat_t& InvOverlap()         const=0; //!< inverse of the overlap metric, cached
};

//! \brief One class that implelements all four abstract interfaces.
//! Design question: Does the client code really need to see this class, or should it just be working with abstract interfaces? 
class Fit_IBS
    : public virtual FIT_CD_NonOrtho
    , public virtual FIT_SF_NonOrtho
    , protected virtual Evaluatable_IBS<double> //Support op(r)
{
public:
    using Integrals_Overlap<double>::Overlap;       // un-hide the metric Overlap() past the Overlap(Sf) override
    using Integrals_Overlap<double>::MakeOverlap;
    //! A Gaussian/Slater/BSpline auxiliary basis is inherently NON-orthonormal (it carries both metric-solve
    //! refinements) -- the single override that satisfies the \c isOrtho contract for BOTH fit faces.
    bool isOrtho() const override {return false;}
    const rsmat_t& Repulsion() const override;
    const  rmat_t& Repulsion(const rFIT_CD_ABS& b) const override;
    const rsmat_t& InvOverlap() const override;
    const rsmat_t& InvRepulsion() const override;

protected:
    //! \brief Build and OWN the fit quadrature mesh (from the structure) AT CONSTRUCTION.
    //!
    //! Every numerical integral this class provides -- \c Norm(), \c Overlap(f) -- runs over that mesh, so
    //! there is no valid state between "constructed" and "has a mesh": a mesh-less fit basis can answer
    //! nothing it exists to answer.  It was a post-ctor \c SetMesh (two-phase construction) whose only guard
    //! was an assert inside each numerical accessor -- i.e. the invariant was re-checked at every use
    //! instead of established once.  The creators (\c CreateCDFitBasisSet / \c CreateVxcFitBasisSet) already
    //! hold the Structure and the MeshParams, so they simply pass them down (R2.10).
    Fit_IBS(const Structure&, const qcMesh::MeshParams&);

public:
    // Numerical (mesh-quadrature) versions -- run over the fit basis's OWN mesh (itsMesh).
    const rvec_t& Norm   ()           const override; //!< 1/sqrt(<f_a|f_a>), cached
    //! \copydoc BasisSet::FIT_SF_ABS::OverlapDiagonal  (\f$1/\mathrm{Norm}_a^2\f$ -- the same numbers,
    //! un-inverted; a Gaussian fit takes the full \f$S^{-1}\f$ solve and never reads just this)
    rvec_t OverlapDiagonal() const override;
    //! \copydoc BasisSet::FIT_SF_ABS::Overlap  (a mesh quadrature over \c itsMesh, in the NORMALISED
    //! convention -- \f$\langle\hat f_a|f\rangle\f$ with \f$\hat f_a=f_a\,\mathrm{Norm}_a\f$, which is what
    //! \c InvOverlap()'s metric and \c Charge() are also in)
    rvec_t        Overlap(const Sf& f) const override; //!< projection <f_a|f> (Vxc fit RHS; NOT cached)
    //! \copydoc BasisSet::FIT_SF_ABS::Charge
    //! ONE override satisfying BOTH fit faces -- \f$\langle f_a|1\rangle\f$ is one quantity, so it is
    //! declared once per face with the same signature and answered once here.
    rvec_t        Charge() const override;

protected:
    virtual  rvec_t MakeCharge      () const=0;
    virtual rsmat_t MakeRepulsion   () const=0;
    virtual  rmat_t MakeRepulsion   (const rFIT_CD_ABS&) const=0;
    virtual rsmat_t MakeInvOverlap  () const;
    virtual rsmat_t MakeInvRepulsion() const;

    virtual  rvec_t MakeNorm   () const; //Numerical, over itsMesh.

private:
    qcMesh::Mesh itsMesh;   //!< the fit basis's own quadrature mesh.
    std::string  itsMeshID; //!< identity of itsMesh (= MeshParams::ID()); the cache key axis for Norm()
                            //!< so the SAME fit basis built with a DIFFERENT mesh gets a distinct Norm.
};

}//namespace
