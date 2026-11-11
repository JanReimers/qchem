// File: BasisSet/Internal/Fit_IBS.C  The shared IMPLEMENTATION base of a Gaussian-family fit basis.
//
// INTERNAL since 2026-08-25 (user), and the code had been asking for it: the class carried the comment
// "Design question: Does the client code really need to see this class, or should it just be working with
// abstract interfaces?"  Measured answer: NO, and it already does not.  qcHamiltonian's one alias to the
// concrete union (`using fbs_t = BasisSet::Fit_IBS`) was DEAD -- every term already takes the narrow
// abstract face (rFIT_CD_ABS / rFIT_SF_ABS) -- and deleting it built the whole tree unchanged.  So the
// only things that name this class are the three concrete bases that DERIVE it.
//
// That makes it exactly the shape Internal/ is for, and it follows the precedent of
// Internal.IrrepBasisSetImp: a shared implementation base imported by the concrete bases in the Atom,
// Molecule and Lattice_3D sub-libraries, and by nothing above them.  The abstract faces it implements
// stay public in qchem.BasisSet.Orbital_DFT_IBS:Fit_IBS -- consumers depend on those and never on this.
//
// What it provides: the cached metric accessors (Repulsion / InvRepulsion / InvOverlap / Norm / Charge), all
// ANALYTIC, and -- on a Vxc fit only -- the projection of a scalar function over the quadrature mesh it owns.
// A Gaussian auxiliary basis is a family of FUNCTIONS.  The Coulomb (rho) fit owns no mesh at all (user
// 2026-11-11: "The rho fit should be 100% analytic, no grid"): every integral it serves is analytic, and the
// one thing that needs points -- projecting a matrix-free field, the SAD seed -- names them in the fitter.
module;
#include <memory>
#include <optional>
#include <string>
export module qchem.BasisSet.Internal.Fit_IBS;
export import qchem.BasisSet.Orbital_DFT_IBS;   // FIT_CD_NonOrtho / FIT_SF_NonOrtho -- the faces it implements
export import qchem.Mesh;                       // qcMesh::Mesh / MeshParams -- a Vxc fit's projection mesh
import qchem.Structure;                         // Structure (the ctor builds that mesh from it)

export namespace qchem::BasisSet
{

//! \brief One class that implements all four abstract interfaces -- the shared Gaussian-family base.
//! The design question this comment used to ask -- "does the client code really need to see this class,
//! or should it just be working with abstract interfaces?" -- is ANSWERED (2026-08-25): it should work
//! with the abstract interfaces, it already did, and the class is now Internal so it cannot do otherwise.
class Fit_IBS
    : public virtual FIT_CD_NonOrtho
    , public virtual FIT_SF_NonOrtho
    //! PUBLIC since V-XC-ONE-TERM step 3 (lane S, 2026-11-10; was protected).  Client: the overlap-metric PROJECTED XC sampler
    //! (\c ChargeDensity::OverlapProjectedDensitySampler), which needs \f$f_a(r_g)\f$ at the XC mesh's points to evaluate \f$P\rho\f$
    //! there and to project \f$v\f$ back with the SAME points and weights -- the exact adjoint.  The basis's own \c Overlap(f) would
    //! not do: it integrates on this class's PRIVATE mesh, so \f$H\f$ would be the adjoint of a different quadrature than \f$E\f$'s.
    //! The client reaches it by the house cross-cast between abstract faces (rFIT_SF_ABS -> Evaluatable_IBS<double>), never by
    //! naming this class.
    , public virtual Evaluatable_IBS<double> //Support op(r)
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
    //! \brief A COULOMB (rho) fit: analytic, no mesh (user 2026-11-11, "The rho fit should be 100% analytic, no grid").
    //! It took the XC MeshParams until then and built a private mesh from them, for \c Norm() and \c Overlap(f); both had
    //! an analytic answer or belong to the field: a field that is only a function on space integrates itself on points it owns
    //! (\c SelfIntegratingField; the SAD seed).
    Fit_IBS();
    //! \brief A Vxc fit: build and OWN its PROJECTION mesh (from the structure) at construction, for \c Overlap(f) -- the
    //! one integral here that is a quadrature, \f$\langle f_a|v\rangle\f$ of a pointwise field (R2.10: the creator holds
    //! the Structure and the MeshParams and passes them down).
    Fit_IBS(const Structure&, const qcMesh::MeshParams&);

public:
    //! \f$1/\sqrt{\langle f_a|f_a\rangle}\f$ of the functions as \c op(r) evaluates them: ANALYTIC, from the metric's
    //! diagonal (1 to roundoff for the analytically normalised Gaussian/Slater/BSpline shells).  It was a mesh quadrature,
    //! keyed by the mesh in the cache after a coarse-mesh value leaked into a fine run (a 585 ppm drift); a norm is a
    //! property of the function, so it has no mesh to key on.
    const rvec_t& Norm   ()           const override;
    //! \brief REFUSED: a non-orthogonal fit has no honest diagonal-only answer.  Doc/CleanupCandidates.md R1.0 (user ruling
    //! 2026-08-23): this is to be DELETED from the metric-neutral face, never corrected -- a corrected value is all ones,
    //! indistinguishable from the plane-wave answer.  It returned \f$1/\mathrm{Norm}_a^2\f$ of the MESH norm until the norm
    //! went analytic (2026-11-11), when it would have become exactly that corrected value; nothing reads it (a Gaussian fit
    //! takes the full \f$S^{-1}\f$ solve), so it throws until R1.0 moves it onto an orthogonal refinement.
    rvec_t OverlapDiagonal() const override;
    //! \copydoc BasisSet::FIT_SF_ABS::Overlap  (a mesh quadrature over the VXC fit's projection mesh, in the NORMALISED
    //! convention, \f$\langle\hat f_a|f\rangle\f$ with \f$\hat f_a=f_a\,\mathrm{Norm}_a\f$, which \c InvOverlap()'s metric and
    //! \c Charge() are also in).  THROWS on a Coulomb fit, which owns no mesh.
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

    virtual  rvec_t MakeNorm   () const; //!< analytic: 1/sqrt(diag of the overlap metric)

private:
    std::optional<qcMesh::Mesh> itsProjectionMesh;   //!< a Vxc fit's projection mesh; EMPTY on a Coulomb fit
};

}//namespace
