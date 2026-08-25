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
// What it provides: the cached metric accessors (Repulsion / InvRepulsion / InvOverlap / Norm / Charge)
// and the numerical integrals over the quadrature mesh it OWNS and does not hand out.  A Gaussian
// auxiliary basis is a family of FUNCTIONS; the mesh is the device it computes its own integrals with.
module;
#include <memory>
#include <string>
export module qchem.BasisSet.Internal.Fit_IBS;
export import qchem.BasisSet.Orbital_DFT_IBS;   // FIT_CD_NonOrtho / FIT_SF_NonOrtho -- the faces it implements
export import qchem.Mesh;                       // qcMesh::Mesh / MeshParams -- the mesh it owns
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
