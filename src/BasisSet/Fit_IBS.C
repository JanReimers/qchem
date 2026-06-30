// File: BasisSet/Fit_IBS.C  Interfaces for a fitting (auxiliary) Basis Set.
module;
#include <string>
export module qchem.BasisSet.Fit_IBS;
export import qchem.BasisSet.IrrepBasisSet;
export import qchem.ScalarFunction;
export import qchem.Mesh;            // qcMesh::Mesh / MeshParams -- the fit quadrature mesh + knobs
import qchem.Structure;               // Structure (SetMesh builds the Becke mesh from it)

export namespace BasisSet
{

//! \brief Auxiliary basis set face for a CHARGE-DENSITY fit (the Coulomb metric): the integrals
//! FittedVee / FittedCD need.  Density fitting solves \f$\min_c \|\rho-\sum_c c_c f_c\|_V\f$ in the
//! Coulomb norm under a charge constraint, so this face serves the Coulomb metric \c Repulsion (the
//! \f$\langle f_a|1/r_{12}|f_b\rangle\f$ system matrix), its inverse, the cross-repulsion against
//! another CD fit basis (self-energy), and the per-function \c Charge.  (ISP sibling of \c FIT_SF_ABS.)
class FIT_CD_ABS
    : public virtual Real_IBS //Real Irrep basis Set
{
public:
    virtual const  rvec_t& Charge      () const=0;  //!< <f_a|1> per fit function (the charge constraint RHS)
    virtual const rsmat_t& Repulsion   () const=0;  //!< Coulomb metric <f_a|1/r12|f_b>, cached
    virtual const  rmat_t& Repulsion   (const FIT_CD_ABS&) const=0; //!< cross Coulomb <f_a|1/r12|g_b>
    virtual const rsmat_t& InvRepulsion() const=0;  //!< inverse of the Coulomb metric, cached
};

//! \brief Auxiliary basis set face for a SCALAR-FUNCTION fit (the overlap metric): the integrals
//! FittedVxc needs.  Potential fitting projects a pointwise field (e.g. \f$v_{xc}(\rho(\vec r))\f$)
//! onto the fit basis in the ordinary overlap norm, so this face serves the overlap metric \c Overlap
//! (inherited \f$\langle f_a|f_b\rangle\f$), its inverse, the projection RHS \c Overlap(Sf)
//! \f$=\langle f_a|f\rangle\f$, and the normalisation.  (ISP sibling of \c FIT_CD_ABS.)
class FIT_SF_ABS
    : public virtual Integrals_Overlap<double>
    , public virtual Real_IBS
{
public:
    using Integrals_Overlap<double>::Overlap;       // the metric <f_a|f_b> (un-hidden past Overlap(Sf))
    using Integrals_Overlap<double>::MakeOverlap;
    typedef ScalarFunction<double> Sf;
    virtual const  rvec_t& Norm   ()            const=0; //!< 1/sqrt(<f_a|f_a>), cached
    virtual        rvec_t  Overlap(const Sf& f) const=0; //!< projection <f_a|f> (the fit RHS; NOT cached)
    virtual const rsmat_t& InvOverlap()         const=0; //!< inverse of the overlap metric, cached
};

//! \brief A fit basis that can do BOTH fits -- the Gaussian auxiliary basis implements all of it.  The
//! concrete-facing union of the two ISP faces; it carries the shared quadrature mesh (built from the
//! Structure in \c SetMesh) and the cached-accessor implementations.  Clients take the narrow face
//! (FIT_CD_ABS for a density fit, FIT_SF_ABS for a potential fit) for type safety; the union exists so
//! one concrete object can be handed to either creator.
class Fit_IBS
    : public virtual FIT_CD_ABS
    , public virtual FIT_SF_ABS
{
public:
    using Integrals_Overlap<double>::Overlap;       // un-hide the metric Overlap() past the Overlap(Sf) override
    using Integrals_Overlap<double>::MakeOverlap;
    const  rvec_t& Charge   () const override;
    const rsmat_t& Repulsion() const override;
    const  rmat_t& Repulsion(const FIT_CD_ABS& b) const override;
    const rsmat_t& InvOverlap() const override;
    const rsmat_t& InvRepulsion() const override;

    //! Build and OWN the fit quadrature mesh (Becke, from the structure).  Called by the
    //! CreateCD/VxcFitBasisSet creators, which already hold the Structure.
    void SetMesh(const Structure&, const qcMesh::MeshParams&);

    // Numerical (mesh-quadrature) versions -- run over the fit basis's OWN mesh (itsMesh).
    const rvec_t& Norm   ()           const override; //!< 1/sqrt(<f_a|f_a>), cached
    rvec_t        Overlap(const Sf& f) const override; //!< projection <f_a|f> (Vxc fit RHS; NOT cached)

protected:
    virtual  rvec_t MakeCharge      () const=0;
    virtual rsmat_t MakeRepulsion   () const=0;
    virtual  rmat_t MakeRepulsion   (const FIT_CD_ABS&) const=0;
    virtual rsmat_t MakeInvOverlap  () const;
    virtual rsmat_t MakeInvRepulsion() const;

    virtual  rvec_t MakeNorm   () const; //Numerical, over itsMesh.

private:
    qcMesh::Mesh itsMesh;   //!< the fit basis's own quadrature mesh.
    std::string  itsMeshID; //!< identity of itsMesh (= MeshParams::ID()); the cache key axis for Norm()
                            //!< so the SAME fit basis built with a DIFFERENT mesh gets a distinct Norm.
};

}//namespace
