// File: BasisSet/Fit_IBS.C  Interfaces for a fitting (auxiliary) Basis Set.
module;
#include <string>
export module qchem.BasisSet.Fit_IBS;
export import qchem.BasisSet.IrrepBasisSet;
export import qchem.ScalarFunction;
export import qchem.Mesh;            // qcMesh::Mesh / MeshParams -- the fit quadrature mesh + knobs
import qchem.Structure;               // Structure (SetMesh builds the Becke mesh from it)

export namespace qchem::BasisSet
{

//! \brief The MINIMAL, metric-neutral face of a CHARGE-DENSITY fit basis: just "I am a density-fit basis" --
//! its fit FUNCTIONS, via \c IrrepBasisSet<T>.  This is what \c CreateCDFitBasisSet returns and what the
//! 3-centre \c Repulsion3C consumes (it needs the functions that define each \f$f_c\f$, not their metric).
//! Templated on the representation \a T so a real (Gaussian) fit basis is \c rFIT_CD_ABS (=FIT_CD_ABS<double>,
//! real \c VectorFunction) and a plane-wave one is \c cFIT_CD_ABS (=FIT_CD_ABS<dcmplx>, the complex \f$e^{iG
//! \cdot r}\f$ functions -- honestly complex, no NA-stub).  The two design axes are ORTHOGONAL: this T axis is
//! the representation; the \c FIT_CD_NonOrtho refinement below is the metric axis.  (ISP sibling of \c
//! FIT_SF_ABS.)
template <class T> class FIT_CD_ABS
    : public virtual IrrepBasisSet<T>
{
public:
    //! \brief Is this fit basis ORTHONORMAL (metric = I)?  The metric axis, declared as a mandatory contract
    //! so the fitter Factory can pick the right fitter WITHOUT interrogating concrete identity: \c false selects
    //! the metric-solve (non-ortho) fitter and GUARANTEES the object IS-A \c FIT_CD_NonOrtho; \c true selects the
    //! orthonormal (plane-wave, projection-IS-the-fit) fitter.  Every fit basis must declare its metric.
    virtual bool isOrtho() const=0;
};
using rFIT_CD_ABS = FIT_CD_ABS<double>;  //!< real (Gaussian/Slater/BSpline) density-fit basis
using cFIT_CD_ABS = FIT_CD_ABS<dcmplx>;  //!< complex (plane-wave, G-space) density-fit basis

//! \brief A NON-orthonormal (Gaussian/Slater/BSpline) density-fit basis: adds the Coulomb metric-solve inputs
//! the least-squares density fit needs.  Density fitting solves \f$\min_c \|\rho-\sum_c c_c f_c\|_V\f$ in the
//! Coulomb norm under a charge constraint, so this face serves the Coulomb metric \c Repulsion (the
//! \f$\langle f_a|1/r_{12}|f_b\rangle\f$ system matrix), its inverse, the cross-repulsion against another CD
//! fit basis (self-energy), and the per-function \c Charge.  SOLE consumer: the non-ortho \c ConstrainedFF
//! density fitter.  It refines \c rFIT_CD_ABS -- a non-orthonormal fit basis is inherently REAL (there are no
//! complex non-ortho fit bases); an orthonormal plane-wave basis omits ALL of this (the projection IS the fit).
class FIT_CD_NonOrtho
    : public virtual rFIT_CD_ABS
{
public:
    virtual const  rvec_t& Charge      () const=0;  //!< <f_a|1> per fit function (the charge constraint RHS)
    virtual const rsmat_t& Repulsion   () const=0;  //!< Coulomb metric <f_a|1/r12|f_b>, cached
    virtual const  rmat_t& Repulsion   (const rFIT_CD_ABS&) const=0; //!< cross Coulomb <f_a|1/r12|g_b> (arg = functions)
    virtual const rsmat_t& InvRepulsion() const=0;  //!< inverse of the Coulomb metric, cached
};

//! \brief The MINIMAL, metric-neutral face of a SCALAR-FUNCTION (potential) fit basis: just "I am a
//! potential-fit basis" -- its fit FUNCTIONS, via \c IrrepBasisSet<T>.  This is what \c CreateVxcFitBasisSet
//! returns.  Templated on the representation \a T (mirror of \c FIT_CD_ABS): real \c rFIT_SF_ABS (Gaussian)
//! and complex \c cFIT_SF_ABS (plane-wave).  The T axis is the representation; the \c FIT_SF_NonOrtho
//! refinement below is the overlap-metric axis.  (ISP sibling of \c FIT_CD_ABS -- identical shape.)
template <class T> class FIT_SF_ABS
    : public virtual IrrepBasisSet<T>
{
public:
    //! \brief Is this fit basis ORTHONORMAL (metric = I)?  Mirror of \c FIT_CD_ABS::isOrtho on the overlap-metric
    //! axis: \c false selects the overlap metric-solve fitter and GUARANTEES the object IS-A \c FIT_SF_NonOrtho;
    //! \c true selects the orthonormal (plane-wave) scalar fitter.  Every fit basis must declare its metric.
    virtual bool isOrtho() const=0;
};
using rFIT_SF_ABS = FIT_SF_ABS<double>;  //!< real (Gaussian/Slater/BSpline) potential-fit basis
using cFIT_SF_ABS = FIT_SF_ABS<dcmplx>;  //!< complex (plane-wave, G-space) potential-fit basis

//! \brief A NON-orthonormal (Gaussian/Slater/BSpline) potential-fit basis: adds the overlap metric-solve
//! inputs the least-squares potential fit needs -- the projection RHS \c Overlap(Sf) \f$=\langle f_a|f\rangle\f$
//! (the field \a f is always the real \f$v_{xc}(\vec r)\f$), the normalisation, the overlap matrix
//! (\c Integrals_Overlap), and the inverse metric \f$S^{-1}\f$ (the fit is \f$c=S^{-1}\langle f|v\rangle\f$).
//! SOLE consumer: the non-ortho \c FunctionFitterImp scalar fitter.  It refines \c rFIT_SF_ABS -- a
//! non-orthonormal fit basis is inherently REAL; an orthonormal plane-wave basis omits ALL of this (\f$S=I\f$,
//! the projection IS the fit -- \c cFIT_SF_ABS stays the empty marker).  Mirror of \c FIT_CD_NonOrtho.
class FIT_SF_NonOrtho
    : public virtual rFIT_SF_ABS
    , public virtual Integrals_Overlap<double>
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
//! (rFIT_CD_ABS for a density fit, FIT_SF_ABS for a potential fit) for type safety; the union exists so
//! one concrete object can be handed to either creator.
class Fit_IBS
    : public virtual FIT_CD_NonOrtho
    , public virtual FIT_SF_NonOrtho
{
public:
    using Integrals_Overlap<double>::Overlap;       // un-hide the metric Overlap() past the Overlap(Sf) override
    using Integrals_Overlap<double>::MakeOverlap;
    //! A Gaussian/Slater/BSpline auxiliary basis is inherently NON-orthonormal (it carries both metric-solve
    //! refinements) -- the single override that satisfies the \c isOrtho contract for BOTH fit faces.
    bool isOrtho() const override {return false;}
    const  rvec_t& Charge   () const override;
    const rsmat_t& Repulsion() const override;
    const  rmat_t& Repulsion(const rFIT_CD_ABS& b) const override;
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
