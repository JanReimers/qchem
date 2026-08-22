// File: BasisSet/DeltaFit_IBS.C  The delta (identity) v_xc fit basis over a real-space mesh.
//
// The OTHER representation a periodic Vxc fit basis can take (PlaneWaveFit_IBS being the {G} one): the
// delta basis, whose n_pts functions are the mesh's own delta functions.  Everything it knows is the
// quadrature -- points, weights, and the orbit fold the run's imposed ops gave that mesh -- because that
// IS what a delta basis is: W_g[h] = w_g delta_gh, a family of weight vectors with nothing left over.
//
// It exists so that CreateVxcFitBasisSet can make BOTH decisions -- which representation AND which points
// -- and return ONE object (doc/OpenWork.md, "separation of concerns in the XC terms").  Before this the
// mesh came back through a SECOND factory (CreateXCQuadrature) as a bare struct, and the Hamiltonian had
// to re-decide which term class to build in order to know which of the two answers to use; the fit basis
// and the grid were separate returns of one question.
//
// It ANSWERS, it does not hand out (user ruling 2026-08-22).  The mesh and its orbit fold are private:
// what leaves this class is an integral, a per-site integral, a sampled field, or a symmetrized field --
// operations over value arrays in its own point ORDER.  No consumer knows whether these points are an
// atom-centred Becke build or a uniform cell grid, which is exactly the fit/grid separation the design
// item is about.  (The one survivor is Quadrature::Mesh(); see the \todo on FIT_SF_Delta.)
module;
#include <cassert>
#include <ostream>
#include <string>
#include <utility>
export module qchem.BasisSet.DeltaFit_IBS;
export import qchem.BasisSet.Fit_IBS;                    // cFIT_SF_Delta (the face) + FitQuadrature (the bundle)
import qchem.BasisSet.Internal.IrrepBasisSetImp;         // GetSymmetry/GetSymt/GetIrrep + itsSymmetry
import qchem.Symmetry;                                   // sym_t (the Bloch irrep)
import qchem.Types;                                      // dcmplx, rvec3_t, vec_t
import qchem.Mesh.Quadrature;                            // qcMesh::Integrate / SiteIntegrals (my own quadrature)
import qchem.Symmetry.Lattice_3D.Fold;                   // SymmetrizeValues / SymmetrizeValuesSigned (my orbit fold)
import qchem.Reporting;                                  // EmitFold -- this basis announces its own star-average

export namespace qchem::BasisSet
{

//! \brief The \f$\delta\f$-basis realization of the periodic \f$v_{xc}\f$ fit face: it holds the finished
//! XC quadrature and answers nothing else.
//!
//! Symmetry is the \f$\Gamma\f$ Bloch irrep (user ruling): the mesh is cell-periodic and carries no
//! crystal momentum of its own, exactly as \c PlaneWaveFit_IBS's fit block is built at \f$k=0\f$.
//! \c BasisSetID folds in the point count and the fold's orbit count -- the only things that distinguish
//! two of these -- but nothing keys an integral cache on a \f$\delta\f$ basis today: it has no 3-centre
//! tensor, because the quadrature route never forms one.
class DeltaFit_IBS
    : public virtual cFIT_SF_Delta            // FIT_SF_Delta<dcmplx> : cFIT_SF_ABS + Quadrature
    , public         IrrepBasisSetImp<dcmplx> // GetSymmetry/GetSymt/GetIrrep
{
public:
    //! \a q is the FINISHED quadrature (mesh + fold + Shubnikov tags) from the owning basis's
    //! \c CreateXCQuadrature, and \a sym the \f$\Gamma\f$ Bloch irrep.
    DeltaFit_IBS(FitQuadrature q, const sym_t& sym)
        : IrrepBasisSetImp<dcmplx>(sym)
        , itsQuad(std::move(q))
    {
        assert(itsQuad.mesh && itsQuad.mesh->size()>0 && "DeltaFit_IBS: a delta basis IS its mesh -- it cannot be empty");
        // The invariants of the bundle, checked HERE, once, where it becomes an object -- they used to be
        // re-checked by the consumer that received the loose struct.
        assert(itsQuad.fold.owner.empty()   || itsQuad.fold.owner.size()==itsQuad.mesh->size());
        assert(itsQuad.sigmas.empty()       || !itsQuad.fold.owner.empty());
        assert(itsQuad.flipFixed.empty()    || itsQuad.flipFixed.size()==itsQuad.mesh->size());
        // Providers self-report (user ruling 2026-08-16): the star-average this basis will apply to rho on
        // EVERY iteration announces itself at birth, armed or not.  nOps is known only for a magnetic fold
        // (sigmas run parallel to the ops); 0 = "not reported", NOT "not armed".
        qchem::report::EmitFold("XC mesh", itsQuad.sigmas.size(), itsQuad.mesh->size(),
                                itsQuad.fold.owner.empty() ? itsQuad.mesh->size() : itsQuad.fold.Reps(),
                                itsQuad.sigmas.empty() ? std::string() : std::string("magnetic (Shubnikov)"));
    }

    //! \copydoc BasisSet::Quadrature::Mesh  (the one remaining getter -- see the \todo on \c FIT_SF_Delta)
    const qcMesh::Mesh& Mesh() const override {return *itsQuad.mesh;}

    // ---- the quadrature operations (see FIT_SF_Delta) -------------------------------------------------
    double Integrate    (const rvec_t& f) const override {return qcMesh::Integrate   (*itsQuad.mesh, f);}
    rvec_t SiteIntegrals(const rvec_t& f) const override
        {return itsQuad.mesh->NSites()==0 ? rvec_t() : qcMesh::SiteIntegrals(*itsQuad.mesh, f);}
    rvec_t Sample(const ScalarFunction<double>& f) const override {return f(itsQuad.mesh->Points());}

    void Symmetrize(rvec_t& f) const override
    {
        if (itsQuad.fold.owner.empty()) return;              // free run: the projector is the identity
        Symmetry::Lattice_3D::SymmetrizeValues(itsQuad.fold, f);
    }
    void SymmetrizeSpin(rvec_t& rho, rvec_t& m) const override
    {
        if (itsQuad.fold.owner.empty()) return;
        if (itsQuad.sigmas.empty()) {                        // grey/free: the two channels are independent
            Symmetry::Lattice_3D::SymmetrizeValues(itsQuad.fold, rho);
            Symmetry::Lattice_3D::SymmetrizeValues(itsQuad.fold, m);
            return;
        }
        // Shubnikov S3: rho is EVEN under sigma, m is ODD -- and m must vanish exactly at points a
        // sigma=Flip op maps onto themselves, where the exact projector annihilates it.  Per-CHANNEL
        // averaging would be wrong here: a Flip op maps rho_up onto rho_dn, not onto itself.
        for (size_t g=0; g<itsQuad.flipFixed.size(); ++g) if (itsQuad.flipFixed[g]) m[g]=0.0;
        Symmetry::Lattice_3D::SymmetrizeValues      (itsQuad.fold, rho);
        Symmetry::Lattice_3D::SymmetrizeValuesSigned(itsQuad.fold, itsQuad.sigmas, m);
    }

    //! One \f$\delta\f$ function per mesh point.  ~100k of them: see the sizing warning on \c FIT_SF_Delta.
    size_t GetNumFunctions() const override {return itsQuad.mesh->size();}

    // NO op(r)/Gradient, and no stub either: since 2026-08-22 IrrepBasisSet does not promise pointwise
    // evaluation (doc/CleanupCandidates.md R1.0 step 1), so a basis whose functions are DISTRIBUTIONS
    // simply does not offer it.  This class used to carry two throwing overrides for exactly that reason
    // -- the interface demanded an answer it could not honestly give.  Nothing asks: the delta route's
    // fit coefficients ARE the field's values at the mesh points, and it has no expansion to plot, so it
    // does not implement FieldEvaluator either.

    std::string   Name      () const override {return "DeltaFit";}
    std::string   BasisSetID() const override
        {return Name()+"|npts="+std::to_string(itsQuad.mesh->size())
                      +"|orbits="+std::to_string(itsQuad.fold.owner.empty() ? itsQuad.mesh->size() : itsQuad.fold.Reps());}
    std::ostream& Write     (std::ostream& os) const override
        {return os << Name() << " fit IBS: " << itsQuad.mesh->size() << " delta functions on the XC mesh";}

private:
    FitQuadrature itsQuad;   //!< the mesh + its orbit fold + the Shubnikov spin tags
};

} //namespace
