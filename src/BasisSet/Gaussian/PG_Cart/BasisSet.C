// File BasisSet/Gaussian/PG_Cart/BasisSet.C
module;
#include <vector>
#include <memory>
#include <functional>
export module qchem.BasisSet.Gaussian.PG_Cart;
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.Internal.Block;
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.PGData;
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.BasisSet.Gaussian.Point.Reader;
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD;  // NR_Evaluator: the IBS IS-A evaluator (base subobject)
import qchem.BasisSet.Gaussian.Point.IBS;                     // Gaussian::Orbital_{1E,DFT,HF}_IBS<E> templated mixins
import qchem.BasisSet.Gaussian.Lattice.LatticeSum1E;            // Gaussian::LatticeSum1E (the GPW periodic-1E seam)
import qchem.UnitCell;                                 // UnitCell (CollocateDensity grid<->cell map)

import qchem.BasisSet.Internal.BasisSetImp;
import qchem.BasisSet.Internal.IrrepBasisSetImp;
import qchem.Structure;
import qchem.Types;
import qchem.BasisSet.Orbital_DFT_IBS;
import qchem.BasisSet.Internal.Fit_IBS;   // the shared Gaussian-family implementation base we DERIVE

export namespace qchem::BasisSet::Gaussian::PG_Cart
{
using namespace ::qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD;  // Cartesian glue moved out to PG_Cart_MnD

rsmat_t MakeOverlap2C  (const PGData* ab);   // EFit 2-centre fit integrals (named radial kernels)
rsmat_t MakeRepulsion2C(const PGData* ab);

class IrrepBasisSet
        : public virtual Evaluatable_IBS<double>,   // a Gaussian AO basis IS evaluatable (was the plain Real_IBS,
                                                     // which since 2026-08-22 no longer promises op(r))
          public IrrepBasisSetImp<double>,
          public Evaluators::PG_Cart_MnD::NR_Evaluator   // IS-A evaluator, which IS-A PGData
    {
    public:
        typedef std::vector<std::unique_ptr<Block>> bv_t;
    
        IrrepBasisSet(Reader *, const Structure *);
        IrrepBasisSet(const rvec_t &exponents, size_t L, const Structure *);
        IrrepBasisSet(const rvec_t &exponents, size_t L);
        virtual ~IrrepBasisSet(); //g++ 15.2 BUG Compiler generated, or inline destructor does instance std::vector templates destructor.

        virtual size_t  GetNumFunctions() const {return PGData::size();}
        virtual size_t  size() const {return PGData::size();}
        virtual rvec_t     operator() (const rvec3_t&) const;
        // d/dR, forward (GeometryOptimizationPlan §13.2; the user's ruling 2026-11-08: "Use Gradient(r)"): the jets come from this block's
        // values and its hand Gradient(r), since ∂φ/∂R_centre = -∇_r φ.  So they need no Enzyme, and the block never refuses.
        virtual bool CanGiveJets() const override {return true;}
        virtual Outcome<size_t,AD::NoGradient> JetsOn(const rvec3vec_t& pts, size_t owner, const Structure* st, qchem::BasisSet::PointJets<double>& out) const override;
        virtual rvec3vec_t Gradient   (const rvec3_t&) const;

        virtual std::string BasisSetID() const {return PGData::BasisSetID();} // geometry-aware cache identity (molecular: no RadialID/AngularID)
        virtual std::string Name     () const;
        virtual std::ostream &Write(std::ostream &) const;

    private:
        
        IrrepBasisSet(const IrrepBasisSet *bs, const bv_t&);
    
        bv_t itsBlocks;
    protected:
    };
// All of the 1E / 3C / 4C integral building is inherited from the Molecule-generic, evaluator-templated
// mixins (instantiated with PG_Evaluator -- the IBS IS-A PG_Evaluator base subobject).  The loops there
// are the basis-agnostic ones; nothing PG-specific remains in the IBS itself.  MakeKinetic returns the
// <p^2>=<-nabla^2> building block (no 1/2 -- the Hamiltonian's; no centrifugal -- Cartesian).
class Orbital_IBS
    : public Gaussian::EOrbital_1E_IBS<Evaluators::PG_Cart_MnD::NR_Evaluator>
    , public Gaussian::Orbital_ERI4_IBS <Evaluators::PG_Cart_MnD::NR_Evaluator>
    , public Gaussian::Orbital_DFT_IBS<Evaluators::PG_Cart_MnD::NR_Evaluator>
    , public Gaussian::EOrbital_PP_IBS<Evaluators::PG_Cart_MnD::NR_Evaluator>
    , public virtual Gaussian::Periodic_Gaussian_IBS   // the GPW periodic seam: all four faces (ISP split 2026-09-08)
    , public virtual Gaussian::TransitionCollocation   // the (k+q, k) transition density + adjoint (LinearResponsePlan B1)
    , public virtual Gaussian::ToleranceBinding        // the builder binds the pair-loop tolerances (D-ENV-TOL option A)
    , public IrrepBasisSet
{
public:
    Orbital_IBS(Reader *, const Structure *);
    Orbital_IBS(const rvec_t& exponents, size_t L, const Structure *);
    Orbital_IBS(const rvec_t& exponents, size_t L);

    virtual rFIT_CD_ABS* CreateCDFitBasisSet(const Structure *) const override;
    virtual rFIT_SF_ABS* CreateVxcFitBasisSet(const Structure *, const qcMesh::MeshParams&) const override;

    //! This basis's AO shells in Cartesian-monomial form (delegates to ExtractAoShells on its own PGData).
    virtual std::vector<Symmetry::Molecule::AoShell> GetAoShells() const override;   // Gaussian::Orbital_1E_IBS

    // Gaussian::LatticeSum1E -- forward to the NR_Evaluator base subobject's lattice-sum kernels (it owns the
    // radials/pols/ns and ENUMERATES the offsets internally per shell pair: there is no cut in R).
    // (The (phase,A) overloads are distinct from the mixin's finite MakeOverlap()/MakeKinetic()/MakeNuclear(cl).)
    virtual chmat_t MakeOverlap(const cellphase_t& phase, const UnitCell& A) const override;
    virtual cvec_t  MakeOverlap(const cellphase_t& phase, const UnitCell& A,
                                const Gaussian::LatticeSum1E::GaussianFunction& g) const override;
    virtual cvec_t  MakeOverlap(const Gaussian::LatticeSum1E::GaussianFunction& g) const override;
    virtual dcmplx  MakeOverlap(const cellphase_t& phase, const UnitCell& A, const Gaussian::LatticeSum1E::GaussianFunction& g,
                                const Gaussian::LatticeSum1E::GaussianFunction& h) const override;
    virtual double  MakeOverlap(const Gaussian::LatticeSum1E::GaussianFunction& g, const Gaussian::LatticeSum1E::GaussianFunction& h) const override;
    virtual chmat_t MakeKinetic(const cellphase_t& phase, const UnitCell& A) const override;
    //! Their derivatives (LatticeSum1E's doc), the gradient and the image virial in one walk.  ⚠ TEMPORARY: through the one MnD/libcint
    //! seam (Lattice.LibcintPullbackWorkaround, user ruling 2026-11-10) until Enzyme can differentiate the MnD kernel (EnzymeAD/Enzyme#3747).
    virtual GeometryPullback MakeLatticeOverlapPullback(const chmat_t& Abar, const cellphase_t& phase, const UnitCell& A, const Structure* st) const override;
    virtual GeometryPullback MakeLatticeKineticPullback(const chmat_t& Abar, const cellphase_t& phase, const UnitCell& A, const Structure* st) const override;
    virtual chmat_t MakeNuclear(const cellphase_t& phase, const UnitCell& A, const Structure* cl) const override;
    //! The Bloch orbitals on a point set (the XC-mesh \f$\Phi\f$ table) -- see LatticeSum1E for WHY the
    //! image sum belongs on this side of the seam rather than in the periodic caller.
    virtual void    BlochPointValues(const rvec3vec_t& pts, const cellphase_t& phase, const UnitCell& A,
                                     mat_t<dcmplx>& Phi) const override;
    virtual chmat_t MakeLocalGaussian(const cellphase_t& phase, const UnitCell& A, const Structure* cl,
                                      const std::function<Gaussian::LatticeSum1E::GaussianFunction(int)>& opForZ) const override;
    virtual chmat_t MakeLocalGaussian(const Structure* cl,
                                      const std::function<Gaussian::LatticeSum1E::GaussianFunction(int)>& opForZ) const override;
    virtual double  MaxExponent() const override;   // finest exponent -> the GPW density-grid cutoff floor
    virtual double  MinExponent() const override;   // coarsest exponent -> the multi-grid coarsest level
    virtual const GPWTolerances& Tolerances() const override;
    virtual void    BindTolerances(const GPWTolerances& t) override;   // ToleranceBinding (the builder's door)
    // CP2K analytic multi-grid collocation + its exact adjoint (the GPW density/KS bridge).
    virtual std::vector<rvec_t> CollocateDensity(const chmat_t& D, const cellphase_t& phase, const UnitCell& A,
                                                 const std::vector<ivec3_t>& N_L,
                                                 const std::vector<double>& ecut_L,
                                                 const LatticeScreener& screener,
                                                 double relFieldSharp=-1.0) const override;
    virtual chmat_t IntegratePotential(const std::vector<rvec_t>& V_L, const cellphase_t& phase, const UnitCell& A,
                                       const std::vector<ivec3_t>& N_L,
                                       const std::vector<double>& ecut_L,
                                       const LatticeScreener& screener, double relCutoffScale,
                                       const chmat_t* screenD, double fieldSharpness=0.0,
                                       double relFieldSharp=-1.0,
                                       const std::vector<size_t>* pairLevels=nullptr) const override;
    virtual rvec3vec_t IntegratePotentialPullback(const std::vector<rvec_t>& V_L, const chmat_t& D, const cellphase_t& phase,
                                                  const UnitCell& A, const std::vector<ivec3_t>& N_L,
                                                  const std::vector<double>& ecut_L, const LatticeScreener& screener,
                                                  double relFieldSharp, const Structure* st) const override;
    virtual std::vector<size_t> StaticFieldPairLevels(const std::vector<double>& ecut_L,
                                                      double beta, double lnEps) const override;
    virtual void EmitLatticeSumReport(const UnitCell& A) const override;             // lattice-sum economy readout
    virtual size_t SetStreamSymmetryOps(const std::vector<Symmetry::Lattice_3D::DirectOp>& ops,
                                        const UnitCell& A,
                                        const rvec3_t& kFrac=rvec3_t(0,0,0)) const override;  // T3 stream fold (§6b)
    virtual size_t StreamFoldOrder() const override;                                  // fold-state cache-key input
    virtual std::vector<cvec_t> CollocateTransition(const mat_t<dcmplx>& dD, const cellphase_t& ketPhase, const rvec3_t& q,
                                                    const UnitCell& A, const std::vector<ivec3_t>& N_L,
                                                    const std::vector<double>& ecut_L, double relFieldSharp=-1.0) const override;
    virtual mat_t<dcmplx> IntegrateTransition(const std::vector<cvec_t>& v_L, const cellphase_t& ketPhase, const rvec3_t& q,
                                              const UnitCell& A, const std::vector<ivec3_t>& N_L,
                                              const std::vector<double>& ecut_L, double relFieldSharp=-1.0) const override;
private:
    GeometryPullback LatticePullback(int which, const chmat_t& Abar, const cellphase_t& phase, const UnitCell& A, const Structure* st) const;
};
// Use E prefix to avoid name clash with the interface class Fit_IBS
class EFit_IBS
    : public virtual Fit_IBS
    , public Symmetry::IdentityRepresentation<double> 
    , public IrrepBasisSet
{
public:
    // A fit set has no matrix-image client: its representation is the identity at E only, refusing any other op (S3).
    //! \a name = the auxiliary set this was read from ("A1_coul"), a fact only the creating factory knows.
    EFit_IBS(Reader *, const Structure *, const std::string& name);                            //!< the Coulomb fit: no mesh
    EFit_IBS(Reader *, const Structure *, const qcMesh::MeshParams&, const std::string& name); //!< the Vxc fit: its projection mesh

    //! The run report's name for this fit.  Not the cache identity: \c BasisSetID stays the shells' own.
    virtual std::string Name() const override {return itsName;}
    virtual rsmat_t MakeOverlap() const {return MakeOverlap2C(this);}
    // d/dR (GeometryOptimizationPlan §13.2): no geometry derivative for this fit block yet (the Dunlap increment adds the fit faces').
    virtual bool CanPullbackOverlap() const override {return false;}
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeOverlapPullback(const hmat_t<double>&, const Structure*) const override
        {return qchem::BasisSet::NoGeometryDerivativeYet<rvec3vec_t>("the PG_Cart fit block", "overlap pullback");}
    // The Coulomb metric's pullback (the Dunlap force, §13.2): libcint's int2c2e_ip1 on this block's own PGData (Imp/Fit_IBS.C).
    virtual bool CanPullbackRepulsion() const override {return true;}
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeRepulsionPullback(const rsmat_t& A, const Structure* st) const override;
    virtual  rvec_t MakeCharge   () const;
    virtual rsmat_t MakeRepulsion() const {return MakeRepulsion2C(this);}
    virtual  rmat_t MakeRepulsion(const rFIT_CD_ABS& b) const;
private:
    std::string itsName;
};
class BasisSet 
    : public virtual ::qchem::BasisSet::tBasisSet<double>
    , public ::qchem::BasisSet::BasisSetImp<double>
    , public virtual Gaussian::ToleranceBinding   // the builder binds every block's pair-loop tolerances (D-ENV-TOL option A)
{
public:
    BasisSet() {};
    BasisSet( Reader*, const Structure*);
    virtual void Insert(obs_t* bs);
    virtual void BindTolerances(const GPWTolerances& t) override;   //!< each block that can take them

};

} //namespace


