// File: BasisSet/Internal/SymmetryAdapted_IBS.C  Symmetry-adapted (SALC) decorator over an orbital IBS.
//
// Stage 4 of the molecular-symmetry plan (doc/MolecularSymmetryPlan.md): one point-group
// irrep of a molecular basis, presented as a normal IrrepBasisSet so the SCF/accelerator
// stack is none the wiser.  It wraps the raw (whole-molecule) AO basis and an O_Gamma column
// block (nAO x dGamma) of the SALC transform, and returns every 1-electron matrix transformed
// to that irrep block:  M_Gamma = O_Gamma^T M_raw O_Gamma.  The raw matrix is the public,
// already-cached accessor (computed once, shared by all irreps); the transformed block is
// memoized by the existing cache under an irrep-specific AngularID -- both cached, no new
// cache code.  (The 2-electron Fock path -- build F_AO, slice per irrep -- comes next.)
// ★ WHY `.Internal.` (V1.20, ruled 2026-09-09).  Its ONLY consumers are
// `qchem.BasisSet.Gaussian.Point.SymmetryAdaptedBasisSet` and a molecular unit test — both inside the
// `src/BasisSet/` tree, but in a different CMake target (`qcMolecule_BS`).  The user ruled that the
// **qcBasisSet* family counts as ONE library** for CLAUDE.md's "do not import internals across library
// BOUNDARIES" rule: `.Internal.` marks the family boundary, not the target boundary.  So this decorator —
// which no consumer outside the basis-set family has any business naming — moves in, and the rule keeps
// its teeth where it matters.  (The genuine violations the ruling exposes are the
// `qchem.BasisSet.GMap` imports from `qcChargeDensity` and `qcFitting`, which are unrelated
// libraries; filed as V1.20b.)
module;
#include <memory>
#include <vector>
#include <string>
#include <iosfwd>
#include <map>
#include <memory>
export module qchem.BasisSet.Internal.SymmetryAdapted_IBS;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Orbital_HF_IBS;         // HF 2-electron CONTRACTION face (no ERI4 -- R1.7)
export import qchem.BasisSet.Orbital_PP_IBS;         // species fields, forwarded to the raw basis
export import qchem.BasisSet.Orbital_DFT_IBS;        // DFT 3-centre mixin (fitted Coulomb / Vxc)
import qchem.BasisSet.Internal.IrrepBasisSetImp;
import qchem.Structure;
import qchem.Types;

export namespace qchem::BasisSet
{

// SymFockCache is GONE (V1.31).  It memoized the AO Coulomb/exchange per cd-irrep, keyed by BasisSetID and
// invalidated by an ELEMENTWISE density compare, to stop the composite's pair loop from doing ~N^2
// whole-molecule AO builds.  The loop is the wrong shape for this basis, not the builds: with the
// WholeSystemFock_IBS route below the composite assembles ONE AO density, this class builds ONCE, and every
// irrep slices that -- so there is nothing left to memoize.  See doc/CleanupCandidates.md V1.31.

class SymmetryAdapted_IBS
    : public virtual Orbital_1E_IBS<double>
    , public virtual Orbital_HF_IBS<double>  // HF Coulomb/exchange CONTRACTION only -- no ERI4 substrate
    , public virtual WholeSystemFock_IBS<double>  // ...and it builds the AO Fock ONCE, then slices (V1.31)
    , public virtual Orbital_DFT_IBS<double> // DFT 3-centre fitted Coulomb / Vxc
    , public virtual Orbital_PP_IBS<double>  // species fields: the raw basis's matrix, transformed
    , public IrrepBasisSetImp<double>        // provides GetSymmetry / GetSymt / GetIrrep
    , public Symmetry::IdentityRepresentation<double>   // a 1-D irrep block: M_g = χ(g)·1, |χ|² = 1, so its bilinear images are
                                                         // identities -- for the ops of the group it was adapted under (S3)
{
public:
    // raw: the whole-molecule AO basis (not owned).  Oblock: this irrep's SALC columns
    // (nAO x dGamma).  label: Mulliken irrep label (used for the cache key + display).
    // group: the ops the SALCs were adapted under (SALCs::ops), SHARED by every block of the basis; {} = {E}.
    SymmetryAdapted_IBS(const Orbital_1E_IBS<double>* raw, const rmat_t& Oblock,
                        const std::string& label, const sym_t& sym,
                        std::shared_ptr<const std::vector<Symmetry::SymOp>> group={});

    virtual size_t GetNumFunctions() const {return itsO.columns();}
    const rmat_t&  GetO() const {return itsO;}        // this irrep's SALC columns

    // 1-electron integrals in the irrep basis (O^T M_raw O).
    virtual rsmat_t MakeOverlap()                 const;
    virtual rsmat_t MakeKinetic()                 const;
    virtual rsmat_t MakeNuclear(const Structure* cl) const;
    // d/dR (GeometryOptimizationPlan §13.2).  M_Γ = Oᵀ M_raw O with the SALC columns O FIXED (they are symmetry-operation images, not
    // functions of R), so Tr(D ∂M_Γ/∂R) = Tr(O D Oᵀ ∂M_raw/∂R): every pullback below is the RAW basis's, handed this block's density
    // folded up to AO (AddAODensity), and a block can give one exactly when its raw basis can.  Summed over the irreps the AO densities
    // add to the whole-molecule D, so the force is the AO force of the same density -- along a symmetry-breaking direction too, since a
    // symmetric SCF density is stationary against the off-irrep rotations the blocking left out (F has no off-irrep block).
    virtual bool CanPullbackOverlap() const override;
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeOverlapPullback(const hmat_t<double>& D, const Structure* st) const override;
    virtual bool CanPullbackKinetic() const override;
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeKineticPullback(const hmat_t<double>& D, const Structure* st) const override;
    //! The virials forward through the same fold (CellHalf_Design §5): the SALC map is geometry- and strain-free.
    virtual bool CanVirialOverlap() const override;
    virtual Outcome<rmat3d_t,AD::NoGradient> MakeOverlapVirial(const hmat_t<double>& D, const Structure* st) const override;
    virtual bool CanVirialKinetic() const override;
    virtual Outcome<rmat3d_t,AD::NoGradient> MakeKineticVirial(const hmat_t<double>& D, const Structure* st) const override;
    virtual bool CanPullbackNuclear() const override;
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeNuclearPullback(const hmat_t<double>& D, const Structure* st) const override;
    // Species fields (O^T M_raw O), when the raw basis can assemble them.
    virtual bool    CanAssembleField     (const SpeciesRadialField&)  const;
    virtual bool    CanAssembleProjectors(const SpeciesProjectorSet&) const;
    virtual rsmat_t MakeSpeciesFieldMatrix(const Structure*, const SpeciesRadialField&, FieldRange) const;
    virtual rsmat_t MakeProjectorMatrix   (const Structure*, const SpeciesProjectorSet&) const;
    // d/dR (GeometryOptimizationPlan §13.2): the raw basis's, on the AO-folded D (as the 1E pullbacks above).
    virtual bool CanPullbackSpeciesField(const SpeciesRadialField&)  const override;
    virtual bool CanPullbackProjectors  (const SpeciesProjectorSet&) const override;
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeSpeciesFieldPullback(const Structure*, const SpeciesRadialField&, FieldRange, const hmat_t<double>&) const override;
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeProjectorPullback(const Structure*, const SpeciesProjectorSet&, const hmat_t<double>&) const override;

    // 2-electron (HF): build the AO Coulomb/exchange from the cd-irrep's density block (linear
    // in D, so no 4-index ERI transform) and slice to this irrep.  Summed over the cd irreps by
    // the charge density, this yields O^T J_AO(D_total) O.
    virtual void AccumulateDirect  (rsmat_t& Jab, const rsmat_t& Dcd, const Orbital_HF_IBS<double>* bs_cd) const;
    virtual void AccumulateExchange(rsmat_t& Kab, const rsmat_t& Dcd, const Orbital_HF_IBS<double>* bs_cd) const;
    //! The SALC path builds the AO Fock and slices it -- there are no per-irrep-pair ERI4 blocks at all
    //! (this class has no ERI4 substrate face -- R1.7), so the ERI4 bra-ket fusion does not apply.  Fall
    //! back to the two independent AO slices (the whole-AO build already banks the full 8-fold symmetry).
    //! See doc/ERI4Rework.md §5.4.
    virtual void AccumulateDirectBoth(rsmat_t& Ji, rsmat_t& Jj, const rsmat_t& Di, const rsmat_t& Dj,
                                      const Orbital_HF_IBS<double>* cd) const;
    virtual void AccumulateExchangeBoth(rsmat_t& Ki, rsmat_t& Kj, const rsmat_t& Di, const rsmat_t& Dj,
                                        const Orbital_HF_IBS<double>* cd) const;
    // d/dR (GeometryOptimizationPlan §13.2): the raw basis's ERI pullback on BOTH blocks' AO-folded densities, the partner folding its own
    // block through the abstract whole-system face (as PairAOFock does).  Correct for any pair, but each call is one whole-AO pullback: the
    // composite's whole-system route (MakeAOPullback below) folds every block first and calls ONE, so a run never walks these N² times.
    virtual bool CanPullbackERI4(const Orbital_HF_IBS<double>* other) const override;
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeDirectPullback  (const rsmat_t& Da, const rsmat_t& Db, const Orbital_HF_IBS<double>* other, const Structure* st) const override;
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeExchangePullback(const rsmat_t& Da, const rsmat_t& Db, const Orbital_HF_IBS<double>* other, const Structure* st) const override;

    // --- WholeSystemFock_IBS (V1.31): the route the composite density actually takes for a SALC basis ---
    virtual size_t  AODimension() const {return itsO.rows();}
    virtual void    AddAODensity(rsmat_t& Dao, const rsmat_t& D) const;   // Dao += O D O^T
    virtual rsmat_t MakeAOFock  (const rsmat_t& Dao, bool exchange) const;// the ONE whole-AO build
    virtual void    SliceAOFock (rsmat_t& Fab, const rsmat_t& Fao) const; // Fab += O^T Fao O
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeAOPullback(const rsmat_t& Dao, bool exchange, const Structure* st) const; // ∂Tr(Dao X[Dao])/∂R

    // DFT 3-centre fitted Coulomb / Vxc.  The cached accessors Overlap3C/Repulsion3C are inherited
    // from Orbital_DFT_IBS unchanged: they key the *transformed* block under this irrep's
    // AngularID and, on a miss, call our MakeXxx3C below.  MakeXxx3C transforms the raw basis's
    // *cached* 3C (now safe -- the integral cache is re-entrant), so the raw 3C is computed once
    // and shared by every irrep.  Fit bases are atom-centred (geometry, not symmetry), so creation
    // delegates to the raw DFT basis.
    virtual rFIT_CD_ABS* CreateCDFitBasisSet (const Structure*, const qcMesh::MeshParams&) const;
    virtual rFIT_SF_ABS* CreateVxcFitBasisSet(const Structure*, const qcMesh::MeshParams&) const;
protected:
    virtual PairProjection<double> MakeOverlap3C  (const rFIT_SF_ABS& c) const;
    virtual PairProjection<double> MakeRepulsion3C(const rFIT_CD_ABS& c) const;
public:
    // d/dR (GeometryOptimizationPlan §13.2): the raw DFT basis's, on the AO-folded D (the fit basis is the raw one's, atom-centred).
    virtual bool CanPullbackRepulsion3C(const rFIT_CD_ABS& c) const override;
    virtual Outcome<rvec3vec_t,AD::NoGradient> MakeRepulsion3CPullback(const rFIT_CD_ABS& c, const hmat_t<double>& D, const rvec_t& a, const Structure* st) const override;

    // The per-irrep cache distinction lives in BasisSetID() (raw id + "[label]"), so each irrep's
    // transformed blocks key separately.  (No RadialID/AngularID -- those are the atom identity face;
    // a SALC basis is molecular, and keys on the raw basis's own BasisSetID.)
    virtual std::string BasisSetID() const;
    virtual std::string Name()       const;

    // VectorFunction: a SALC evaluated at r is O^T applied to the raw AO values.
    virtual rvec_t     operator()(const rvec3_t&) const;
    // d/dR, forward (GeometryOptimizationPlan §13.2): the gradient source awaits the user's ruling, so every family refuses for now.
    virtual bool CanGiveJets() const override {return false;}
    virtual Outcome<size_t,AD::NoGradient> JetsOn(const rvec3vec_t&, size_t, const Structure*, BasisSet::PointJets<double>&) const override
        {return BasisSet::NoGeometryDerivativeYet<size_t>("a symmetry-adapted (SALC) block", "function jets");}
    virtual rvec3vec_t Gradient  (const rvec3_t&) const;

    virtual std::ostream& Write(std::ostream&) const;

private:
    //! The pair route's AO Fock: the PARTNER folds its own block up (through the abstract whole-system
    //! face), we build once.  V1.10 -- no reach into a sibling's private SALC columns.
    rsmat_t      PairAOFock(const rsmat_t& Dcd, const Orbital_HF_IBS<double>* bs_cd, bool exchange) const;
    rsmat_t      Transform(const rsmat_t& Mraw) const;          // O^T Mraw O, symmetrized
    rsmat_t      FoldUp(const rsmat_t& D) const;                // O D O^T, symmetrized: the pullbacks' AO density
    PairProjection<double> Transform3C(const PairProjection<double>& raw) const; // Transform each fit-function matrix

    const Orbital_1E_IBS<double>*  itsRaw;            // raw whole-molecule AO basis (not owned)
    const Orbital_HF_IBS<double>*  itsRawHF;          // same object, HF interface (for the AO J/K build)
    const Orbital_DFT_IBS<double>* itsRawDFT;         // same object, DFT interface (raw 3C + fit bases)
    rmat_t                        itsO;               // this irrep's SALC columns (nAO x dGamma)
    std::string                   itsLabel;           // Mulliken irrep label
};

} //namespace
