// File: BasisSet/SymmetryAdapted_IBS.C  Symmetry-adapted (SALC) decorator over an orbital IBS.
//
// Stage 4 of the molecular-symmetry plan (doc/MolecularSymmetryPlan.md): one point-group
// irrep of a molecular basis, presented as a normal IrrepBasisSet so the SCF/accelerator
// stack is none the wiser.  It wraps the raw (whole-molecule) AO basis and an O_Gamma column
// block (nAO x dGamma) of the SALC transform, and returns every 1-electron matrix transformed
// to that irrep block:  M_Gamma = O_Gamma^T M_raw O_Gamma.  The raw matrix is the public,
// already-cached accessor (computed once, shared by all irreps); the transformed block is
// memoized by the existing cache under an irrep-specific AngularID -- both cached, no new
// cache code.  (The 2-electron Fock path -- build F_AO, slice per irrep -- comes next.)
module;
#include <string>
#include <iosfwd>
#include <map>
#include <memory>
export module qchem.BasisSet.SymmetryAdapted_IBS;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Orbital_HF_IBS;         // HF 2-electron CONTRACTION face (no ERI4 -- R1.7)
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
    , public IrrepBasisSetImp<double>        // provides GetSymmetry / GetSymt / GetIrrep
{
public:
    // raw: the whole-molecule AO basis (not owned).  Oblock: this irrep's SALC columns
    // (nAO x dGamma).  label: Mulliken irrep label (used for the cache key + display).
    SymmetryAdapted_IBS(const Orbital_1E_IBS<double>* raw, const rmat_t& Oblock,
                        const std::string& label, const sym_t& sym);

    virtual size_t GetNumFunctions() const {return itsO.columns();}
    const rmat_t&  GetO() const {return itsO;}        // this irrep's SALC columns

    // 1-electron integrals in the irrep basis (O^T M_raw O).
    virtual rsmat_t MakeOverlap()                 const;
    virtual rsmat_t MakeKinetic()                 const;
    virtual rsmat_t MakeNuclear(const Structure* cl) const;

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

    // --- WholeSystemFock_IBS (V1.31): the route the composite density actually takes for a SALC basis ---
    virtual size_t  AODimension() const {return itsO.rows();}
    virtual void    AddAODensity(rsmat_t& Dao, const rsmat_t& D) const;   // Dao += O D O^T
    virtual rsmat_t MakeAOFock  (const rsmat_t& Dao, bool exchange) const;// the ONE whole-AO build
    virtual void    SliceAOFock (rsmat_t& Fab, const rsmat_t& Fao) const; // Fab += O^T Fao O

    // DFT 3-centre fitted Coulomb / Vxc.  The cached accessors Overlap3C/Repulsion3C are inherited
    // from Orbital_DFT_IBS unchanged: they key the *transformed* block under this irrep's
    // AngularID and, on a miss, call our MakeXxx3C below.  MakeXxx3C transforms the raw basis's
    // *cached* 3C (now safe -- the integral cache is re-entrant), so the raw 3C is computed once
    // and shared by every irrep.  Fit bases are atom-centred (geometry, not symmetry), so creation
    // delegates to the raw DFT basis.
    virtual rFIT_CD_ABS* CreateCDFitBasisSet (const Structure*, const qcMesh::MeshParams&) const;
    virtual rFIT_SF_ABS* CreateVxcFitBasisSet(const Structure*, const qcMesh::MeshParams&) const;
protected:
    virtual Projector3<double> MakeOverlap3C  (const rFIT_SF_ABS& c) const;
    virtual Projector3<double> MakeRepulsion3C(const rFIT_CD_ABS& c) const;
public:

    // The per-irrep cache distinction lives in BasisSetID() (raw id + "[label]"), so each irrep's
    // transformed blocks key separately.  (No RadialID/AngularID -- those are the atom identity face;
    // a SALC basis is molecular, and keys on the raw basis's own BasisSetID.)
    virtual std::string BasisSetID() const;
    virtual std::string Name()       const;

    // VectorFunction: a SALC evaluated at r is O^T applied to the raw AO values.
    virtual rvec_t     operator()(const rvec3_t&) const;
    virtual rvec3vec_t Gradient  (const rvec3_t&) const;

    virtual std::ostream& Write(std::ostream&) const;

private:
    //! The pair route's AO Fock: the PARTNER folds its own block up (through the abstract whole-system
    //! face), we build once.  V1.10 -- no reach into a sibling's private SALC columns.
    rsmat_t      PairAOFock(const rsmat_t& Dcd, const Orbital_HF_IBS<double>* bs_cd, bool exchange) const;
    rsmat_t      Transform(const rsmat_t& Mraw) const;          // O^T Mraw O, symmetrized
    Projector3<double> Transform3C(const Projector3<double>& raw) const; // Transform each fit-function matrix

    const Orbital_1E_IBS<double>*  itsRaw;            // raw whole-molecule AO basis (not owned)
    const Orbital_HF_IBS<double>*  itsRawHF;          // same object, HF interface (for the AO J/K build)
    const Orbital_DFT_IBS<double>* itsRawDFT;         // same object, DFT interface (raw 3C + fit bases)
    rmat_t                        itsO;               // this irrep's SALC columns (nAO x dGamma)
    std::string                   itsLabel;           // Mulliken irrep label
};

} //namespace
