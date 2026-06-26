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
export import qchem.BasisSet.Orbital_HF_IBS;         // HF 2-electron mixin + ERI4
export import qchem.BasisSet.Orbital_DFT_IBS;        // DFT 3-centre mixin (fitted Coulomb / Vxc)
import qchem.BasisSet.Internal.IrrepBasisSetImp;
import qchem.Structure;
import qchem.Types;

export namespace BasisSet
{

// Shared by one basis's irrep decorators: memoizes the AO Coulomb/exchange built from each
// cd-irrep density block (raw->AccumulateDirect/Exchange), so the O(nAO^4) AO build happens
// once per cd-irrep per SCF iteration -- not once per (irrep, cd-irrep) pair (was N^2/iter).
// Keyed by the cd-irrep IBS pointer; invalidated when that block's density changes.
class SymFockCache
{
public:
    const rsmat_t& Direct  (const Orbital_HF_IBS<double>* raw, const void* cd,
                            const rmat_t& Ocd, const rsmat_t& Dcd);
    const rsmat_t& Exchange(const Orbital_HF_IBS<double>* raw, const void* cd,
                            const rmat_t& Ocd, const rsmat_t& Dcd);
private:
    struct Entry { rsmat_t D, M; bool valid=false; };
    std::map<const void*, Entry> itsJ, itsK;
};

class SymmetryAdapted_IBS
    : public virtual Orbital_1E_IBS<double>
    , public virtual Orbital_HF_IBS<double>  // HF Coulomb/exchange (the 2-electron path)
    , public virtual Orbital_DFT_IBS<double> // DFT 3-centre fitted Coulomb / Vxc
    , public IrrepBasisSetImp<double>        // provides GetSymmetry / GetSymt / GetIrrep
{
public:
    // raw: the whole-molecule AO basis (not owned).  Oblock: this irrep's SALC columns
    // (nAO x dGamma).  label: Mulliken irrep label (used for the cache key + display).
    // cache: optional, shared across the irreps of one basis (turns the 2-e build N^2 -> N per
    // iteration).  Null -> build the AO Coulomb/exchange directly each call (fine for tests).
    SymmetryAdapted_IBS(const Orbital_1E_IBS<double>* raw, const rmat_t& Oblock,
                        const std::string& label, const sym_t& sym,
                        std::shared_ptr<SymFockCache> cache = nullptr);

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
    // Pure-virtual ERI accessors -- unused here (Accumulate* are overridden); never called.
    virtual ERI4 MakeDirect  (const Orbital_HF_IBS<double>&) const {return ERI4();}
    virtual ERI4 MakeExchange(const Orbital_HF_IBS<double>&) const {return ERI4();}

    // DFT 3-centre fitted Coulomb / Vxc.  The cached accessors Overlap3C/Repulsion3C are inherited
    // from Orbital_DFT_IBS unchanged: they key the *transformed* block under this irrep's
    // AngularID and, on a miss, call our MakeXxx3C below.  MakeXxx3C transforms the raw basis's
    // *cached* 3C (now safe -- the integral cache is re-entrant), so the raw 3C is computed once
    // and shared by every irrep.  Fit bases are atom-centred (geometry, not symmetry), so creation
    // delegates to the raw DFT basis.
    virtual Fit_IBS* CreateCDFitBasisSet (const Structure*, const qcMesh::MeshParams&) const;
    virtual Fit_IBS* CreateVxcFitBasisSet(const Structure*, const qcMesh::MeshParams&) const;
protected:
    virtual ERI3<double> MakeOverlap3C  (const Fit_IBS& c) const;
    virtual ERI3<double> MakeRepulsion3C(const Fit_IBS& c) const;
public:

    // RadialID/AngularID delegate to the raw basis unchanged; the per-irrep cache distinction
    // lives in BasisSetID() (raw id + "[label]"), so each irrep's transformed blocks key separately.
    virtual std::string RadialID()   const;
    virtual std::string AngularID()  const;
    virtual std::string BasisSetID() const;
    virtual std::string Name()       const;

    // VectorFunction: a SALC evaluated at r is O^T applied to the raw AO values.
    virtual rvec_t     operator()(const rvec3_t&) const;
    virtual rvec3vec_t Gradient  (const rvec3_t&) const;

    virtual std::ostream& Write(std::ostream&) const;

private:
    rsmat_t      Transform(const rsmat_t& Mraw) const;          // O^T Mraw O, symmetrized
    ERI3<double> TransformERI3(const ERI3<double>& raw) const;  // Transform each fit-function matrix

    const Orbital_1E_IBS<double>*  itsRaw;            // raw whole-molecule AO basis (not owned)
    const Orbital_HF_IBS<double>*  itsRawHF;          // same object, HF interface (for the AO J/K build)
    const Orbital_DFT_IBS<double>* itsRawDFT;         // same object, DFT interface (raw 3C + fit bases)
    rmat_t                        itsO;               // this irrep's SALC columns (nAO x dGamma)
    std::string                   itsLabel;           // Mulliken irrep label
    std::shared_ptr<SymFockCache> itsCache;           // shared AO J/K cache (null = build directly)
};

} //namespace
