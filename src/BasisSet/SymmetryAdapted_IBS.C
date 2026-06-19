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
export module qchem.BasisSet.SymmetryAdapted_IBS;
export import qchem.BasisSet.Orbital_1E_IBS;
import qchem.BasisSet.Internal.IrrepBasisSetImp;
import qchem.Cluster;
import qchem.Types;

export namespace BasisSet
{

class SymmetryAdapted_IBS
    : public virtual Orbital_1E_IBS<double>
    , public IrrepBasisSetImp<double>        // provides GetSymmetry / GetSymt / GetIrrep
{
public:
    // raw: the whole-molecule AO basis (not owned).  Oblock: this irrep's SALC columns
    // (nAO x dGamma).  label: Mulliken irrep label (used for the cache key + display).
    SymmetryAdapted_IBS(const Orbital_1E_IBS<double>* raw, const rmat_t& Oblock,
                        const std::string& label, const sym_t& sym);

    virtual size_t GetNumFunctions() const {return itsO.columns();}

    // 1-electron integrals in the irrep basis (O^T M_raw O).
    virtual rsmat_t MakeOverlap()                 const;
    virtual rsmat_t MakeKinetic()                 const;
    virtual rsmat_t MakeNuclear(const Cluster* cl) const;

    // Distinct IDs so the global cache keys the transformed blocks per irrep.
    virtual std::string RadialID()  const;
    virtual std::string AngularID() const;
    virtual std::string Name()      const;

    // VectorFunction: a SALC evaluated at r is O^T applied to the raw AO values.
    virtual rvec_t     operator()(const rvec3_t&) const;
    virtual rvec3vec_t Gradient  (const rvec3_t&) const;

    virtual std::ostream& Write(std::ostream&) const;

private:
    rsmat_t Transform(const rsmat_t& Mraw) const;     // O^T Mraw O, symmetrized

    const Orbital_1E_IBS<double>* itsRaw;             // raw whole-molecule AO basis (not owned)
    rmat_t                        itsO;               // this irrep's SALC columns (nAO x dGamma)
    std::string                   itsLabel;           // Mulliken irrep label
};

} //namespace
