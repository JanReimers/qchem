// File: BasisSet/Molecule/PolarizedGaussian/SymmetryAdaptedBasisSet.C
// A molecular basis re-expressed in symmetry-adapted (SALC) blocks: one SymmetryAdapted_IBS
// per point-group irrep, each carrying its Mulliken label.  The SCF iterator / accelerators
// iterate these IBSs exactly as they do the l-channels of an atom -- the molecular case is
// now "atoms with O != identity".  Stage 4 of the molecular-symmetry plan.
module;
export module qchem.BasisSet.Molecule.SymmetryAdaptedBasisSet;
export import qchem.BasisSet;                         // BasisSet<double>, Orbital_1E_IBS
import qchem.BasisSet.Internal.BasisSetImp;           // BasisSetImp (iteration/Insert)
import qchem.BasisSet.SymmetryAdapted_IBS;            // the per-irrep decorator
import qchem.Symmetry.SALC;                           // SALCs (the transform O + labels)

export namespace BasisSet::Molecule
{

class SymmetryAdaptedBasisSet
    : public virtual ::BasisSet::BasisSet<double>
    , public ::BasisSet::BasisSetImp<double>
{
public:
    // raw: the whole-molecule AO basis.  It is REFERENCED by the per-irrep decorators (not
    // owned here), so it must outlive this object.  salc: the SALC transform from BuildSALCs.
    SymmetryAdaptedBasisSet(const ::BasisSet::Orbital_1E_IBS<double>* raw, const Symmetry::SALCs& salc);
};

} //namespace
