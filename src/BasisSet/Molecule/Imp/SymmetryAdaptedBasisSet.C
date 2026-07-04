// File: BasisSet/Molecule/Imp/SymmetryAdaptedBasisSet.C
// Molecule-general: builds one SymmetryAdapted_IBS per irrep block from a precomputed SALC transform.
// No basis-set coupling -- the basis-specific SymmetryAdapt(rawBasis,cl) factory lives in each tree.
module;
#include <string>
#include <memory>
module qchem.BasisSet.Molecule.SymmetryAdaptedBasisSet;
import qchem.Symmetry.Molecule.Irrep;     // Molecule::Irrep (the labelled point-group irrep symmetry)
import qchem.Blaze;                        // submatrix

namespace qchem::BasisSet::Molecule
{

SymmetryAdaptedBasisSet::SymmetryAdaptedBasisSet(const ::qchem::BasisSet::Orbital_1E_IBS<double>* raw,
                                                 const Symmetry::Molecule::SALCs& salc)
{
    size_t nAO = salc.O.rows();
    auto cache = std::make_shared<::qchem::BasisSet::SymFockCache>();  // shared by all irreps (N^2 -> N)
    for (size_t r=0; r+1<salc.blockStart.size(); ++r)        // one non-empty irrep block -> one IBS
    {
        size_t start = salc.blockStart[r], dG = salc.blockStart[r+1]-start;
        if (dG==0) continue;                                 // skip irreps absent from this basis
        rmat_t            Or    = blazem::submatrix(salc.O, 0, start, nAO, dG);
        const std::string label = salc.irrep[start];         // all columns of the block share it
        sym_t             sym(new Symmetry::Molecule::Irrep(label, r));
        this->Insert(new ::qchem::BasisSet::SymmetryAdapted_IBS(raw, Or, label, sym, cache));
    }
}

} //namespace
