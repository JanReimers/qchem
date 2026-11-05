// File: BasisSet/Gaussian/Point/Imp/SymmetryAdaptedBasisSet.C
// Molecule-general: builds one SymmetryAdapted_IBS per irrep block from a precomputed SALC transform.
// No basis-set coupling -- the basis-specific SymmetryAdapt(rawBasis,cl) factory lives in each tree.
module;
#include <vector>
#include <memory>
#include <string>
#include <memory>
module qchem.BasisSet.Gaussian.Point.SymmetryAdaptedBasisSet;
import qchem.Symmetry.Molecule.Irrep;     // Gaussian::Irrep (the labelled point-group irrep symmetry)
import qchem.Blaze;                        // submatrix
import qchem.BasisSet.Internal.SymmetryAdapted_IBS;

namespace qchem::BasisSet::Gaussian
{

SymmetryAdaptedBasisSet::SymmetryAdaptedBasisSet(const ::qchem::BasisSet::Orbital_1E_IBS<double>* raw,
                                                 const Symmetry::Molecule::SALCs& salc)
{
    size_t nAO = salc.O.rows();
    const auto group = std::make_shared<const std::vector<Symmetry::SymOp>>(salc.ops);   // ONE group, shared by every block
    for (size_t r=0; r+1<salc.blockStart.size(); ++r)        // one non-empty irrep block -> one IBS
    {
        size_t start = salc.blockStart[r], dG = salc.blockStart[r+1]-start;
        if (dG==0) continue;                                 // skip irreps absent from this basis
        rmat_t            Or    = blazem::submatrix(salc.O, 0, start, nAO, dG);
        const std::string label = salc.irrep[start];         // all columns of the block share it
        sym_t             sym(new Symmetry::Molecule::Irrep(label, r));
        this->Insert(new ::qchem::BasisSet::SymmetryAdapted_IBS(raw, Or, label, sym, group));
    }
}

} //namespace
