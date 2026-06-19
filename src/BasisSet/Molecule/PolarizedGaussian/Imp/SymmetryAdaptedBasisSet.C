// File: BasisSet/Molecule/PolarizedGaussian/Imp/SymmetryAdaptedBasisSet.C
module;
#include <string>
module qchem.BasisSet.Molecule.SymmetryAdaptedBasisSet;
import qchem.Symmetry.MolecularIrrep;     // MolecularIrrep (the labelled irrep symmetry)
import qchem.Blaze;                        // submatrix

namespace BasisSet::Molecule
{

SymmetryAdaptedBasisSet::SymmetryAdaptedBasisSet(const ::BasisSet::Orbital_1E_IBS<double>* raw,
                                                 const Symmetry::SALCs& salc)
{
    size_t nAO = salc.O.rows();
    for (size_t r=0; r+1<salc.blockStart.size(); ++r)        // one non-empty irrep block -> one IBS
    {
        size_t start = salc.blockStart[r], dG = salc.blockStart[r+1]-start;
        if (dG==0) continue;                                 // skip irreps absent from this basis
        rmat_t            Or    = blazem::submatrix(salc.O, 0, start, nAO, dG);
        const std::string label = salc.irrep[start];         // all columns of the block share it
        sym_t             sym(new MolecularIrrep(label, r));
        this->Insert(new ::BasisSet::SymmetryAdapted_IBS(raw, Or, label, sym));
    }
}

} //namespace
