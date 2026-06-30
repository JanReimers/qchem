// File: BasisSet/Molecule/Evaluators/PG_Cart_MnD/PGData.C  Flattened rep of PG IBS suitable for integral evaluation.
module;
#include <vector>
#include <string>

export module qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Internal.Block;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.Types;

export namespace qchem::BasisSet::Molecule::Evaluators::PG_Cart_MnD
{

struct PGData
{
    std::string BasisSetID() const; // geometry-aware cache identity: radial @ centre : pol per fn

    void Init(std::vector<const Block*>&);
    std::vector<const GaussianRF*> radials; // Flattened radials
    std::vector<Polarization>          pols;    // Flattened polarizations
    rvec_t                             ns;      //Norm constants

    size_t size() const {return radials.size();}

};

}