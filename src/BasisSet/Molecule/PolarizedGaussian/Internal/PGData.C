// File BasisSet/Molecule/PolarizedGaussian/Internal/PGData.C Flattened rep of PG IBS suitable for integral evaluation.
module;
#include <vector>
export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.PGData;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Block;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.IEClient;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.RadialFunction;
import qchem.Types;

export namespace PolarizedGaussian
{

struct PGData
{
    void Init(std::vector<const Block*>&);
    std::vector<const RadialFunction*> radials1; // Flattened radials
    std::vector<Polarization>          pols1;    // Flattened polarizations
    rvec_t                             ns1;      //Norm constants

    size_t size1() const {return radials1.size();}

};

}