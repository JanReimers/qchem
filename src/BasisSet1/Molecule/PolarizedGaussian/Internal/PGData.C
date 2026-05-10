// File BasisSet/Molecule/PolarizedGaussian/Internal/PGData.C Flattened rep of PG IBS suitable for integral evaluation.
module;
#include <vector>
#include <string>

export module qchem.BasisSet1.Molecule.PolarizedGaussian.Internal.PGData;
import qchem.BasisSet1.Molecule.PolarizedGaussian.Internal.Block;
import qchem.BasisSet1.Molecule.PolarizedGaussian.Internal.Polarization;
import qchem.BasisSet1.Molecule.PolarizedGaussian.Internal.RadialFunction;
import qchem.Types;

export namespace BasisSet1::Molecule::PolarizedGaussian
{

struct PGData
{
    std::string RadialID () const;
    std::string AngularID() const;

    void Init(std::vector<const Block*>&);
    std::vector<const RadialFunction*> radials; // Flattened radials
    std::vector<Polarization>          pols;    // Flattened polarizations
    rvec_t                             ns;      //Norm constants

    size_t size() const {return radials.size();}

};

}