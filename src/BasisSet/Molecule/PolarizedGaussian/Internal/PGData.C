// File BasisSet/Molecule/PolarizedGaussian/Internal/PGData.C Flattened rep of PG IBS suitable for integral evaluation.
module;
#include <vector>
#include <string>

export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.PGData;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Block;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.GaussianRF;
import qchem.Types;

export namespace BasisSet::Molecule::PolarizedGaussian
{

struct PGData
{
    std::string RadialID () const;
    std::string AngularID() const;
    std::string BasisSetID() const; // geometry-aware cache identity: radial @ centre : pol per fn

    void Init(std::vector<const Block*>&);
    std::vector<const GaussianRF*> radials; // Flattened radials
    std::vector<Polarization>          pols;    // Flattened polarizations
    rvec_t                             ns;      //Norm constants

    size_t size() const {return radials.size();}

};

}