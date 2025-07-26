// File: ContractedGaussianH3.C   Class for managing a contraction of hermites3's.
module;
#include <cassert>
#include <vector>
module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Radial.ContractedGaussianH3;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;
import oml;

namespace PolarizedGaussian
{
//----------------------------------------------------------------------------------------
//
//  Construction zone.
//

ContractedGaussianH3::ContractedGaussianH3(const Vector<double>& c)
    : TheCoeff(c)
{};

ContractedGaussianH3::~ContractedGaussianH3()
{};

double ContractedGaussianH3::operator()(const Polarization& Pa,const Polarization& Pb,const Polarization& Pc) const
{
    //No UT coverage
    assert(itsH3s.size()==TheCoeff.size());
    double ret=0;
    int i=1;
    for (auto& b:itsH3s) 
    {
        ret+=(*b)(Pa,Pb,Pc)*TheCoeff(i);
        i++;
    }
    return ret;
}

} //namespace PolarizedGaussian
