// File: ContractedGaussianH3.H  Class for managing a contraction of hermites3's.
module;

#include <vector>
#include <memory>
export module qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Radial.ContractedGaussianH3;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.MnD.Hermite3;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Polarization;

import oml;

export namespace PolarizedGaussian
{

//
// Store a list hermite3 pointers.  Contraction coefficients should be absorbed into the scale factors
// of the uncontracted hermite3's before insertion.
//
class ContractedGaussianH3 : public Hermite3
{
public:
    ContractedGaussianH3(const Vector<double>&);
    virtual ~ContractedGaussianH3();

    virtual double operator()(const Polarization& Pa,const Polarization& Pb,const Polarization& Pc) const;

private:
    friend class ContractedGaussianRF;
    void Insert(Hermite3* h3)
    {
        itsH3s.push_back(std::unique_ptr<Hermite3>(h3));
    }

    typedef std::vector<std::unique_ptr<Hermite3>> h3v_t;
    h3v_t itsH3s;
    const Vector<double>& TheCoeff; //Contraction coefficients.
};

} //namespace PolarizedGaussian
