// File: ContractedGaussianH3.H  Class for managing a contraction of hermites3's.
#ifndef _ContractedGaussianH3_H_
#define _ContractedGaussianH3_H_


#include "types.H"
#include "Imp/BasisSet/Molecule/PolarizedGaussian/MnD/Hermite3.H"
#include "Imp/Containers/ptr_vector.h"

namespace PolarizedGaussian
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
        itsH3s.push_back(h3);
    }

    optr_vector1<Hermite3*> itsH3s;
    const Vector<double>& TheCoeff; //Contraction coefficients.
};

} //namespace PolarizedGaussian

#endif
