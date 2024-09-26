// File: ContractedGaussianH3.H  Class for managing a contraction of hermites3's.
#ifndef _ContractedGaussianH3_H_
#define _ContractedGaussianH3_H_



#include "BasisSetImplementation/PolarizedGaussian/Hermite/Hermite3.H"
#include "Misc/ptr_vector.h"

template <class T> class Vector;
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

    ptr_vector <Hermite3*>  itsH3s;
    const Vector<double>& TheCoeff; //Contraction coefficients.
};
#endif
