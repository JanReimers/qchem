// File: ContractedGaussianH3.C   Class for managing a contraction of hermites3's.



#include "BasisSetImplementation/PolarizedGaussian/ContractedGaussian/ContractedGaussianH3.H"
#include "oml/vector.h"
#include <cassert>

//----------------------------------------------------------------------------------------
//
//  Construction zone.
//

ContractedGaussianH3::ContractedGaussianH3(const Vector<double>& c)
    : TheCoeff(c)
{};

ContractedGaussianH3::~ContractedGaussianH3()
{
    ptr_vector <Hermite3*>::iterator b(itsH3s.begin());
    for (; b!=itsH3s.end(); b++) delete &b;
}

double ContractedGaussianH3::operator()(const Polarization& Pa,const Polarization& Pb,const Polarization& Pc) const
{
    assert(itsH3s.size()==TheCoeff.size());
    double ret=0;
    ptr_vector <Hermite3*>::const_iterator b(itsH3s.begin());
    Vector<double>::const_iterator c(TheCoeff.begin());

    for (; b!=itsH3s.end()&&c!=TheCoeff.end(); b++,c++) ret+=(*b)(Pa,Pb,Pc)*(*c);
    return ret;
}
