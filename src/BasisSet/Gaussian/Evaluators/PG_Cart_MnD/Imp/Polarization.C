// File: Polarization.C  Structure describing just the polarization portion of a basis function.
module;
#include <iostream>
module qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.Polarization;

namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD
{

int Polarization::GetMaximumL() const
{
    int t=l>m ? l : m;
    return t>n ? t : n;
}



std::ostream& operator<<(std::ostream& os, const Polarization& p)
{
    static char SPDF[]="SPDFGHIJ";
    int n=1;
    os << SPDF[p.GetTotalL()];
    for (int i=0; i<p.n; i++,n++) os << "x";
    for (int i=0; i<p.l; i++,n++) os << "y";
    for (int i=0; i<p.m; i++,n++) os << "z";
//    for (;n<=5;n++) os << " ";
    os << " ";
    return os;
}


} //namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD
