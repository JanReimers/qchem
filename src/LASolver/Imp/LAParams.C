// File: Imp/LAParams.C Parameter for determing how the linear algebra for Lowden orthogonalization is done.
module;
#include <iostream>
module qchem.LAParams;
std::string PkgStrs[]={"OML","Lapack"};
std::string OrthStrs[]={"Cholsky","Eigen","SVD"};
std::ostream& operator<<(std::ostream& os, const LAParams& lap)
{
    os << PkgStrs[lap.LinearAlgebraPackage] << " " << OrthStrs[lap.BasisOrthoAlgorithm] << " truncate below " << lap.TruncationTolerance << " abstol=" << lap.abstol;
    return os;
}
 
