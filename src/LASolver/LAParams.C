#include <LASolver/LAParams.H>
#include <iostream>

std::string PkgStrs[]={"OML","Lapack"};
std::string OrthStrs[]={"Cholsky","Eigen","SVD"};
std::ostream& operator<<(std::ostream& os, const LAParams& lap)
{
    os << PkgStrs[lap.LinearAlgebraPackage] << " " << OrthStrs[lap.BasisOrthoAlgorithm] << " truncate below " << lap.TruncationTolerance << " abstol=" << lap.abstol;
    return os;
}
 
