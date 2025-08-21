module;
#include <vector>
export module qchem.BasisSet.Internal.ERI3;
// import oml;
export import qchem.Types;
import oml.SMatrix;

export template <class T> using ERI3=std::vector<SMatrix<T>>;
