module;
#include <vector>
export module qchem.BasisSet.Internal.ERI3;
export import qchem.Types;

export template <class T> using ERI3=std::vector<smat_t<T>>;
