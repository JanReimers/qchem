module;
#include <vector>
export module qchem.BasisSet1.Internal.ERI3;
export import qchem.Types;

export template <class T> using ERI3=std::vector<smat_t<T>>;
export template <class T> double fnorm(const ERI3<T>& a, const ERI3<T>& b);
export template <class T> double relative_fnorm(const ERI3<T>& a, const ERI3<T>& b);