// File: Common/Convert.C
module;
#include <vector>
#include <valarray>
export module qchem.Conversions;
import oml.Vector;

export 
{

template <class T> std::valarray<T> to_valarray(const std::vector<T>& v)
{
    std::valarray<T> ret(v.size());
    size_t n=0;
    for (auto iv:v) ret[n++]=iv;
    return ret;
}

template <class T> Vector<T> to_omlVector(const std::valarray<T>& v)
{
    Vector<T> ret(v.size());
    size_t n=0;
    for (auto i:v) ret(++n)=i;
    return ret;
}


}

