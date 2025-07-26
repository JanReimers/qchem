// File: Common/stl_io.C
module;
#include <iostream>
#include <vector>
#include <set>
#include <memory>
export module qchem.stl_io;

export 
{
//
//  W/R size and data
//
template <template<class> class V,class T> std::ostream& Write(std::ostream& os,const V<T>& v)
{
    // os << v.size() << " ";
    for (auto i:v) os << i << " ";
    return os;
}
template <template<class> class V,class T> std::ostream& Write(std::ostream& os,const V<std::unique_ptr<T>>& v)
{
    for (auto& i:v) os << *i.get() << " ";
    return os;
}

//
//  Generic streaming operators.
//
template <class T> std::ostream& operator<<(std::ostream& os,const std::vector<T>& v) {return Write(os,v);}
template <class T> std::ostream& operator<<(std::ostream& os,const std::set   <T>& s) {return Write(os,s);}

} //export block