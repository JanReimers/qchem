#ifndef STL_IO_H_INCLUDED
#define STL_IO_H_INCLUDED

#include "oml/imp/stream.h"
#include "oml/imp/binio.h"
#include <iostream>
#include <typeinfo>
#include <cassert>
#include <vector>
#include <set>
#include <memory>

// Make some stl containers streamable in Pretty, Ascii and Binary modes.

//
//  W/R headers with output mode and typename
//
void WriteHeader(std::ostream& os,c_str type);
StreamableObject::Mode ReadHeader(std::istream& is,c_str type);

template <class T> void WriteHeader(const T& t,std::ostream& os)
{
    WriteHeader(os, typeid(T).name());
}

template <class T> StreamableObject::Mode ReadHeader(T& t,std::istream& is)
{
    return ReadHeader(is,typeid(T).name());
}

//
//  W/R size and data
//
template <template<class> class V,class T> std::ostream& Write(std::ostream& os,const V<T>& v)
{
    WriteHeader(v,os);
    if (StreamableObject::Binary()) 
    {
        BinaryWrite(v.size(),os);
        for (auto i:v) BinaryWrite(i,os);
    }
    else 
    {   
        os << v.size() << " ";
        for (auto i:v) os << i << " ";
    }
    return os;
}
template <template<class> class V,class T> std::ostream& Write(std::ostream& os,const V<std::unique_ptr<T>>& v)
{
    WriteHeader(v,os);
    if (StreamableObject::Binary()) 
    {
        BinaryWrite(v.size(),os);
        for (auto& i:v) BinaryWrite(*i.get(),os);
    }
    else 
    {   
        os << v.size() << " ";
        for (auto& i:v) os << *i.get() << " ";
    }
    return os;
}

template <template<class> class V,class T> std::istream& Read(std::istream& is,V<T>& v);

template <class T> std::istream& Read(std::istream& is,std::vector<T>& v)
{
    StreamableObject::Mode current=ReadHeader(v,is);
    size_t n;
    if (StreamableObject::Binary()) BinaryRead(n,is);
    if (StreamableObject::Ascii())  is >> n;;
    
    assert(is);
    if (v.size()==0)
        v.resize(n);
    else
        assert(v.size()==n);
    
    if (StreamableObject::Binary()) for (auto& i:v) BinaryRead(i,is);
    if (StreamableObject::Ascii())  for (auto& i:v) is >> i;

    StreamableObject::SetOutputMode(current); //Restore to previous state
    return is;
}

// Specialization for std::set
template <class T> std::istream& Read(std::istream& is,std::set<T>& v)
{
    StreamableObject::Mode current=ReadHeader(v,is);
    std::size_t n;
    if (StreamableObject::Binary()) BinaryRead(n,is);
    if (StreamableObject::Ascii())  is >> n;;
    assert(is);
    
    if (v.size()!=0)
        v.clear();
    T vi;
    if (StreamableObject::Binary())
        for (std::size_t i=0;i<n;i++)
        {
            BinaryRead(vi,is);
            v.insert(vi);
        } 
    if (StreamableObject::Ascii())
        for (std::size_t i=0;i<n;i++)
        {
            is >> vi;
            v.insert(vi);
        } 
        
    StreamableObject::SetOutputMode(current); //Restore to previous state
    return is;
}

//
//  Generic streaming operators.
//
template <class T> std::ostream& operator<<(std::ostream& os,const std::vector<T>& v) {return Write(os,v);}
template <class T> std::ostream& operator<<(std::ostream& os,const std::set   <T>& s) {return Write(os,s);}
template <class T> std::istream& operator>>(std::istream& is,      std::vector<T>& v) {return Read (is,v);}
template <class T> std::istream& operator>>(std::istream& is,      std::set   <T>& s) {return Read (is,s);}




#endif // STL_IO_H_INCLUDED
