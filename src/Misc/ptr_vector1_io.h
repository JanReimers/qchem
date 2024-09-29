// File: ptr_vector1_io  io templates for streaming.
#ifndef _ptr_vector1_io_h
#define _ptr_vector1_io_h

#include "oml/imp/binio.h"
#include "oml/imp/stream.h"
#include "Misc/ptr_vector1.h"

template <class T> std::ostream& operator<<(std::ostream& os, const std::vector<T*>& v)
{
    unsigned n=v.size();
    if(StreamableObject::Binary())
        BinaryWrite(n,os);
    if(StreamableObject::Ascii())
        os << n << std::endl;
        
    for (auto i:v) os << *i;
    return os;
}

template <class T> std::istream& operator>>(std::istream& is, std::vector<T*>& v)
{
    v.clear();
    unsigned n;
    if(StreamableObject::Binary())
        BinaryRead(n,is);
    else
        is >> n;

    for(unsigned i=0; i<n; i++)
    {
        v.push_back(T::Factory(is));
        is >> *v.back();
    }
    return is;
}

#endif // _ptr_vector1_io_h
