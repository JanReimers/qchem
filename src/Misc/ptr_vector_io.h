// File: ptrvector_io  io templates for streaming.

#include "oml/imp/binio.h"
#include "oml/imp/stream.h"
#include "Misc/ptr_vector.h"

template <class T> std::ostream& operator<<(std::ostream& os, const ptr_vector<T*>& v)
{
    unsigned n=v.size();
    if(StreamableObject::Binary())
        BinaryWrite(n,os);
    else
        os << n << std::endl;
    typedef typename optr_vector<T*>::const_iterator CITER;
    for (CITER i=v.begin(); i!=v.end(); i++) os << *i;
    return os;
}

template <class T> std::istream& operator>>(std::istream& is, ptr_vector<T*>& v)
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

template <class T> std::ostream& operator<<(std::ostream& os, const optr_vector<T*>& v)
{
    unsigned n=v.size();
    if(StreamableObject::Binary())
        BinaryWrite(n,os);
    else
        os << n << std::endl;
    typedef typename optr_vector<T*>::const_iterator CITER;
    for (CITER i=v.begin(); i!=v.end(); i++) os << *i;
    return os;
}

template <class T> std::istream& operator>>(std::istream& is, optr_vector<T*>& v)
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

