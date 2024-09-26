// File: RadialFunctionReader.h  Abstract interface for a basis set reader.
#ifndef _RadialFunctionReader_h_
#define _RadialFunctionReader_h_

#include <vector>

class Atom;
class RadialFunction;

//-------------------------------------------------------------------------
//
//  Derived classes all read particular formats for the basis sets.
//
class RadialFunctionReader
{
public:
    virtual ~RadialFunctionReader() {};

    virtual RadialFunction*  ReadNext(const Atom&)      =0;
    virtual bool             FindAtom(const Atom&)      =0;
    virtual std::vector<int> GetLs   (           ) const=0;
};

#endif
