// File: Reader.h  Abstract interface for a basis set reader.
#ifndef _Reader_h_
#define _Reader_h_

#include "types.H"
#include <vector>

namespace PolarizedGaussian
{

class RadialFunction;
//-------------------------------------------------------------------------
//
//  Derived classes all read particular formats for the basis sets.
//
class Reader
{
public:
    virtual ~Reader() {};

    virtual RadialFunction*  ReadNext(const Atom&)      =0;
    virtual bool             FindAtom(const Atom&)      =0;
    virtual std::vector<int> GetLs   (           ) const=0;
};

} //namespace PolarizedGaussian

#endif
