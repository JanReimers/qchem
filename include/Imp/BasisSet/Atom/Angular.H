// File: Angular_Sym.H Common interface for various atomic (spherical) symmetries.
#ifndef _Angular_Sym_H_
#define _Angular_Sym_H_

#include <Symmetry.H>

class ElCounts;
class ElCounts_l;

//---------------------------------------------------------------------------------
//
// Angular_Sym for atoms
//
class Angular_Sym
    : public virtual Symmetry
{
public:
    virtual int GetPrincipleOffset() const {return GetL();} //Add to principle QN.  For atoms this is just l.
    // New atom specific functions.
    virtual ElCounts_l GetN(const ElCounts&) const=0;
    virtual int GetL() const=0;
};

#endif //_Angular_Sym_H_