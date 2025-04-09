// File: Angular_Sym.H Common interface for various atomic (spherical) symmetries.
#ifndef _Angular_Sym_H_
#define _Angular_Sym_H_

#include <Symmetry.H>
#include <Spin.H>   

class ElCounts;
struct ElCounts_l
{
    int n;  //# valance electrons for a given l
    int nu; //# of un paired electrons (for a given l)
    int GetN(Spin s) {return s==Spin::Up ? (n+nu)/2 : (n-nu)/2;  }
};

//---------------------------------------------------------------------------------
//
// Angular_Sym for atoms
//
class Angular_Sym
    : public virtual Symmetry
{
public:
    virtual ElCounts_l GetN(const ElCounts&) const {return {0,0};}

    virtual int GetL() const=0;
    virtual Angular_Sym* Clone() const=0;
};

#endif //_Angular_Sym_H_