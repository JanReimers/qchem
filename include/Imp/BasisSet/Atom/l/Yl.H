// File: Yl_Sym.H  Non magnetic (m-degenerate) spherical harmonic Y_l(theta,phi) symmetry
#ifndef _Yl_Sym_H_
#define _Yl_Sym_H_

#include "Imp/BasisSet/Atom/Angular.H"

//---------------------------------------------------------------------------------
//
//  The relevant quantum number for spherical symmetry is the magnitude of the
//  orbital angular momentum L.
//

class Yl_Sym
    : public virtual Angular_Sym
{
public:
    Yl_Sym(        );
    Yl_Sym(int theL);

    virtual size_t SequenceIndex() const; //Used for op<
    virtual int GetDegeneracy() const;
    
    virtual std::ostream&  Write(std::ostream&) const;
    virtual Angular_Sym* Clone(        ) const;

    virtual int     GetL() const {return itsL;}
    virtual ElCounts_l GetN(const ElCounts&) const;
protected:
    int itsL;
    static const int LMax=3;
};

#endif //_Yl_Sym_H_
