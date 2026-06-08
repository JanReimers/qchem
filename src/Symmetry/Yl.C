// File: Symmetry/Yl.C  Non magnetic (m-degenerate) spherical harmonic Y_l(theta,phi) symmetry
module;
#include <iosfwd>
#include <cassert>

export module qchem.Symmetry.Yl;
export import qchem.Symmetry.Angular;
import qchem.Symmetry.AtomEC;


export class Yl_Sym
    : public virtual Angular_Sym
{
public:
    Yl_Sym(int theL);

    virtual size_t SequenceIndex() const; //Used for op<
    virtual int GetDegeneracy() const;
    
    virtual std::ostream&  Write(std::ostream&) const;
   
    virtual int     GetL() const {return itsL;}
protected:
    int itsL;
    static const int LMAX=4;
};


