// File: Symmetry/Ylm_Sym.H  Magnetic spherical harmonic Y_lm(theta,phi) symmetry
#ifndef _Symmetry_Ylm_Sym_H_
#define _Symmetry_Ylm_Sym_H_

#include "Imp/BasisSet/Atom/l/Yl.H"
#include <vector>


//---------------------------------------------------------------------------------
//
//  angular momentum l and ml.
//

class Ylm_Sym
    : public virtual Symmetry
    , public Yl_Sym
{
public:
    Ylm_Sym(             );
    Ylm_Sym(int l, const std::vector<int>& ml);

    virtual size_t     SequenceIndex() const; //Used for op<
    virtual int        GetDegeneracy() const;
    virtual ElCounts_l GetN(const ElCounts&) const;

    virtual std::ostream&  Write(std::ostream&) const;
    virtual Angular_Sym*   Clone(        ) const;

    using Yl_Sym::GetL;
protected:
    std::vector<int> ml;
};

#endif //_Symmetry_Ylm_Sym_H_
