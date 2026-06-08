// File: Symmetry/Yl.C  Non magnetic (m-degenerate) spherical harmonic Y_l(theta,phi) symmetry
module;
#include <iosfwd>
#include <cassert>

export module qchem.Symmetry.Yl;
export import qchem.Symmetry.Angular;


export class Yl_Sym
    : public virtual SphericalSym
{
public:
    Yl_Sym(size_t l) : itsL(l)
{};

    virtual size_t SequenceIndex() const {return itsL;} //Used for op<
    virtual size_t GetDegeneracy() const {return 2*itsL+1;}
    virtual size_t Getl         () const {return itsL;}
    virtual ivec_t Getmls       () const {return {};}
    virtual std::ostream&  Write(std::ostream&) const;
   
protected:
    size_t itsL;
    static const size_t LMAX=4;
};


