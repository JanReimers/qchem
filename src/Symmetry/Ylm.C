// File: Slater_m/QuantumNumer.H  Spherical harmonic Ylm symmetry.
module;
#include <vector>
#include <iosfwd>

export module qchem.Symmetry.Ylm;
export import qchem.Symmetry.Yl;

//---------------------------------------------------------------------------------
//
//  angular momentum l and ml.
//

export class Ylm_Sym
    : public virtual SphericalSym
    , private Yl_Sym
{
public:
    Ylm_Sym(size_t l, const ivec_t& mls);

    virtual size_t SequenceIndex() const; //Used for op<
    virtual size_t GetDegeneracy() const;
    using Yl_Sym::Getl;
    virtual ivec_t Getmls() const {return mls;}

    virtual std::ostream&  Write(std::ostream&) const;
   
protected:
    ivec_t mls;
};

