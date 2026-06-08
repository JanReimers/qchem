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
    : public virtual Symmetry
    , public Yl_Sym
{
public:
    Ylm_Sym(size_t l, const std::vector<int>& ml);

    virtual size_t SequenceIndex() const; //Used for op<
    virtual size_t GetDegeneracy() const;
    using Yl_Sym::GetL;
    const std::vector<int>& Getmls() const {return ml;}

    virtual std::ostream&  Write(std::ostream&) const;
   
protected:
    std::vector<int> ml;
};

