// File: Slater_m/QuantumNumer.H  Spherical harmonic Ylm symmetry.
module;
#include <vector>
#include <cstddef>
#include <iosfwd>

export module qchem.Symmetry.Ylm;
export import qchem.Symmetry.Yl;
import qchem.Symmetry.AtomEC;

//---------------------------------------------------------------------------------
//
//  angular momentum l and ml.
//

export class Ylm_Sym
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
   
    using Yl_Sym::GetL;
protected:
    std::vector<int> ml;
};

