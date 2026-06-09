// File: Irrep.C  Combine Symmetry with Spin.
module;
#include <iosfwd>
export module qchem.Symmetry.Irrep;
export import qchem.Symmetry;
export import qchem.Symmetry.Spin;
import qchem.Streamable;

export struct Irrep
    : public virtual Streamable
{   
    Irrep() : ms(Spin::None), sym(0) {};
    Irrep(Spin _ms,const sym_t& _sym);
    ~Irrep();

    virtual size_t  SequenceIndex() const; //Used for op<

    friend bool operator<(const Irrep& a, const Irrep& b)
    {
        return a.SequenceIndex()<b.SequenceIndex();
    }
    virtual size_t GetDegeneracy() const;

    std::ostream& Write(std::ostream&) const;
    
    Spin  ms;
    sym_t sym;

    static const size_t ms_max; //Used for calculating sequenxe indexes.
};


