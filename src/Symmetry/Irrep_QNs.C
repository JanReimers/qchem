// File: Irrep_QNs.C  Combine Symmetry with Spin.
module;
#include <memory>

export module qchem.Symmetry.Irrep;
export import qchem.Symmetry;
export import qchem.Symmetry.Spin;
import qchem.Streamable;

export struct Irrep_QNs
    : public virtual Streamable
{   
    typedef std::shared_ptr<const Symmetry> sym_t;
    Irrep_QNs() : ms(Spin::None), sym(0) {};
    Irrep_QNs(Spin _ms,const sym_t& _sym);
    ~Irrep_QNs();

    virtual size_t  SequenceIndex() const; //Used for op<

    friend bool operator<(const Irrep_QNs& a, const Irrep_QNs& b)
    {
        return a.SequenceIndex()<b.SequenceIndex();
    }
    virtual size_t GetDegeneracy() const;

    std::ostream& Write(std::ostream&) const;
    
    Spin  ms;
    sym_t sym;

    static const size_t ms_max; //Used for calculating sequenxe indexes.
};
