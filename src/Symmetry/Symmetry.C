// File: Symmetry.C  Abstract interface for symmetries that do not include spin.
module;
#include <string>
#include <memory>

export module qchem.Symmetry;
import qchem.Streamable;
export import qchem.Types;


export namespace qchem::Symmetry
{
class Symmetry
    : public virtual Streamable
{
public:
    virtual ~Symmetry() {};
    virtual size_t SequenceIndex() const=0; //Used for op<
    //! Does not include spin degeneracy which is handled separately
    virtual size_t GetDegeneracy     () const=0;
    virtual size_t GetPrincipleOffset() const=0; //Add to principle QN.  For atoms this is just l.
    //! \brief Brillouin-zone integration weight of this symmetry's block (\f$\sum w=1\f$ over a BZ
    //! sampling).  1 for ordinary symmetries; a Bloch k-point carries its k-mesh weight \f$w_k\f$ so a
    //! charge density built per-irrep sums to the BZ average \f$\sum_k w_k\rho_k\f$.
    virtual double GetWeight         () const {return 1.0;}
    std::string    GetLabel          () const;
};

} //namespace
// use this type to pass around polymorkpic Symmetry instances.
export using sym_t=std::shared_ptr<const qchem::Symmetry::Symmetry>;
