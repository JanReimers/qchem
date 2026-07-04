// File: Symmetry.C  Abstract interface for a symmetry-group irrep label.  Most concretes are spatial (spin
// is added on top by Irrep), but the relativistic double group (SphericalSpinor, Ωκ) is spin-orbit coupled --
// see CarriesSpin().
module;
#include <string>
#include <memory>

export module qchem.Symmetry;
import qchem.Streamable;
export import qchem.Types;


export namespace qchem::Symmetry
{
//! \brief Abstract SPATIAL symmetry label of an orbital/irrep -- the base of the quantum-number hierarchy
//! (spin is added on top by \c Irrep).  Pure quantum-number identity: the concretes (\c Molecule::Irrep for
//! point-group irreps, \c Yl / \c Ylm and \c Ωκ / \c Ωκmj for atomic/Dirac shells, \c BlochQN for crystal
//! k-points, \c UnitQN for no symmetry) carry only their labelling data.  The class does no physics -- it
//! provides an ORDERING key (\c SequenceIndex, for sorting/caching irreps) plus block metadata (degeneracy,
//! BZ weight).  Passed around polymorphically as \c sym_t.
class Symmetry
    : public virtual Streamable
{
public:
    virtual ~Symmetry() {};
    virtual size_t SequenceIndex() const=0;         //!< monotone ordering key (defines \c operator< on irreps)
    //! Spatial degeneracy of the block.  Does NOT include spin degeneracy -- that is handled by \c Irrep/\c Spin.
    virtual size_t GetDegeneracy     () const=0;
    virtual size_t GetPrincipleOffset() const=0;    //!< offset added to the principal QN (for atoms this is \a l)
    //! \brief Brillouin-zone integration weight of this symmetry's block (\f$\sum w=1\f$ over a BZ
    //! sampling).  1 for ordinary symmetries; a Bloch k-point carries its k-mesh weight \f$w_k\f$ so a
    //! charge density built per-irrep sums to the BZ average \f$\sum_k w_k\rho_k\f$.
    virtual double GetWeight         () const {return 1.0;}
    //! \brief Does this irrep label already carry spin (spin-orbit coupled)?  \c false for ordinary spatial
    //! symmetries -- spin is layered on top by \c Irrep, whose \a ms is the non-relativistic spin channel.
    //! \c true only for the relativistic double group (\c SphericalSpinor), where \f$\kappa\f$ encodes
    //! \f$j=l\pm\tfrac12\f$; there the layer-2 \a ms is \c Spin::None (spin is already in \f$\kappa\f$).
    virtual bool   CarriesSpin        () const {return false;}
    std::string    GetLabel          () const;      //!< human-readable label (streams \c Write to a string)
};

} //namespace
//! Polymorphic handle for any concrete \c Symmetry (shared, const) -- the spatial part of an \c Irrep.
export using sym_t=std::shared_ptr<const qchem::Symmetry::Symmetry>;
