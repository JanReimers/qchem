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
    //! \brief Brillouin-zone integration weight of this symmetry's block.  1 for ordinary symmetries; a
    //! Bloch k-point returns the UNIFORM per-point sampling weight \f$1/N_\mathrm{mesh}\f$ -- its k-star
    //! multiplicity lives in \c GetDegeneracy (the ATOM SHELL CONVENTION: one stored representative per
    //! star, like an atom l-shell), so the BZ average is \f$\sum_\mathrm{blocks} w\,\rho_\mathrm{block}\f$
    //! with the star already inside each block's occupations/density.
    virtual double GetWeight         () const {return 1.0;}
    //! \brief Can this irrep's basis functions be chosen REAL (making S, T, V real-symmetric)?
    //! The basis-type half of doc/RealComplexPlan.md's rule: block is real ⇔ irrep.IsReal() ∧ every
    //! term.PreservesReal().  Default \c true -- point-group irreps, atomic \c Yl / \c Ylm shells and
    //! \c UnitQN all admit real functions (complex point-group irreps are a CONVENTION, not a necessity:
    //! real combinations always exist for a real H, and the SALC machinery already takes them).  The one
    //! override is \c BlochQN: \f$\chi^k=\sum_R e^{ikR}\chi\f$ is real iff \f$2k\f$ is a reciprocal
    //! lattice vector (a TRIM point) -- an EXACT integer fact there, no tolerance.
    virtual bool   IsReal             () const {return true;}
    //! \brief Does this irrep label already carry spin (spin-orbit coupled)?  \c false for ordinary spatial
    //! symmetries -- spin is layered on top by \c Irrep, whose \a ms is the non-relativistic spin channel.
    //! \c true only for the relativistic double group (\c SphericalSpinor), where \f$\kappa\f$ encodes
    //! \f$j=l\pm\tfrac12\f$; there the layer-2 \a ms is \c Spin::None (spin is already in \f$\kappa\f$).
    virtual bool   CarriesSpin        () const {return false;}
    //! \brief May energy levels of DIFFERENT irreps of this symmetry kind form ONE degenerate shell when
    //! their eigenvalues coincide (within the caller's MergeTol)?  \c false for ordinary spatial symmetries:
    //! an accidental cross-irrep degeneracy is NOT a shell, and atom/molecule configuration labels must stay
    //! distinct.  \c true for crystal k-points (\c BlochQN): a k-STAR of band levels is one
    //! symmetry-degenerate shell of the crystal group, split across k-blocks only by bookkeeping -- the
    //! \c EnergyLevels reporting layer merges them so an equal-eigenvalue star displays as one level (and
    //! the SCF cfg change-flag stops flapping on the ULP ordering of exact cross-k ties).
    virtual bool   MergeAcrossIrreps  () const {return false;}
    std::string    GetLabel          () const;      //!< human-readable label (streams \c Write to a string)
};

} //namespace
//! Polymorphic handle for any concrete \c Symmetry (shared, const) -- the spatial part of an \c Irrep.
export using sym_t=std::shared_ptr<const qchem::Symmetry::Symmetry>;
