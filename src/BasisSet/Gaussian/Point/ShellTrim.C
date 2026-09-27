// File: BasisSet/Gaussian/Point/ShellTrim.C  The VET-STAGE basis trim, as a VALUE, and the reader decorator that
// applies it (doc/Pins.md pin 22; doc/OpenWork.md §2 "Diffuse-basis ACTUATOR").
//
// WHAT IT REPLACES.  A near-dependent diffuse shell used to be removed in one of two ways, both refuted: by
// hand-editing (or script-rewriting) the committed .bsd file -- the run could not say which span it used --
// or by the ORTHO-TIME pivoted Cholesky, which drops individual AOs PER k-BLOCK, greedily: on NiO AFM-II
// (2026-09-27) it dropped Ni1's alpha=0.06 s in some k-blocks, Ni2's in others, both or neither elsewhere --
// a basis that differs from k to k and treats symmetry-equivalent sites inequivalently.
//
// WHAT IT IS.  A trim names SHELLS the way a basis file does: (element, l, exponents).  It is applied WHILE
// the basis is READ, so it is decided once, before anything downstream is built (pin 22 ruling a), holds for
// every k-block alike (ruling b), and -- being per ELEMENT -- removes the shell from EVERY site of that
// species, a union of whole symmetry orbits under any group (ruling c).  This is exactly the hand practice
// it automates ("copy the atom's block, remove the diffuse function"), minus the edited file.
// ⚠ Per ELEMENT is coarser than per ORBIT: two crystallographically INEQUIVALENT sites of one species share
// the trim.  That is the conservative direction (both keep the same span), and it is what the .bsd format
// itself can express; a per-orbit trim is a later refinement if a material ever needs it.
module;
#include <iosfwd>
#include <vector>
export module qchem.BasisSet.Gaussian.Point.ShellTrim;
export import qchem.BasisSet.Gaussian.Point.Reader;
import qchem.Types;

export namespace qchem::BasisSet::Gaussian
{

//! One shell, named the way a basis file names it.
struct TrimmedShell
{
    int    Z=0;          //!< the element
    int    l=0;          //!< the angular momentum removed (a multi-l shell, e.g. SP, loses only this l)
    rvec_t exponents;    //!< the shell's primitive exponents -- its identity within the element's block
};

//! \brief The set of shells a basis is read WITHOUT.  Empty = the file as written.
struct ShellTrim
{
    std::vector<TrimmedShell> shells;
    bool empty() const {return shells.empty();}
    //! Does this trim remove angular momentum \a l of the shell with \a exponents on element \a Z?
    //! Exponents match to 1e-10 relative (they are read from the same file text).
    bool Removes(int Z, int l, const rvec_t& exponents) const;
    //! "Ni l=0 {0.06}; ..." -- the decision reported as a BASIS, never as AO indices (pin 22).
    std::ostream& Write(std::ostream&) const;
};

//! \brief A \c Reader DECORATOR that skips the trimmed shells of the reader it wraps.  The basis classes
//! consume it through the unchanged \c Reader face, so no basis construction learns about trimming.
class TrimmingReader : public Reader
{
public:
    //! \a inner and \a trim must outlive this reader (it is a construction-time adaptor).
    TrimmingReader(Reader& inner, const ShellTrim& trim) : itsInner(inner), itsTrim(trim) {}
    virtual GaussianRF*      ReadNext(const Atom&) override;
    virtual bool             FindAtom(const Atom& a) override {return itsInner.FindAtom(a);}
    virtual std::vector<int> GetLs   () const override {return itsLs;}
private:
    Reader&          itsInner;
    const ShellTrim& itsTrim;
    std::vector<int> itsLs;   //!< the Ls of the last shell returned, minus the trimmed ones
};

} // namespace
