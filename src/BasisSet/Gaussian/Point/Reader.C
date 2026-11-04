// File: BasisSet/Gaussian/Point/Reader.C  Abstract interface for a molecular basis-set reader.
module;

#include <vector>
#include <memory>
#include <algorithm>
export module qchem.BasisSet.Gaussian.Point.Reader;
export import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.GaussianRF;
export import qchem.Structure;
 
export namespace qchem::BasisSet::Gaussian
{
using namespace ::qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD;  // Cartesian glue moved out to PG_Cart_MnD

//-------------------------------------------------------------------------
//
//  Derived classes all read particular formats for the basis sets.
//
class Reader
{
public:
    virtual ~Reader() {};

    virtual GaussianRF*  ReadNext(const Site&)      =0;
    //! Position the reader at \a atom's entry.  THROWS when the file has no entry for that element: every
    //! caller goes on to ReadNext, so an absent element used to read the NEXT bytes as a shell (a garbage
    //! "angular-momentum label", 2026-10-21) -- an input error no basis construction can act on.
    virtual void             FindAtom(const Site&)      =0;
    virtual std::vector<int> GetLs   (           ) const=0;
};

struct ReadShell { std::unique_ptr<GaussianRF> rf; std::vector<int> Ls; };
//! Every shell \a bsr holds for every atom of \a cl -- the ONE copy of the loop the three Gaussian IrrepBasisSet readers
//! (PG_Cart, PG_Spherical, PG_LibCint) each used to carry verbatim.  When two shells on one atom share a radial (same centre,
//! same exponents within 0.1 % -- \c GaussianRF::operator==, which ignores L), a duplicate is the same radial AND the same
//! l-set; a same-exponent shell of ANOTHER l is its own block.  For EVERY basis, orbital and fit alike:
//!   - orbital bases since 2026-10-26 (OpenWork "Two shells of DIFFERENT l sharing one exponent ... MERGED -- silently": NiO's
//!     Ni p 0.30 beside the Ni d window starting at 0.30 loaded +2 functions, not +12);
//!   - the FIT bases since 2026-10-31 (OpenWork "The molecular FIT bases drop their D shells"): they kept a historical merge
//!     that never added a new l, and every A1 file puts SP and D shells on the same exponents, so N2's A1_coul loaded 32
//!     functions, not 68 -- the Coulomb fit 4-62 mHa low against PySCF's own A1_coul density fit.  The "Si2 size-consistency"
//!     that once argued for keeping it was an ill-posed SCF the d-less (spherical) Hartree had masked: a degenerate 3p^2 at
//!     integer occupation does not converge in PySCF either, with exact integrals or with the full fit.
inline std::vector<ReadShell> ReadShells(Reader& bsr, const Structure& cl)
{
    std::vector<ReadShell> out;
    for (const Site& site : cl.Sites())
    {
        bsr.FindAtom(site);
        while (std::unique_ptr<GaussianRF> rf{bsr.ReadNext(site)})
        {
            const std::vector<int> newLs=bsr.GetLs();
            bool duplicate=false;
            for (auto& s : out)
                if (*s.rf==*rf && s.Ls==newLs) { duplicate=true; break; }   // the same shell twice: drop the repeat
            if (!duplicate) out.push_back({std::move(rf), newLs});          // another l on a shared radial: its own block
        }
    }
    return out;
}

} //namespace qchem::BasisSet::Gaussian

