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

    virtual GaussianRF*  ReadNext(const Atom&)      =0;
    //! Position the reader at \a atom's entry.  THROWS when the file has no entry for that element: every
    //! caller goes on to ReadNext, so an absent element used to read the NEXT bytes as a shell (a garbage
    //! "angular-momentum label", 2026-10-21) -- an input error no basis construction can act on.
    virtual void             FindAtom(const Atom&)      =0;
    virtual std::vector<int> GetLs   (           ) const=0;
};

//! \brief How a basis file's shells become (radial, l-set) blocks when two shells on one atom share a radial (same centre, same
//! exponents within 0.1 % -- \c GaussianRF::operator==, which ignores L).
//!   \c KeepEachL   -- ORBITAL bases: a duplicate is the same radial AND the same l-set; a same-exponent shell of ANOTHER l is
//!                    its own block (OpenWork §3 "Two shells of DIFFERENT l sharing one exponent ... MERGED -- silently",
//!                    fixed 2026-10-26: NiO's Ni p 0.30 beside the Ni d window starting at 0.30 loaded +2 functions, not +12).
//!   \c LegacyMerge -- the DENSITY-FIT bases: the historical merge, bit-identical, quirk and all (a genuinely new l is NOT
//!                    added -- the condition is inverted -- and the higher-Lmax radial replaces the old).  Flipping it moved the
//!                    fit bases' Si2 size-consistency 0.0026 -> 0.024 Ha: that path is non-additive and not yet understood.
enum class ShellMerge { KeepEachL, LegacyMerge };
struct ReadShell { std::unique_ptr<GaussianRF> rf; std::vector<int> Ls; };
//! Every shell \a bsr holds for every atom of \a cl, de-duplicated by \a m -- the ONE copy of the loop the three Gaussian
//! IrrepBasisSet readers (PG_Cart, PG_Spherical, PG_LibCint) each used to carry verbatim.
inline std::vector<ReadShell> ReadShells(Reader& bsr, const Structure& cl, ShellMerge m)
{
    std::vector<ReadShell> out;
    for (auto atom : cl)
    {
        bsr.FindAtom(*atom);
        while (std::unique_ptr<GaussianRF> rf{bsr.ReadNext(*atom)})
        {
            const std::vector<int> newLs=bsr.GetLs();
            bool duplicate=false;
            for (auto& s : out)
            {
                if (!(*s.rf==*rf)) continue;
                if (m==ShellMerge::KeepEachL)
                {
                    if (s.Ls==newLs) { duplicate=true; break; }   // the same shell twice: drop the repeat
                    continue;                                        // another l on this radial: its own block
                }
                duplicate=true;                                      // LegacyMerge, verbatim (see above)
                const bool useNewRF = *std::max_element(newLs.begin(),newLs.end()) > *std::max_element(s.Ls.begin(),s.Ls.end());
                for (auto l:newLs)
                    if (std::find(s.Ls.begin(),s.Ls.end(),l)!=s.Ls.end()) s.Ls.push_back(l);
                if (useNewRF) s.rf=std::move(rf);
                break;
            }
            if (!duplicate) out.push_back({std::move(rf), newLs});
        }
    }
    return out;
}

} //namespace qchem::BasisSet::Gaussian

