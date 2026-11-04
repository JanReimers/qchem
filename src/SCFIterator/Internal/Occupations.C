// File: SCFIterator/Internal/Occupations.C  The occupied configuration the cfg '*' compares, and the comparison.
module;
#include <cmath>
#include <cstddef>
#include <utility>
#include <vector>
#include <algorithm>

export module qchem.SCFIterator.Internal.Occupations;

export namespace qchem::SCFIterator
{
//! \brief THE OCCUPIED CONFIGURATION, as the aufbau decided it: per irrep block (in the wave function's own block order), each
//! orbital's occupation in the block's own orbital order -- (Irrep::SequenceIndex, occupations in block order).
//!
//! The cfg '*' compares two of these (AnElectronMoved), so it fires only when an electron moves -- between blocks, or to a
//! different orbital of its block.  It was an energy-ordered label string, (1a₁)²(2a₁)²...: exactly degenerate levels of DIFFERENT
//! irreps (CH₄'s t2 = a₁ + b₁ + b₂ in C2v; O₂'s ↑ and ↓ levels crossing) reorder on noise, so the string changed on most iterations
//! while no occupation did (OpenWork cfg row, 2026-11-03).  The merged level table is the wrong source too: a MergeTol host is
//! whichever level came first.
using OccupiedConfiguration=std::vector<std::pair<std::size_t,std::vector<double>>>;

//! \brief AN ELECTRON MOVED when some orbital's occupation changed by half an electron or more.
//!
//! The threshold is on the CHANGE, not on each occupation.  Rounding each occupation (lround) and comparing the integers put the
//! flag's edge AT a half-integer occupation, so a share sitting at exactly ½ flipped (1,0) <-> (0,1) on 1e-13 noise and flagged
//! while nothing moved (lane A, 2026-11-04).  Half an electron is the scale the rounding used; a smeared run's continuous drift
//! below it is not a move.  Blocks are matched by irrep; a block or orbital present on one side only counts as occupation 0.
inline bool AnElectronMoved(const OccupiedConfiguration& a, const OccupiedConfiguration& b)
{
    static const std::vector<double> none;
    auto blockOf=[](const OccupiedConfiguration& c, std::size_t irrep) -> const std::vector<double>&
    {
        for (const auto& [i,occ] : c) if (i==irrep) return occ;
        return none;
    };
    auto moved=[](const std::vector<double>& x, const std::vector<double>& y)
    {
        for (std::size_t k=0; k<std::max(x.size(),y.size()); ++k)
            if (std::fabs((k<x.size() ? x[k] : 0.0)-(k<y.size() ? y[k] : 0.0))>=0.5) return true;
        return false;
    };
    for (const auto& [i,occ] : a) if (moved(occ, blockOf(b,i))) return true;
    for (const auto& [i,occ] : b) if (moved(occ, blockOf(a,i))) return true;
    return false;
}
} // namespace qchem::SCFIterator
