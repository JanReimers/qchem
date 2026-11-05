// File: BasisSet/Conditioning.C  The basis conditioning report: one explicit phase, the same on every facade that has an
// overlap (atom, molecule, Gaussian solid; plane waves are orthonormal).  doc/Pins.md pin 7.
module;
#include <cstddef>
export module qchem.BasisSet.Conditioning;
import qchem.BasisSet;

export namespace qchem::BasisSet
{

//! \brief Report how conditioned \a bs is, block by block, and judge the whole: \c perIrrep and \c removed rows and a
//! \c conditioning summary one level into the caller's context (its \c basis section), plus ONE \c [basis trim] line, printed
//! whether or not a run report is open.  Reads the analytic S only, so a facade can run it before any grid exists.  Edits
//! nothing.  S does not depend on spin, so there is one row per block.
//! \return the redundant-function count from the spectral-gap detector (0 == well conditioned).
template <class T> size_t ReportConditioning(const tBasisSet<T>& bs);

//! \brief The smallest eigenvalue of S over every block of \a bs: how strongly the orthonormal frame amplifies an error in a
//! matrix element (\f$\delta H' \lesssim \delta H/\lambda_{\min}\f$).  The solid facade's restart rule for a pilot run whose
//! first iterate was refused reads it (\c SolidCalculation::Converge: the start raster's alias must not exceed it).
template <class T> double MinOverlapEigenvalue(const tBasisSet<T>& bs);

} // namespace qchem::BasisSet
