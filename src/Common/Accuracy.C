// File: Common/Accuracy.C  THE ONE ACCURACY TARGET eps (doc/Pins.md pin 31), as a facade-neutral value (V-FACADE-ACCURACY, 2026-11-03).
//
// Pin 31 (extended 2026-11-03): one target eps governs every numerical approximation, for every structure type -- atom, molecule and
// solid alike.  This is the part every facade shares: the number and what it means for the run's client.  A term that consumes it
// derives its own tolerance from it (the GPW rasters: BasisSet::Gaussian::GridAccuracy, which adds the GPW calibration).  It lives
// this LOW (qcCommon) because the term libraries consume it and cannot import the facades.
module;
#include <cstddef>
export module qchem.Accuracy;

export namespace qchem
{

//! \brief WHAT THE RUN'S CLIENT CONSUMES of the Hamiltonian -- the measure the ONE accuracy target is judged on (user ruling
//! 2026-10-31, doc/GridEpsilonPlan.md §1): a single-point ENERGY is judged on the energy error (eps per atom), and anything that
//! produces forces, a geometry or a surface on the matrix elements (eps per element, max|dH_ij|: the energy error cancels inside
//! Tr D dH, M2c verdict 3).  The run ANSWERS it (who built it: \c Calculation::Consumes), the user never picks a metric.
enum class ClientConsumes { Energy, MatrixElements };

//! \brief The ONE accuracy target eps (deck `accuracy`): Ha per atom on an energy run, Ha per matrix element on a geometry run.
struct AccuracyTarget
{
    double eps = 1.0e-5;
    //! The target one term's delivered error is held to: eps x the atom count on an energy run, eps per element on a geometry run.
    double Target(ClientConsumes c, size_t nAtoms) const {return c==ClientConsumes::Energy ? eps*double(nAtoms) : eps;}
    bool operator==(const AccuracyTarget&) const = default;
};

} // namespace qchem
