// File: UnitTests/TestUtils.C  Free test-utility helpers shared by facade-based molecular tests.
//
// RelativeError is the (Eref-E)/Eref check with the ppm/ppb/ppt reporter, lifted out of
// QchemTester::RelativeError so that tests driving qchem::Calculation (which has no tester) can use the
// same relative-energy assertion + console report.  Pure helper -- no calculation state.  The oracle
// helpers (RelativeHFError/DFTError/DHFError vs NIST atomic energies) stay Z-keyed on the atom scaffold.
module;
#include <cmath>
#include <iostream>
export module qchem.Unittests.TestUtils;
import qchem.PeriodicTable;   // thePeriodicTable(): NIST/Dirac atomic reference energies (the oracle)

export namespace qchem
{
//! Signed relative error (Eref - E)/Eref of a computed energy E against the reference Eref.  Prints a
//! human ppm/ppb/ppt report (suppress with quiet=true).  Caller asserts on fabs(RelativeError(...)).
inline double RelativeError(double E, double Eref, bool quiet = false)
{
    const double error = (Eref - E) / Eref;
    if (!quiet)
    {
        std::cout.precision(9);
        std::cout << "E relative error=" << error * 100.0 << "%, ";
        std::cout.precision(2);
        if (std::fabs(error) > 1e-7)       std::cout << error * 1e6  << "(ppm)" << std::endl;
        else if (std::fabs(error) > 1e-10) std::cout << error * 1e9  << "(ppb)" << std::endl;
        else                               std::cout << error * 1e12 << "(ppt)" << std::endl;
    }
    return error;
}

// Z-keyed oracle checks: the computed atomic energy E vs the stored NIST (HF/DFT) or Dirac (DHF) reference
// for element Z.  Lifted out of QchemTester::RelativeHF/DFT/DHFError so AtomCalculation-driven atom tests
// share them.  All return the SIGNED relative error (the atom tests bound it as the scaffold did).
inline double RelativeHFError (double E, int Z, bool quiet = false)
    {return RelativeError(E, thePeriodicTable().GetEnergyHF (Z), quiet);}
inline double RelativeDFTError(double E, int Z, bool quiet = false)
    {return RelativeError(E, thePeriodicTable().GetEnergyDFT(Z), quiet);}
inline double RelativeDHFError(double E, int Z, bool quiet = false)
    {return RelativeError(E, thePeriodicTable().GetEnergyDHF(Z), quiet);}
} // namespace qchem
