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
} // namespace qchem
