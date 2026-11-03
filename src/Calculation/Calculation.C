// File: Calculation/Calculation.C  The abstract base of the three SCF facades (atom, molecule, solid).
//
// USER RULING (2026-11-02, V-FACADE-PARITY): the abstract base takes the name \c Calculation, the concrete molecular facade became
// \c MoleculeCalculation.  The base is where the TEMPLATE METHOD lives -- the high-level policy setup -> Prepare -> Converge ->
// Restart as a sequence of virtual calls, each step pure virtual or defaulted here and overridden only where the GEOMETRY differs --
// and where the shared state contract goes (V-FACADE-STATE: a virtual SaveState writing what every kind shares, each concrete
// facade adding its own).  It is born EMPTY: the steps move up one at a time, each flipping cells of the facade matrix
// (UTCalculation Facades/FacadeMatrix.Cell/<Kind>_<Capability>), so nothing lands here that two facades do not already share.
module;
export module qchem.Calculation;

export namespace qchem
{

//! \brief The abstract base of \c AtomCalculation, \c MoleculeCalculation and \c SolidCalculation (the template-method home).
//!
//! Data-free, derived from virtually (the project's interface convention).  A facade never opens a report run: it writes into the
//! run its CALLER opened (user ruling D, 2026-11-02).
class Calculation
{
public:
    virtual ~Calculation() = default;
};

} // namespace qchem
