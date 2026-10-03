//! \file WaveFunction/EigenRows.C
//! \brief The ROWS of the eigenvalue/occupation level tables, as pure values (D10).
//!
//! WHY THIS EXISTS.  `DisplayEigen` used to build each row and render it in ONE pass, so nothing about a row
//! was assertable -- and "printing those tables has been a constant source of bugs" (user, 2026-08-10).  The
//! bug list is the evidence: doubly-empty levels dropped BELOW the frontier, `setprecision(0)` rounding a
//! smeared 0.996 to "1/1", and an ABSENT channel's energy filled from the other channel (a row that read as a
//! hole, with the exchange splitting printing exactly 0.00000000).  Every one is a property of the ROWS; every
//! one was found by reading a run log.  These builders take the \c EnergyLevels and return the cell strings a
//! reader sees, so the next one is a unit test (src/WaveFunction/tests/EigenRowsUT.C).  The renderer
//! (tCompositeWF::DisplayEigen*) only draws them.
//!
//! KNOWN LIMIT, kept visible rather than hidden: polarized rows are still PAIRED by (n, symmetry), and \c n indexes
//! a DEGENERATE GROUP whose grouping can differ between spin channels once the exchange splitting is nonzero --
//! pairing by energy/character is the open remainder of D10 (doc/CleanCode.md).
module;
#include <string>
#include <vector>
export module qchem.WaveFunction.EigenRows;
export import qchem.EnergyLevel;     // EnergyLevels, EnergyLevel, Orbital_QNs

export namespace qchem::WaveFunction
{

//! \brief One row of the UNPOLARIZED table: occupation/degeneracy, energy, "<n><symmetry>"; \c l colours the row.
struct EigenRow
{
    std::string occ, e, label;
    size_t      l=0;
};

//! \brief One row of the POLARIZED table.  A channel that does not carry the level prints "--" in every one of its
//! cells and contributes no splitting -- ABSENT IS NOT EMPTY.  \c dnDim is true when the down channel is empty or
//! absent (the renderer greys it).
struct PolarizedEigenRow
{
    std::string occUp, eUp, label, occDn, eDn, dE;
    size_t      l=0;
    bool        dnDim=false;
};

//! \brief The unpolarized rows of \a els: every level up to the HIGHEST OCCUPIED one (never stopping short of it -- a
//! MOM-pinned run can leave a level empty BELOW an occupied one, and that hole is the row a reader needs), then
//! the first empty level past it ends the table.  Stops by OCCUPATION, not by the sign of the energy.
std::vector<EigenRow> UnpolarizedEigenRows(const Orbitals::EnergyLevels& els);

//! \brief The polarized rows.  \a combined supplies the order and the frontier; \a up / \a dn the per-channel levels,
//! looked up by (n, symmetry).  A doubly-empty level is dropped only when it sits ABOVE the highest occupied
//! energy over both channels.
std::vector<PolarizedEigenRow> PolarizedEigenRows(const Orbitals::EnergyLevels& combined,
                                                  const Orbitals::EnergyLevels& up,
                                                  const Orbitals::EnergyLevels& dn);

//! \brief "n/deg" with the occupation as an integer when it is one (gapped insulator) and with two decimals otherwise
//! (Fermi-smeared) -- a smeared 0.996 must NOT print as "1/1".
std::string OccupationCell(double occ, int degen);

} // namespace qchem::WaveFunction
