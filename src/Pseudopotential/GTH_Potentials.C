// File: BasisSet/Lattice/GTH_Potentials.C  Reader for the CP2K GTH/HGH pseudopotential database.
//
// The CP2K GTH_POTENTIALS database is transcoded once, offline, into hierarchical JSON
// (doc/scripts/ParseGTH.py -> Data/gth_potentials.json); this module reads that JSON and builds the
// analytic local + separable-nonlocal potentials for a given (element, functional, valence) -- so the
// whole periodic table is available without hardcoded per-element factories.  See LocalPotential.C /
// SeparablePotential.C for the analytic HGH forms these populate.
module;
#include <memory>
#include <string>
#include <utility>
#include <vector>
export module qchem.Pseudopotential.GTH_Potentials;
export import qchem.Pseudopotential.LocalPotential;
export import qchem.Pseudopotential.SeparablePotential;

export namespace qchem::Pseudopotential
{

//! \brief A complete GTH pseudopotential for one species: the ion (valence) charge plus the local and
//! separable-nonlocal parts.  \a zion is the natural charge for the ion-ion (Ewald) and electron-count.
struct GTH_PP
{
    int                    zion;
    HGH_LocalPotential     local;
    HGH_SeparablePotential nonlocal;
    //! The XC functional this potential was FITTED with, as CP2K's primary name spells it: "PADE" (Goedecker-Teter-Hutter Padé LDA) for
    //! every entry looked up as "LDA", else the lookup key ("PBE", "BLYP", ...).  A run comparing its own functional against the PP's reads
    //! this (OpenWork §3 "A GTH-PADE pseudopotential runs under VWN5"); the lookup key names only the family.
    std::string            fitted;
};

//! \brief Look up a GTH/HGH pseudopotential from the database (Data/gth_potentials.json) by element
//! symbol, exchange-correlation \a functional ("LDA","PBE","BLYP","BP","OLYP","HCTH120","HCTH407",
//! "PBESOL"), and valence \a q (0 => the database's default valence for that element/functional).
//! Throws std::runtime_error if the (element, functional, q) is not tabulated.  The functional MUST
//! match the SCF's XC functional (a pseudopotential is generated for a specific functional).
GTH_PP GetGTH(const std::string& element, const std::string& functional="LDA", int q=0);

//! The valence variants (q values) tabulated for \a element / \a functional, ASCENDING -- the smallest is the
//! valence-only PP (largest core), larger ones are semicore.  Empty if the (element,functional) is absent.
//! Lets a driver choose valence-vs-semicore by policy (the DATABASE default is often the semicore, not the
//! valence, so a valence-first tool must select explicitly).  Throws nothing.
std::vector<int> GetGTHValences(const std::string& element, const std::string& functional="LDA");

//! \brief The multi-species LOCAL pseudopotential of a cell: one GTH local per (element, valence) in \a species, keyed by atomic number.
//! PUBLIC (it was Hamiltonians.C's file-local builder until 2026-11-03) because TWO clients must hold the SAME local: the periodic
//! Hamiltonian's PP terms, and the GPW facade, which reads the local's sharpness off it to size the raster ladder (V-BS-RASTER:
//! V_loc long's audit row is a pure function of the ladder and that sharpness).  One builder, so the two agree by construction.
std::shared_ptr<const LocalPotential> MakeMultiSpeciesLocal(const std::vector<std::pair<std::string,int>>& species,
                                                            const std::string& functional);

} // namespace
