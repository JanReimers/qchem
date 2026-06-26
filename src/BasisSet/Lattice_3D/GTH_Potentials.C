// File: BasisSet/Lattice_3D/GTH_Potentials.C  Reader for the CP2K GTH/HGH pseudopotential database.
//
// The CP2K GTH_POTENTIALS database is transcoded once, offline, into hierarchical JSON
// (doc/scripts/ParseGTH.py -> Data/gth_potentials.json); this module reads that JSON and builds the
// analytic local + separable-nonlocal potentials for a given (element, functional, valence) -- so the
// whole periodic table is available without hardcoded per-element factories.  See LocalPotential.C /
// SeparablePotential.C for the analytic HGH forms these populate.
module;
#include <string>
export module qchem.BasisSet.Lattice_3D.GTH_Potentials;
export import qchem.BasisSet.LocalPotential;
export import qchem.BasisSet.SeparablePotential;

export namespace BasisSet::Lattice_3D
{

//! \brief A complete GTH pseudopotential for one species: the ion (valence) charge plus the local and
//! separable-nonlocal parts.  \a zion is the natural charge for the ion-ion (Ewald) and electron-count.
struct GTH_PP
{
    int                    zion;
    HGH_LocalPotential     local;
    HGH_SeparablePotential nonlocal;
};

//! \brief Look up a GTH/HGH pseudopotential from the database (Data/gth_potentials.json) by element
//! symbol, exchange-correlation \a functional ("LDA","PBE","BLYP","BP","OLYP","HCTH120","HCTH407",
//! "PBESOL"), and valence \a q (0 => the database's default valence for that element/functional).
//! Throws std::runtime_error if the (element, functional, q) is not tabulated.  The functional MUST
//! match the SCF's XC functional (a pseudopotential is generated for a specific functional).
GTH_PP GetGTH(const std::string& element, const std::string& functional="LDA", int q=0);

} // namespace
