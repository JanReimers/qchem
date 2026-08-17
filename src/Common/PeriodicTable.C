// File: PeriodicTable.C  Implement a periodic table.
module;
#include <string>
#include <cassert>
#include <cmath>      // std::fabs (the oracle helpers below); a std header, NOT qchem.Math -- qcCommon stays math-free
#include <iostream>   // the oracle helpers' ppm/ppb/ppt reporter (module TUs must include <ostream>/<iostream> THEMSELVES to stream literals -- the R1.9 lesson)
#include <vector>
// #include <fstream>
#include <nlohmann/json.hpp>


export module qchem.PeriodicTable;
// (no math import: the periodic table is pure data -- string/double/json -- so qcCommon stays math-free.)

export
{
    const int N_Elements=110;

        struct OrbitalRecordSaito
    {
        OrbitalRecordSaito(nlohmann::json& j);
        std::string Symbol;
        double  Energy_HF;
        double  Energy_DFT;
        std::vector<double> r_moments; //<r^2>, <r^1>, <r^-1>, <r^-2>, <r^-3>,
    };

    // Relativistic Dirac-Hartree-Fock orbital eigenvalue, labelled by the spin-orbit
    // split symbol, e.g. "1s+", "2p-", "2p+".  Source: doc/DHF_GS_Energies_rel.json.
    struct DHFOrbitalRecord
    {
        std::string Label;
        double      Energy;
    };

    struct ElementRecordSaito
    {
        ElementRecordSaito(nlohmann::json& j);

        size_t      Z;
        std::string Symbol;
        std::string ValConfigString;
        std::string Term;
        size_t      NUnpaired;
        size_t      MaxL;
        size_t      ValConfig[4]; //spdf
        double      Energy_HF;  //Saito, Shiro L. Hartree–Fock–Roothaan energies and expectation values for the neutral atoms He to Uuo: The B-spline expansion method, Atomic Data and Nuclear Data Tables, 95,6, 836--870
        double      Energy_DFT; //NIST https://math.nist.gov/DFTdata/atomdata/tables/ptable.html
        double      Energy_DHF=0.0; //Relativistic DHF total energy (TE_rel) from doc/DHF_GS_Energies_rel.json
        std::vector<OrbitalRecordSaito> Orbitals;
        std::vector<DHFOrbitalRecord>   DHFOrbitals;
    };

    std::ostream& operator<<(std::ostream& os, const OrbitalRecordSaito& o);
    std::ostream& operator<<(std::ostream& os, const ElementRecordSaito& e);

    class PeriodicTableSaito
    {
        public:
        PeriodicTableSaito();
        std::string GetSymbol(size_t Z) const {return get(Z).Symbol;}
        size_t      GetZ     (const std::string& symbol) const;        //!< reverse symbol -> Z lookup (0 if not found)
        size_t      GetNumElements() const {return elements.size();}   //!< tabulated element count (max valid Z)
        double GetSlaterAlpha         (size_t Z) const;                //!< Schwarz X-alpha optimized exchange parameter (0.70 default)
        double GetElectronegativity   (size_t Z) const;                //!< Pauling electronegativity (0.0 = noble gas / un-tabulated)
        double GetEnergyHF            (size_t Z) const {return get(Z).Energy_HF;} //Saito, Shiro L. Hartree–Fock–Roothaan energies and expectation values for the neutral atoms He to Uuo: The B-spline expansion method, Atomic Data and Nuclear Data Tables, 95,6, 836--870
        double GetEnergyDFT           (size_t Z) const {return get(Z).Energy_DFT;} //NIST https://math.nist.gov/DFTdata/atomdata/tables/ptable.html
        double GetEnergyDHF           (size_t Z) const {return get(Z).Energy_DHF;} //Tatewaki et.al. ACS Omega 2017, 2, 9, 6072–6080
        const std::vector<DHFOrbitalRecord>& GetDHFOrbitals(size_t Z) const {return get(Z).DHFOrbitals;}
        double GetNumUnpairedElectrons(size_t Z) const {return get(Z).NUnpaired;}
        int    GetMaxL                (size_t Z) const {return get(Z).MaxL;}
        const size_t*   GetValanceConfiguration(size_t Z) const {return &(get(Z).ValConfig[0]);}
    private:
        const ElementRecordSaito& get(size_t Z) const
        {
            assert(Z>=1);
            assert(Z<=elements.size());
            return elements[Z-1];
        }
        std::vector<ElementRecordSaito> elements;
    };

    //! The one shared periodic table.  Construct-on-first-use (function-local static): the ctor reads several
    //! JSON files, so a single instance replaces the ~handful of per-TU copies that each reloaded everything;
    //! the function-local static also side-steps the static-initialization-order fiasco (see common_data_dir).
    //! Returned by const reference -- the table is immutable reference data, safe to share everywhere.
    const PeriodicTableSaito& thePeriodicTable();
} //export block

// Oracle helpers (moved from qchem.Unittests.TestUtils, R2.20): thin wrappers over the reference-energy
// tables above.  They are PRODUCTION data consumers -- CLIapps/scfrun reports against them -- so they live
// beside the tables they wrap, not in a test module.
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

//! Z-keyed oracle checks: the computed atomic energy E vs the stored NIST (HF/DFT) or Dirac (DHF)
//! reference for element Z.  All return the SIGNED relative error (callers bound fabs of it).
inline double RelativeHFError (double E, int Z, bool quiet = false)
    {return RelativeError(E, thePeriodicTable().GetEnergyHF (Z), quiet);}
inline double RelativeDFTError(double E, int Z, bool quiet = false)
    {return RelativeError(E, thePeriodicTable().GetEnergyDFT(Z), quiet);}
inline double RelativeDHFError(double E, int Z, bool quiet = false)
    {return RelativeError(E, thePeriodicTable().GetEnergyDHF(Z), quiet);}
} // namespace qchem

