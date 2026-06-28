// File: PeriodicTable.C  Implement a periodic table.
module;
#include <string>
#include <cassert>
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
        double GetSlaterAlpha         (size_t Z) const;                //!< Schwarz X-alpha optimized exchange parameter (0.70 default)
        double GetEnergyHF            (size_t Z) const {return get(Z).Energy_HF;}
        double GetEnergyDFT           (size_t Z) const {return get(Z).Energy_DFT;}
        double GetEnergyDHF           (size_t Z) const {return get(Z).Energy_DHF;}
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
} //export block

