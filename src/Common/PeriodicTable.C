// File: PeriodicTable.C  Implement a periodic table.
module;
#include <string>
#include <cassert>
#include <vector> 
// #include <fstream>
#include <nlohmann/json.hpp>


export module Common.PeriodicTable;
import qchem.Types;

export 
{
    const int N_Elements=110;

    class PeriodicTable
    {
    public:
        const char* GetSymbol(int              Z) const;
        int         GetZ     (const char* symbol) const;

        double GetEnergyHF            (int Z) const;
        double GetEnergyDFT           (int Z) const;
        double GetSlaterAlpha         (int Z) const;
        double GetNumUnpairedElectrons(int Z) const;
        int    GetMaxL                (int Z) const;
        int*   GetValanceConfiguration(int Z) const;
    private:
        static char theSymbols[N_Elements][3];
        static double EnergyHF   [N_Elements];
        static double EnergyDFT  [N_Elements];
        static double SlaterAlpha[N_Elements];
        static double NumUnpaired[N_Elements];
        static int    MaxL       [N_Elements];
        static int    ValConfig  [N_Elements][4];

    };

    struct ElementRecord
    {
        size_t      Z;
        std::string Symbol;
        size_t      NUnpaired;
        size_t      MaxL;
        size_t      ValConfig[4]; //spdf
        double      EnergyHF;     //Saito, Shiro L. Hartree–Fock–Roothaan energies and expectation values for the neutral atoms He to Uuo: The B-spline expansion method, Atomic Data and Nuclear Data Tables, 95,6, 836--870
        double      EnergyDFT;    //NIST https://math.nist.gov/DFTdata/atomdata/tables/ptable.html
    };


    struct OrbitalRecordSaito
    {
        OrbitalRecordSaito(nlohmann::json& j);
        std::string Symbol;
        double  Energy_HF;
        double  Energy_DFT;
        std::vector<double> r_moments; //<r^2>, <r^1>, <r^-1>, <r^-2>, <r^-3>, 
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
        std::vector<OrbitalRecordSaito> Orbitals;
    };

    std::ostream& operator<<(std::ostream& os, const OrbitalRecordSaito& o);
    std::ostream& operator<<(std::ostream& os, const ElementRecordSaito& e);

    class PeriodicTableSaito
    {
        public:
        PeriodicTableSaito();
        std::string GetSymbol(size_t Z) const {return get(Z).Symbol;}
        double GetEnergyHF            (size_t Z) const {return get(Z).Energy_HF;}
        double GetEnergyDFT           (size_t Z) const {return get(Z).Energy_DFT;}
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

