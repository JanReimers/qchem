// File: Common/Imp/PeriodicTable.C
module;
#include <string>
#include <cassert>
#include <vector> 
#include <filesystem>
#include <fstream>
#include <iostream>
#include <set>
#include <nlohmann/json.hpp>

module qchem.PeriodicTable;


std::ostream& operator<<(std::ostream& os, const OrbitalRecordSaito& o)
{
    os << std::fixed;
    os << "    " << o.Symbol << " " << std::setw(14) << std::setprecision(11) << o.Energy_HF;
    os << std::fixed;
    for (auto r:o.r_moments)
        os << " " << std::setw(11) << std::setprecision(6) << r;
    return os;
}

std::ostream& operator<<(std::ostream& os, const ElementRecordSaito& e)
{
    os << e.Symbol << "(" << e.Z << ") " << e.ValConfigString << " " << e.Term
    << " E_HF=" << std::setw(14) << std::setprecision(6) << e.Energy_HF
    << " E_DHF=" << std::setw(14) << std::setprecision(6) << e.Energy_DHF
    << " NUnpaired=" << e.NUnpaired << std::endl;
    // for (auto o:e.Orbitals)
    //     os << o << std::endl;
    return os;
}

OrbitalRecordSaito::OrbitalRecordSaito(nlohmann::json& j) : Symbol(j["name"]),Energy_HF(j["e"])
{
    for (auto r:j["rs"])
        r_moments.push_back(r);
}

std::set<size_t> fullShells({2,10,18,36,54,86,118});

// Relativistic orbital labels in physical (n,l,j) order, matching doc/DHF_GS_Energies_rel.json.
const std::vector<std::string> theDHFOrbitalLabels({
    "1s+",
    "2s+", "2p-", "2p+",
    "3s+", "3p-", "3p+", "3d-", "3d+",
    "4s+", "4p-", "4p+", "4d-", "4d+", "4f-", "4f+",
    "5s+", "5p-", "5p+", "5d-", "5d+", "5f-", "5f+",
    "6s+", "6p-", "6p+", "6d-", "6d+",
    "7s+",
});

#ifndef COMMON_DATA_PATH
#error "COMMON_DATA_PATH must be defined by CMake"
#endif

// Construct-on-first-use: this path is read from the ctors of static PeriodicTable* objects that may
// live in other translation units (e.g. QchemTester::itsPT).  A namespace-scope global would be a
// static-initialization-order fiasco -- it could still be zero-initialized when such a ctor runs,
// giving a null path string and a crash.  A function-local static is initialized on first call.
static const std::filesystem::path& common_data_dir()
{
    static const std::filesystem::path p = COMMON_DATA_PATH;
    return p;
}

ElementRecordSaito::ElementRecordSaito(nlohmann::json& j) : Z(j["Z"]), Symbol(j["symbol"]),  ValConfigString(j["valance"]), Term(j["term"]), MaxL(0), Energy_HF(j["HFEnergy"])
{
    for (auto o:j["Orbitals"])
        Orbitals.push_back(OrbitalRecordSaito(o));
    NUnpaired=std::atoi(&Term[0])-1;
    assert(NUnpaired<=8);
    size_t spos=ValConfigString.find('s'),ppos=ValConfigString.find('p'),dpos=ValConfigString.find('d'),fpos=ValConfigString.find('f');
    bool bs=spos!=std::string::npos,bp=ppos!=std::string::npos,bd=dpos!=std::string::npos,bf=fpos!=std::string::npos;
    std::string nss="0",nps="0",nds="0",nfs="0";

    if (bs)
        nss=ValConfigString.substr(spos+1,1); //number of s electron is always one digit.
    if (bp)
        nps=ValConfigString.substr(ppos+1,1); //number of p electron is always one digit.
    if (bd)
    {
        if (bp)
            nds=ValConfigString.substr(dpos+1,ppos-dpos-2); //d is always just before p.
        else
            nds=ValConfigString.substr(dpos+1,ValConfigString.length()-dpos-1);
    }
    if (bf)
    {
        if (bd)
            nfs=ValConfigString.substr(fpos+1,dpos-fpos-2); //f is usually just before d.
        else if (bp)
            nfs=ValConfigString.substr(fpos+1,ppos-fpos-2); //No d, so f is just before p.
        else
            nfs=ValConfigString.substr(fpos+1,ValConfigString.length()-fpos-1);
    }
    ValConfig[0]=stoi(nss);
    ValConfig[1]=stoi(nps);
    ValConfig[2]=stoi(nds);
    ValConfig[3]=stoi(nfs);

    if (Z>4) MaxL=1;
    if (Z>20) MaxL=2;
    if (Z>57) MaxL=3;

    for (size_t l=0; l<4; l++)   // (plain loop: PeriodicTable no longer pulls qchem.Types for iv_t)
    {
        assert(ValConfig[l]<=2*(2*l+1));
        if (fullShells.find(Z)!=fullShells.end()) ValConfig[l]=0; //Clear out full shells.
    }
}


PeriodicTableSaito::PeriodicTableSaito()
{
    // Read in Saito HF data.
    {
        std::ifstream file(common_data_dir() / "saito.json");
        assert(file);
        nlohmann::json jsondata;
        file >> jsondata;
        for (auto e:jsondata)
            elements.push_back(ElementRecordSaito(e));
    }
    // Read in NIST DFT data
    {
        std::ifstream file(common_data_dir() / "nistLDA.json");
        assert(file);
        nlohmann::json jsondata;
        file >> jsondata;
        for (auto& jse:jsondata)
        {
            size_t Z=jse["Z"].template get<size_t>();
            ElementRecordSaito& e=elements[Z-1];
            e.Energy_DFT=jse["Etot"].template get<double>();
            for (auto& o:e.Orbitals)
            {
                o.Energy_DFT=jse[o.Symbol].template get<double>();
            }
        }
    }
    // Read in relativistic DHF data (total energy + spin-orbit split orbital eigenvalues).
    {
        std::ifstream file(common_data_dir() / "DHF_GS_Energies_rel.json");
        assert(file);
        nlohmann::json jsondata;
        file >> jsondata;
        for (auto& jse:jsondata)
        {
            size_t Z=jse["Z"].template get<size_t>();
            ElementRecordSaito& e=elements[Z-1];
            e.Energy_DHF=jse["TE_rel"].template get<double>();
            for (const std::string& label:theDHFOrbitalLabels)
            {
                const auto& val=jse[label];
                if (val.is_null()) continue;
                e.DHFOrbitals.push_back({label,val.template get<double>()});
            }
        }
    }

}

const PeriodicTableSaito& thePeriodicTable()
{
    static const PeriodicTableSaito theTable;   // loaded once, on first use
    return theTable;
}

size_t PeriodicTableSaito::GetZ(const std::string& symbol) const
{
    for (size_t i=0; i<elements.size(); i++)
        if (elements[i].Symbol==symbol) return i+1;   // Z = index+1
    return 0;                                         // not found
}

// Schwarz X-alpha optimized exchange parameters (J.C. Slater / K. Schwarz).  Tabulated only for
// H..Ca; everything else defaults to 0.70 via GetSlaterAlpha.  Z-indexed (slot 0 is a filler).
namespace
{
    const double theSlaterAlpha[] =
    {
        0.0,     //filler
        0.77679, //H
        0.77224, //He
        0.79118, //Li
        0.79526, //Be
        0.78744, //B
        0.77657, //C
        0.76654, //N
        0.76454, //O
        0.75954, //F
        0.73100, //Ne
        0.75110, //Na
        0.74942, //Mg
        0.74797, //Al
        0.74521, //Si
        0.74309, //P
        0.74270, //S
        0.74183, //Cl
        0.722,   //Ar
        0.7227,  //K
        0.7220,  //Ca
    };
    const size_t nSlaterAlpha = sizeof(theSlaterAlpha)/sizeof(theSlaterAlpha[0]);
}

double PeriodicTableSaito::GetSlaterAlpha(size_t Z) const
{
    assert(Z>0);
    double ret = (Z<nSlaterAlpha) ? theSlaterAlpha[Z] : 0.0;
    if (ret==0.0) ret=0.70;   // default for un-tabulated / zero entries
    return ret;
}

// Pauling electronegativity, Z-indexed (slot 0 filler).  Standard textbook scale, H..Rn.  Noble gases (and
// a few unmeasured species) are left 0.0 -- they do not form ions, so the IonicSAD heuristic skips a 0.
namespace
{
    const double thePaulingEN[] =
    {
        0.0,                                                     // filler
        2.20, 0.0,                                               // H  He
        0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, 0.0,           // Li Be B  C  N  O  F  Ne
        0.93, 1.31, 1.61, 1.90, 2.19, 2.58, 3.16, 0.0,           // Na Mg Al Si P  S  Cl Ar
        0.82, 1.00,                                              // K  Ca
        1.36, 1.54, 1.63, 1.66, 1.55, 1.83, 1.88, 1.91, 1.90, 1.65,  // Sc..Zn
        1.81, 2.01, 2.18, 2.55, 2.96, 3.00,                      // Ga Ge As Se Br Kr
        0.82, 0.95,                                              // Rb Sr
        1.22, 1.33, 1.60, 2.16, 1.90, 2.20, 2.28, 2.20, 1.93, 1.69,  // Y..Cd
        1.78, 1.96, 2.05, 2.10, 2.66, 2.60,                      // In Sn Sb Te I  Xe
        0.79, 0.89,                                              // Cs Ba
        1.10, 1.12, 1.13, 1.14, 1.13, 1.17, 1.20, 1.20, 1.10, 1.22, 1.23, 1.24, 1.25, 1.10, 1.27,  // La..Lu
        1.30, 1.50, 2.36, 1.90, 2.20, 2.20, 2.28, 2.54, 2.00,    // Hf Ta W  Re Os Ir Pt Au Hg
        1.62, 2.33, 2.02, 2.00, 2.20, 0.0,                       // Tl Pb Bi Po At Rn
    };
    const size_t nPaulingEN = sizeof(thePaulingEN)/sizeof(thePaulingEN[0]);
}

double PeriodicTableSaito::GetElectronegativity(size_t Z) const
{
    assert(Z>0);
    return (Z<nPaulingEN) ? thePaulingEN[Z] : 0.0;   // 0.0 = un-tabulated / noble gas (no ion)
}

