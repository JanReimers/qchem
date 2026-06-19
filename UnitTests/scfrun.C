// File: UnitTests/scfrun.C  Command-line SCF driver for tuning accelerators without recompiling.
//
// Example:
//   scfrun --Z 5 --model HF  --pol P --basis Slater     --acc Low    --accel Ladder --floor 1e-4 --stall 5
//   scfrun --Z 10 --model DHF --pol U --basis Slater_RKB --acc Medium --accel GDM    --emax 1e-3
//
#include <iostream>
#include <iomanip>
#include <string>
#include <map>
#include <nlohmann/json.hpp>

import qchem.BasisSet.Atom.Factory;
import qchem.Unittests.QchemTester;       // QchemTester, TestAtom, TestDiracAtom, BasisSetAccuracy
import qchem.Hamiltonian.Factory;         // Model, Pol, Factory
import qchem.BasisSet.Internal.DB_Cache_RAM; // theGlobalCache (the integrals cache)

using std::cout;
using std::endl;
using std::string;
using namespace qchem::Hamiltonian;       // Model, Pol, Factory
using namespace BasisSet::Atom;

// Config-driven concrete fixtures: the Hamiltonian model/polarization come from the CLI.
class CliAtom : public TestAtom
{
    Model m; Pol p;
public:
    CliAtom(int Z,int q,Model _m,Pol _p) : TestAtom(Z,q), m(_m), p(_p) {}
    virtual Hamiltonian* GetHamiltonian(cl_t& c) const override { return Factory(m,p,c); }
};
class CliDiracAtom : public TestDiracAtom
{
    Model m; Pol p;
public:
    CliDiracAtom(int Z,int q,Model _m,Pol _p) : TestDiracAtom(Z,q), m(_m), p(_p) {}
    virtual Hamiltonian* GetHamiltonian(cl_t& c) const override { return Factory(m,p,c); }
};

int main(int argc, char** argv)
{
    // ---- defaults ----
    int    Z=2, q=0, maxiter=50;
    string model="HF", pol="U", basis="", acc="Low", accel="DIIS";
    nlohmann::json accj;   // accelerator config passed to QchemTester
    // SCF convergence criteria (-1 => Z-scaled default); raise precision with --minfd/--virial.
    double minro=-1, minde=1e-5, virial=5e-1, minfd=-1, relax=0.5;

    // ---- help text ----
    auto usage=[&](std::ostream& os){
        os <<
        "Usage: scfrun [options]\n"
        "  Runs one SCF calculation for a single atom.  With no options it runs a quick\n"
        "  default (He, HF, DIIS).  Each option takes a value, written either as\n"
        "  \"--flag value\" or \"--flag=value\".\n"
        "\n"
        " System / method:\n"
        "  --Z <int>          atomic number                          (default 2)\n"
        "  --q <int>          net charge (electrons = Z - q)          (default 0)\n"
        "  --model <name>     HF | DHF | E1 | DE1                     (default HF)\n"
        "  --pol <U|P>        Unpolarized or Polarized                (default U)\n"
        "  --basis <name>     Slater|Gaussian|BSpline6|BSpliner6|Slater_RKB|Gaussian_RKB\n"
        "                       (default Slater, or Slater_RKB for Dirac models)\n"
        "  --acc <name>       N3|N5|Low|Medium|High pool accuracy     (default Low)\n"
        "  --accel <name>     DIIS | GDM | Ladder | directmin         (default DIIS)\n"
        "  --maxiter <int>    max SCF iterations                      (default 50)\n"
        "\n"
        " Accelerator tuning (forwarded as JSON to the accelerator factory):\n"
        "  --nproj <int>      DIIS subspace size\n"
        "  --emax <float>     GDM: start GDM steps once [F,D] < emax\n"
        "  --trust <float>    GDM: trust radius (max geodesic angle per step)\n"
        "  --ethresh <float>  Ladder: hand off only while |dE/E| > ethresh\n"
        "  --stall <int>      Ladder: # no-improve steps before a stall hand-off\n"
        "  --floor <float>    Ladder: [F,D] noise floor (never hand off below it)\n"
        "  --switchat <float> Ladder: hand off to the direct-min polisher once [F,D] < switchat\n"
        "\n"
        " Convergence criteria (-1 => Z-scaled default):\n"
        "  --minfd <float>    converge when [F,D] < minfd             (default Z*2e-5)\n"
        "  --minde <float>    converge when |dE/E| < minde            (default 1e-5)\n"
        "  --virial <float>   converge when |2+V/K| < virial          (default 0.5)\n"
        "  --minro <float>    converge when charge-density change < minro (default Z*1e-4)\n"
        "  --relax <float>    initial density relaxation factor       (default 0.5)\n"
        "\n"
        "  -h, --help         show this help and exit\n";
    };

    // ---- parse "--flag value" / "--flag=value" pairs ----
    for (int i=1;i<argc;i++)
    {
        string a=argv[i], inlineVal; bool hasInline=false;
        auto eq=a.find('=');                       // accept --flag=value as well as --flag value
        if (eq!=string::npos){ inlineVal=a.substr(eq+1); a.erase(eq); hasInline=true; }
        auto need=[&](int& i){
            if (hasInline) return inlineVal;
            if (i+1>=argc){cout<<"missing value for "<<argv[i]<<endl; exit(1);}
            return string(argv[++i]);
        };
        if      (a=="-h"||a=="--help") { usage(cout); return 0; }
        else if (a=="--Z")       Z=std::stoi(need(i));
        else if (a=="--q")       q=std::stoi(need(i));
        else if (a=="--model")   model=need(i);
        else if (a=="--pol")     pol=need(i);
        else if (a=="--basis")   basis=need(i);
        else if (a=="--acc")     acc=need(i);
        else if (a=="--accel")   accel=need(i);
        else if (a=="--maxiter") maxiter=std::stoi(need(i));
        else if (a=="--floor")   accj["floor"]=std::stod(need(i));
        else if (a=="--ethresh") accj["ethresh"]=std::stod(need(i));
        else if (a=="--stall")   accj["stall"]=std::stoi(need(i));
        else if (a=="--switchat")accj["switchat"]=std::stod(need(i));
        else if (a=="--trust")   accj["Trust"]=std::stod(need(i));
        else if (a=="--emax")    accj["EMax"]=std::stod(need(i));
        else if (a=="--nproj")   accj["NProj"]=(size_t)std::stoi(need(i));
        else if (a=="--minfd")   minfd=std::stod(need(i));
        else if (a=="--minde")   minde=std::stod(need(i));
        else if (a=="--virial")  virial=std::stod(need(i));
        else if (a=="--minro")   minro=std::stod(need(i));
        else if (a=="--relax")   relax=std::stod(need(i));
        else { cout<<"unknown option "<<argv[i]<<endl; usage(cout); return 1; }
    }
    bool dirac = (model=="DHF" || model=="DE1");
    if (basis.empty()) basis = dirac ? "Slater_RKB" : "Slater";
    accj["type"]=accel;

    BasisSet::theGlobalCache = new BasisSet::IntegralsCache_RAM<double>(true); //integrals cache (cf. gtestmain)

    // ---- string -> enum maps ----
    std::map<string,Model> models={{"HF",Model::HF},{"DHF",Model::DHF},{"E1",Model::E1},{"DE1",Model::DE1}};
    Pol pp = (pol=="P"||pol=="Polarized") ? Pol::Polarized : Pol::UnPolarized;
    using BT=BasisSet::Atom::Type;
    std::map<string,BT> bases={{"Slater",BT::Slater},{"Gaussian",BT::Gaussian},{"BSpline6",BT::BSpline6},
                               {"BSpliner6",BT::BSpliner6},{"Slater_RKB",BT::Slater_RKB},{"Gaussian_RKB",BT::Gaussian_RKB}};
    std::map<string,BasisSetAccuracy> accs={{"N3",BasisSetAccuracy::N3},{"N5",BasisSetAccuracy::N5},
                               {"Low",BasisSetAccuracy::Low},{"Medium",BasisSetAccuracy::Medium},{"High",BasisSetAccuracy::High}};
    if (!models.count(model)||!bases.count(basis)||!accs.count(acc)){cout<<"bad model/basis/acc"<<endl;return 1;}

    cout << "scfrun: Z="<<Z<<" q="<<q<<" model="<<model<<" pol="<<pol
         << " basis="<<basis<<" acc="<<acc<<" accel="<<accel<<" : "<<accj.dump()<<endl;

    // ---- run ----
    QchemTester* t = dirac ? (QchemTester*) new CliDiracAtom(Z,q,models[model],pp)
                           : (QchemTester*) new CliAtom     (Z,q,models[model],pp);
    t->SetAcceleratorConfig(accj);
    t->Init(accs[acc], bases[basis], false);
    if (minro<0) minro=Z*1e-4;
    if (minfd<0) minfd=Z*2e-5;
    //       NMaxIter   MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo MergeTol verbose
    t->Iterate({(size_t)maxiter, minro, minde, virial, minfd, relax, 1e-7, true});

    // ---- report ----
    double E=t->TotalEnergy();
    cout << "RESULT  E="<<std::setprecision(10)<<E
         << "  iters="<<t->GetIterationCount()
         << "  converged="<<(t->Converged()?"yes":"no") << endl;
    if (model=="HF")  t->RelativeHFError();
    if (model=="DHF") t->RelativeDHFError();
    delete t;
    return 0;
}
