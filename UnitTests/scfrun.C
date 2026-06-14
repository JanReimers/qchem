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

import qchem.Unittests.QchemTester;       // QchemTester, TestAtom, TestDiracAtom, BasisSetAccuracy
import qchem.Hamiltonian.Factory;         // Model, Pol, Factory
import qchem.BasisSet.Internal.DB_Cache_RAM; // theGlobalCache (the integrals cache)

using std::cout;
using std::endl;
using std::string;
using namespace qchem::Hamiltonian;       // Model, Pol, Factory

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

    // ---- parse "--flag value" pairs ----
    auto need=[&](int& i){ if (i+1>=argc){cout<<"missing value for "<<argv[i]<<endl; exit(1);} return string(argv[++i]); };
    for (int i=1;i<argc;i++)
    {
        string a=argv[i];
        if      (a=="--Z")       Z=std::stoi(need(i));
        else if (a=="--q")       q=std::stoi(need(i));
        else if (a=="--model")   model=need(i);
        else if (a=="--pol")     pol=need(i);
        else if (a=="--basis")   basis=need(i);
        else if (a=="--acc")     acc=need(i);
        else if (a=="--accel")   accel=need(i);
        else if (a=="--maxiter") maxiter=std::stoi(need(i));
        else if (a=="--floor")   accj["floor"]=std::stod(need(i));
        else if (a=="--stall")   accj["stall"]=std::stoi(need(i));
        else if (a=="--trust")   accj["Trust"]=std::stod(need(i));
        else if (a=="--emax")    accj["EMax"]=std::stod(need(i));
        else if (a=="--nproj")   accj["NProj"]=(size_t)std::stoi(need(i));
        else { cout<<"unknown option "<<a<<endl; return 1; }
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
    //       NMaxIter   MinDeltaRo MinDelE MinVirial MinError StartingRelaxRo MergeTol verbose
    t->Iterate({(size_t)maxiter, Z*1e-4, 1e-5,  5e-1,    Z*2e-5,  0.5,            1e-7,    true});

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
