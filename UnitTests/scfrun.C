// File: UnitTests/scfrun.C  Command-line SCF driver for tuning accelerators without recompiling.
//
// Example:
//   scfrun --Z 5 --model HF  --pol P --basis Slater     --acc Low    --accel Ladder --floor 1e-4 --stall 5
//   scfrun --Z 10 --model DHF --pol U --basis Slater_RKB --acc Medium --accel GDM    --emax 1e-3
//
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <algorithm>
#include <nlohmann/json.hpp>

import qchem.BasisSet.Atom.Factory;
import qchem.Unittests.QchemTester;       // QchemTester, TestAtom, TestDiracAtom, BasisSetAccuracy, DM_CD
import qchem.Hamiltonian.Factory;         // Model, Pol, Factory
import qchem.Hamiltonian.Internal.Hamiltonians;  // Ham_DFTcorr_U (real LSDA: Dirac exchange + VWN5), Ham_PP_U
import qchem.Structure;                    // Atom (the pseudo-ion for --model PP)
import qchem.ElectronConfiguration.AtomNR;  // PseudoAtom_EC (valence-only atomic config for the pseudo-atom)
import qchem.Math;                        // Pi
import qchem.Types;                       // rvec3_t
import qchem.ScalarFunction;              // ScalarFunction<double> (density / orbital evaluation)
import qchem.Mesh.Quadrature;             // qcMesh::RadialMesh, MakeRadial, Integrate (radial quadrature)

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
    virtual Hamiltonian* GetHamiltonian(st_t& c) const override { return Factory(m,p,c); }
};
class CliDiracAtom : public TestDiracAtom
{
    Model m; Pol p;
public:
    CliDiracAtom(int Z,int q,Model _m,Pol _p) : TestDiracAtom(Z,q), m(_m), p(_p) {}
    virtual Hamiltonian* GetHamiltonian(st_t& c) const override { return Factory(m,p,c); }
};
// DFT atom for generating the SAD atomic-density file.  --xc LDA = real LSDA (Dirac exchange + VWN5
// correlation, Ham_DFTcorr_U); --xc Xalpha = Slater X-alpha (alpha from --alpha, else the per-Z optimized
// value).  Unpolarized: an open-shell atom comes out spherical-ish (fractional shell occupation), which is
// exactly the spherically-averaged density a SAD seed wants.
class CliDFTAtom : public TestAtom
{
    bool lda; double alpha;
public:
    CliDFTAtom(int Z,int q,bool _lda,double _alpha) : TestAtom(Z,q), lda(_lda), alpha(_alpha) {}
    virtual Hamiltonian* GetHamiltonian(st_t& c) const override
    {
        if (lda) return new Ham_DFTcorr_U(c, GetMeshParams(), itsBasisSet);
        double a = alpha>0 ? alpha : QchemTester::itsPT.GetSlaterAlpha(GetZ());
        return Factory(Pol::UnPolarized, c, a, GetMeshParams(), itsBasisSet);
    }
};

// LOCAL-pseudopotential pseudo-atom for --model PP.  Build the FULL element's orbital basis (TestAtom(Z,0) so
// the basis factory sees the real Z), then SWAP in (a) aufbau over only the VALENCE electrons (Molecule_EC), and
// (b) the pseudo-ion structure (Atom with charge Z-valence, i.e. `valence` electrons), so the valence states
// fill the right angular channels over the pseudized (core-free) spectrum.  Hamiltonian = Ham_PP_U: kinetic +
// V_loc(r) (the GTH local model, mesh-quadratured) + Hartree + Dirac exchange + VWN5.  No Ven, no ion-ion.
class CliPPAtom : public TestAtom
{
    string element; int valence;
public:
    // TestAtom(Z, Z-val): itsZ = val, so the orbital-basis factory builds a VALENCE-range exponent pool
    // (no tight core functions -- those collapse/linear-depend under the soft pseudopotential), and the
    // structure is the pseudo-ion Atom(Z, Z-val) with exactly `val` electrons.  Only the electron config
    // needs swapping: PseudoAtom_EC(Z) gives the valence shells (e.g. Si s^2 p^2) over the atomic irreps.
    CliPPAtom(int Z, int val) : TestAtom(Z,Z-val), element(QchemTester::itsPT.GetSymbol(Z)), valence(val)
    {
        delete itsEC;
        itsEC = new PseudoAtom_EC(Z);
    }
    virtual Hamiltonian* GetHamiltonian(st_t& c) const override
    {
        return new Ham_PP_U(c, element, valence, GetMeshParams(), itsBasisSet);    // GTH local model for `element`
    }
};

// A valence density rho_val(r) = sum over the outermost orbitals of occ*|phi(r)|^2.  Built from a converged
// spherical all-electron atom (the round-hole solver) by keeping only the top-Nval-electron orbitals -- the
// pseudo-VALENCE source for the plane-wave SAD seed without a (square-peg) cubic plane-wave box.  The core
// region is too peaked vs a true pseudo density, but a seed only needs the bonding region in the ballpark.
class ValenceDensity : public virtual ScalarFunction<double>
{
    std::vector<std::pair<double,const ScalarFunction<double>*>> itsOrbs;   // (occupation, phi) -- phi owned by the WF
public:
    void Add(double occ, const ScalarFunction<double>* phi) { itsOrbs.push_back({occ,phi}); }
    virtual double  operator()(const rvec3_t& r) const
    { double rho=0; for (const auto& [occ,phi]:itsOrbs){ double v=(*phi)(r); rho+=occ*v*v; } return rho; }
    virtual rvec3_t Gradient(const rvec3_t&) const { return rvec3_t(0,0,0); }   // unused by the radial sampler
};

// Take the converged atom's occupied orbitals, sort by energy (descending), and accumulate whole orbitals
// until \a Nval electrons -- i.e. the outermost (valence) shells.
static ValenceDensity BuildValenceDensity(const QchemTester& t, int Nval)
{
    std::vector<std::tuple<double,double,const ScalarFunction<double>*>> occ;   // (energy, occupation, phi)
    for (const Irrep& irr : t.GetIrreps(Spin::None))
        for (const Orbital* o : t.GetOrbitals(irr)->Iterate())
            if (o->IsOccupied())
            {
                const auto* phi = dynamic_cast<const ScalarFunction<double>*>(o);
                assert(phi && "ValenceDensity: orbital is not a ScalarFunction");
                occ.push_back({o->GetEigenEnergy(), o->GetOccupation(), phi});
            }
    std::sort(occ.begin(), occ.end(), [](const auto& a, const auto& b){ return std::get<0>(a) > std::get<0>(b); });
    ValenceDensity v; double cum=0;
    for (const auto& [e,oc,phi] : occ) { if (cum >= Nval-1e-6) break; v.Add(oc,phi); cum+=oc; }
    return v;
}

// Sample a density rho(r) on a log-radial grid (spherically averaged over a fixed direction set), then MERGE
// the element's entry into the JSON database at `out` (saito.json-style array, one object per element x
// functional).  Prints the integral 4*pi*int r^2 rho dr as a charge sanity check.  `extra` adds fields to
// the entry (e.g. kind=valence).
static void DumpRadialDensity(const ScalarFunction<double>& cd, int Z, int Nelec, const string& functional,
                              const string& out, double rmin, double rmax, int ngrid,
                              const nlohmann::json& extra = nlohmann::json::object())
{
    // 26 directions (cube faces+edges+corners), normalized -> a cheap spherical average.
    std::vector<rvec3_t> dirs;
    const double b=1.0/std::sqrt(2.0), c=1.0/std::sqrt(3.0);
    const int s[]={-1,1};
    for (int x:s){ dirs.push_back(rvec3_t(x,0,0)); dirs.push_back(rvec3_t(0,x,0)); dirs.push_back(rvec3_t(0,0,x)); }
    for (int x:s) for (int y:s){ dirs.push_back(rvec3_t(b*x,b*y,0)); dirs.push_back(rvec3_t(b*x,0,b*y)); dirs.push_back(rvec3_t(0,b*x,b*y)); }
    for (int x:s) for (int y:s) for (int z:s) dirs.push_back(rvec3_t(c*x,c*y,c*z));

    // Sample on the shared Log radial mesh (nodes r_i = rmin*q^i; weights fold in the r^2 jacobian).
    qcMesh::RadialMesh mesh = qcMesh::MakeRadial({.radial=qcMesh::RadialKind::Log, .nRadial=ngrid,
                                                  .logStart=rmin, .logStop=rmax});
    const rvec_t& R=mesh.R();
    std::vector<double> rho(ngrid);
    for (int i=0;i<ngrid;i++)
    {
        double avg=0;
        for (const auto& u:dirs) avg += cd(u*R[i]);
        rho[i]=avg/dirs.size();
    }
    double charge = 4.0*Pi*qcMesh::Integrate(mesh, rvec_t(rho.size(), rho.data()));  // 4*pi*int r^2 rho dr

    nlohmann::json entry = {
        {"Z", Z}, {"symbol", QchemTester::itsPT.GetSymbol(Z)}, {"Nelec", Nelec},
        {"functional", functional},
        {"grid", {{"kind","log"},{"rmin",rmin},{"rmax",rmax},{"N",ngrid}}},
        {"charge", charge},          // 4*pi*int r^2 rho dr (~ Nelec)
        {"rho", rho},
    };
    for (const auto& [k,v] : extra.items()) entry[k]=v;

    nlohmann::json db = nlohmann::json::array();
    { std::ifstream in(out); if (in) in >> db; }       // load existing database (if any)
    if (!db.is_array()) db = nlohmann::json::array();
    // replace any existing entry for this (Z, functional)
    for (auto it=db.begin(); it!=db.end(); )
        if ((*it).value("Z",-1)==Z && (*it).value("functional",string())==functional) it=db.erase(it);
        else ++it;
    db.push_back(entry);
    std::sort(db.begin(), db.end(), [](const nlohmann::json& a, const nlohmann::json& b2)
              { return a.value("Z",0) < b2.value("Z",0); });
    std::ofstream(out) << db.dump(1,'\t') << endl;

    cout << "density dump: " << functional << " Z=" << Z << " -> " << out
         << "  (grid log " << rmin << ".." << rmax << " N=" << ngrid
         << ", charge=" << std::setprecision(6) << charge << " vs Nelec=" << Nelec << ")" << endl;
}

int main(int argc, char** argv)
{
    // ---- defaults ----
    int    Z=2, q=0, maxiter=50;
    string model="HF", pol="U", basis="", acc="Low", accel="DIIS";
    nlohmann::json accj;   // accelerator config passed to QchemTester
    // SCF convergence criteria (-1 => Z-scaled default); raise precision with --minfd/--virial.
    double minro=-1, minde=1e-5, virial=5e-1, minfd=-1, relax=0.5;
    // SAD atomic-density generation (DFT models only): dump rho(r) on a log grid to a JSON database.
    string out=""; double rmin=1e-4, rmax=20.0, alpha=-1; int ngrid=400;
    int valence=0;   // >0: dump only the outermost `valence` electrons' density (pseudo-VALENCE for PW SAD)

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
        "  --model <name>     HF | DHF | E1 | DE1 | LDA | Xalpha | PP  (default HF)\n"
        "                       (LDA/Xalpha are atom-DFT, for SAD density generation;\n"
        "                        PP = local-pseudopotential pseudo-atom, needs --valence)\n"
        "  --alpha <float>    Xalpha exchange parameter (default: per-Z optimized)\n"
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
        " SAD atomic-density generation (LDA/Xalpha models only):\n"
        "  --out <file>       merge rho(r) for this element into JSON database <file>\n"
        "  --rmin <float>     log-grid inner radius (bohr)            (default 1e-4)\n"
        "  --rmax <float>     log-grid outer radius (bohr)            (default 20)\n"
        "  --ngrid <int>      number of log-grid points               (default 400)\n"
        "  --valence <int>    dump only the outermost <n> electrons' density (pseudo-valence\n"
        "                       for the plane-wave SAD seed); 0 => full all-electron (default 0)\n"
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
        else if (a=="--alpha")   alpha=std::stod(need(i));
        else if (a=="--out")     out=need(i);
        else if (a=="--rmin")    rmin=std::stod(need(i));
        else if (a=="--rmax")    rmax=std::stod(need(i));
        else if (a=="--ngrid")   ngrid=std::stoi(need(i));
        else if (a=="--valence") valence=std::stoi(need(i));
        else { cout<<"unknown option "<<argv[i]<<endl; usage(cout); return 1; }
    }
    bool dirac = (model=="DHF" || model=="DE1");
    bool dft   = (model=="LDA" || model=="DFT" || model=="Xalpha");
    bool ppmodel = (model=="PP"); // local-pseudopotential pseudo-atom (Ham_PP_U); --valence = valence e- count
    if (basis.empty()) basis = dirac ? "Slater_RKB" : "Slater";
    accj["type"]=accel;

    // ---- string -> enum maps ----
    std::map<string,Model> models={{"HF",Model::HF},{"DHF",Model::DHF},{"E1",Model::E1},{"DE1",Model::DE1}};
    Pol pp = (pol=="P"||pol=="Polarized") ? Pol::Polarized : Pol::UnPolarized;
    using BT=BasisSet::Atom::Type;
    std::map<string,BT> bases={{"Slater",BT::Slater},{"Gaussian",BT::Gaussian},{"BSpline6",BT::BSpline6},
                               {"BSpliner6",BT::BSpliner6},{"Slater_RKB",BT::Slater_RKB},{"Gaussian_RKB",BT::Gaussian_RKB}};
    std::map<string,BasisSetAccuracy> accs={{"N3",BasisSetAccuracy::N3},{"N5",BasisSetAccuracy::N5},
                               {"Low",BasisSetAccuracy::Low},{"Medium",BasisSetAccuracy::Medium},{"High",BasisSetAccuracy::High}};
    // DFT models (LDA/Xalpha) and PP bypass the Model enum (they build the Hamiltonian directly, not Factory(Model,...)).
    if ((!dft && !ppmodel && !models.count(model))||!bases.count(basis)||!accs.count(acc)){cout<<"bad model/basis/acc"<<endl;return 1;}
    if (ppmodel && valence<=0){cout<<"--model PP requires --valence <n> (valence electron count)"<<endl;return 1;}

    cout << "scfrun: Z="<<Z<<" q="<<q<<" model="<<model<<" pol="<<pol
         << " basis="<<basis<<" acc="<<acc<<" accel="<<accel<<" : "<<accj.dump()<<endl;

    // ---- run ----
    QchemTester* t = ppmodel ? (QchemTester*) new CliPPAtom   (Z,valence)
                   : dft     ? (QchemTester*) new CliDFTAtom  (Z,q,model!="Xalpha",alpha)
                   : dirac   ? (QchemTester*) new CliDiracAtom(Z,q,models[model],pp)
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
         << "  Nelec="<<t->TotalCharge()
         << "  converged="<<(t->Converged()?"yes":"no") << endl;
    if (model=="HF")  t->RelativeHFError();
    if (model=="DHF") t->RelativeDHFError();

    // ---- SAD density file (DFT models only) ----
    if (dft && !out.empty())
    {
        const string functional = (model=="Xalpha"?"Xalpha":"LDA");
        if (valence>0)   // pseudo-VALENCE: only the outermost `valence` electrons (for the plane-wave SAD seed)
        {
            ValenceDensity vd = BuildValenceDensity(*t, valence);
            DumpRadialDensity(vd, Z, valence, functional, out, rmin, rmax, ngrid, {{"kind","valence"}});
        }
        else             // full all-electron density (for the molecular SAD seed)
        {
            std::unique_ptr<qchem::ChargeDensity::DM_CD> cd(t->GetChargeDensity());
            DumpRadialDensity(*cd, Z, Z-q, functional, out, rmin, rmax, ngrid);
        }
    }

    delete t;
    return 0;
}
