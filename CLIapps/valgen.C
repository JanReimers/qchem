// File: CLIapps/valgen.C  Command-line valence-basis generator (GPWPlan1 step 1).
//
// A thin arg-parser over the qchem.ValenceBasisGen library: point it at an element + its GTH
// pseudopotential and a set of even-tempered per-l windows, and it runs the spherical pseudo-atom SCF to
// VALIDATE the window (single source of truth -- validate exactly what you emit) and prints the Gaussian94
// element block ready to drop into a BasisSetData/*.bsd file.  The sibling to scfrun.C: scfrun tunes SCF
// accelerators, valgen mints GPW valence bases.
//
// Examples (reproducing the ValenceBasisGen_UT recipes, then an Al metal basis):
//   valgen --element F  --electrons 8 --shell 0:8:0.12:40 --shell 1:6:0.14:12   (F-: q defaults to valence 7)
//   valgen --element Na --shell 0:5:0.03:2   --shell 1:2:0.05:0.3 --seed        (neutral Na, valence q1)
//   valgen --element Al --shell 0:6:0.06:4   --shell 1:6:0.05:3   --out al.bsd  (neutral Al, valence q3)
//
// The CLI dogfoods the qchem.Reporting framework for its console output: report::SetConsole points the
// incremental renderer at stdout, then valgen brackets its own run and Emits recipe/result/seed sections
// (the nested pseudo-atom SCF -- an AtomCalculation, which now opens its own run -- stays console-quiet
// below depth 1, so the CLI shows the generator's summary, not the raw SCF internals).
//
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <optional>
#include <sstream>
#include <exception>
#include <cmath>
#include <map>
#include <algorithm>
#include "nlohmann/json.hpp"      // read the recorded basis.exponents / basis.usage back out of the report

import qchem.ValenceBasisGen;   // ValenceBasisRecipe, GenerateValenceBasis/GenerateSeedDensity/AssembleBasisFile
import qchem.PeriodicTable;     // thePeriodicTable(): symbol <-> Z (for the recipe report)
import qchem.Reporting;         // report::SetConsole / Begin / EmitSection -- the console output framework
import qchem.Pseudopotential.GTH_Potentials;   // GetGTHValences (enumerate valence/semicore PP variants)

using namespace qchem;
using std::cout;
using std::endl;
using std::string;
namespace rpt = qchem::report;

// Parse a "l:N:emin:emax" shell spec into an even-tempered window (l, exponents).  N exponents geometric
// from emin to emax inclusive.  Exits with a message on a malformed spec.
static std::pair<int,std::vector<double>> ParseShell(const string& spec)
{
    std::vector<string> f;
    std::stringstream ss(spec); string tok;
    while (std::getline(ss, tok, ':')) f.push_back(tok);
    if (f.size() != 4)
    {
        cout << "bad --shell \"" << spec << "\" (want l:N:emin:emax, e.g. 0:8:0.12:40)" << endl;
        exit(1);
    }
    const int    l    = std::stoi(f[0]);
    const int    N    = std::stoi(f[1]);
    const double emin = std::stod(f[2]);
    const double emax = std::stod(f[3]);
    return { l, EvenTemperedWindow(N, emin, emax) };
}

// A monospace heat-bar (ASCII '#') scaled to the largest population, so the eye catches which primitives
// carry the density.  Negatives (Mulliken overlap artefacts) clamp to empty.
static std::string HeatBar(double v, double vmax, int width=20)
{
    if (vmax<=0.0) return "";
    const int n=std::max(0,std::min(width,(int)std::lround(std::max(0.0,v)/vmax*width)));
    return std::string(n,'#') + std::string(width-n,' ');
}
static std::string g3(double x) { std::ostringstream os; os<<std::setprecision(3)<<x;                       return os.str(); }
static std::string e1(double x) { std::ostringstream os; os<<std::setprecision(1)<<std::scientific<<x;      return os.str(); }

// valgen's OWN heat-map renderer (the point of running at Normal): pull the basis's recorded exponent list and
// the WF's recorded populations out of the just-completed pseudo-atom run, JOIN them (index -> exponent), and
// render one grouped table [irrep, exponent, pop, bar].  The two blocks are recorded by different owners (the
// basis serializes exponents; the WF serializes populations) and are Verbose-only in the generic renderer, so
// at Normal this friendly joined view is all the user sees.
static void RenderUsageHeatMap()
{
    // Find the pseudo-atom run that carries both blocks (the GenerateValenceBasis AtomCalculation run).
    const rpt::json* run=nullptr;
    for (auto it=rpt::GlobalReport().begin(); it!=rpt::GlobalReport().end(); ++it)
        if (it.value().contains("basis") && it.value()["basis"].contains("usage")
                                         && it.value()["basis"].contains("exponents")) run=&it.value();
    if (!run) return;
    const rpt::json& basis=(*run)["basis"];

    std::map<string,std::vector<double>> exps;                          // irrep -> exponent list
    for (const auto& e : basis["exponents"])
        if (e.contains("irrep") && e.contains("values"))
            exps[e["irrep"].get<string>()] = e["values"].get<std::vector<double>>();

    double pmax=0.0;
    for (const auto& u : basis["usage"]) pmax=std::max(pmax, std::abs(u["pop"].get<double>()));

    rpt::json rows=rpt::json::array();
    for (const auto& u : basis["usage"])
    {
        const string irr=u["irrep"].get<string>();
        const size_t idx=u["index"].get<size_t>();
        const double pop=u["pop"].get<double>();
        const double e  =(exps.count(irr) && idx<exps[irr].size()) ? exps[irr][idx] : 0.0;
        rows.push_back(rpt::json{{"irrep",irr},{"exponent",e},{"pop",pop},{"bar",HeatBar(pop,pmax)}});
    }
    cout << "\nbasis usage  (occupation-weighted Mulliken population per primitive; Sum = electron count)\n";
    rpt::RenderConsole(rows, cout, rpt::Detail::Verbose);   // grouped table (rule between irreps), not gated here

    // ---- ADVICE (rules 1-3): per-shell verdicts, translated to the even-tempered knobs the user actually has
    //      (emin/emax/N/β) -- diffuse end must be well-populated; sharp low-pop is fine (shape fns); NEVER
    //      alternating ±.  Also surfaces the pivoted-Cholesky redundant list + the per-shell condition number. ----
    struct Shell { string irrep; std::vector<std::pair<double,double>> ep; double cond=0; std::vector<size_t> redund; };
    std::vector<Shell> sh;
    auto shellOf=[&](const string& irr)->Shell&{ for (auto& s:sh) if (s.irrep==irr) return s; sh.push_back({irr,{},0,{}}); return sh.back(); };
    for (const auto& u : basis["usage"])
    {
        const string irr=u["irrep"].get<string>(); const size_t idx=u["index"].get<size_t>();
        const double e=(exps.count(irr)&&idx<exps[irr].size())?exps[irr][idx]:0.0;
        shellOf(irr).ep.emplace_back(e, u["pop"].get<double>());
    }
    if (basis.contains("perIrrep")) for (const auto& pr:basis["perIrrep"])
        if (pr.contains("irrep")&&pr.contains("cond")) shellOf(pr["irrep"].get<string>()).cond=pr["cond"].get<double>();
    if (basis.contains("removed")) for (const auto& rm:basis["removed"])
        if (rm.contains("irrep")&&rm.contains("index")) shellOf(rm["irrep"].get<string>()).redund.push_back(rm["index"].get<size_t>());

    cout << "\nbasis usage - advice  (even-tempered knobs: raise beta = fewer N or wider emin..emax; shift emin for the diffuse end)\n";
    for (auto& s : sh)
    {
        std::sort(s.ep.begin(), s.ep.end());                         // ascending exponent: front=diffuse, back=sharp
        const size_t n=s.ep.size(); if (!n) continue;
        double pmax=0; for (auto&[e,p]:s.ep) pmax=std::max(pmax,std::abs(p));
        const double beta = n>1 ? std::pow(s.ep.back().first/s.ep.front().first, 1.0/(n-1)) : 0.0;
        std::vector<string> msg;
        // (1) diffuse end must be well-populated
        const double pdiff=s.ep.front().second;
        if (pmax>0 && std::abs(pdiff) < 0.10*pmax)
            msg.push_back("diffuse alpha="+g3(s.ep.front().first)+" barely populated -> RAISE emin (wasted function)");
        else if (pmax>0 && pdiff >= 0.85*pmax)
            msg.push_back("diffuse end IS the peak (alpha="+g3(s.ep.front().first)+") -> density spills past -> LOWER emin");
        // (3) overcompleteness: alternating O(1) pops, or a high condition number, or too-tight beta
        int alt=0; for (size_t i=1;i<n;i++)
            if (s.ep[i].second*s.ep[i-1].second<0 && std::abs(s.ep[i].second)>0.25*pmax && std::abs(s.ep[i-1].second)>0.25*pmax) alt++;
        if (alt>0)                    msg.push_back("pops ALTERNATE +/- -> overcomplete -> raise beta (fewer N or wider range)");
        else if (s.cond>1e3)          msg.push_back("cond="+e1(s.cond)+" high -> near-dependent -> raise beta");
        else if (beta>0 && beta<2.0)  msg.push_back("beta="+g3(beta)+" tight -> overcomplete risk");
        // sharp low-pop shape tail (informational: E needs them, but trimmable)
        size_t tail=0; for (size_t i=n; i-->0;) { if (pmax>0 && std::abs(s.ep[i].second)<0.10*pmax) tail++; else break; }
        if (tail>=1) msg.push_back(std::to_string(tail)+" sharp low-pop shape fn(s) (up to alpha="+g3(s.ep.back().first)
                                   +"): hold energy but carry little density -> trim emax/N for a leaner basis (costs a few mHa)");
        // (2) pivoted-Cholesky redundant list -> exponents (with the even-tempered caveat)
        if (!s.redund.empty())
        {
            string r; for (size_t idx:s.redund){ if(!r.empty())r+=", "; r+="alpha="+g3(idx<exps[s.irrep].size()?exps[s.irrep][idx]:0.0); }
            msg.push_back("pivoted-Cholesky flags {"+r+"} linearly REDUNDANT (even-tempered: can't drop one -> raise beta / cut N)");
        }
        cout << "  " << s.irrep << ":  N=" << n << " beta=" << g3(beta) << " cond=" << e1(s.cond);
        if (msg.empty()) cout << "   -- looks healthy.\n";
        else { cout << "\n"; for (auto& m:msg) cout << "      - " << m << "\n"; }
    }
}

int main(int argc, char** argv)
{
    // ---- defaults ----
    string element, functional = "LDA", out, seedOut, name, detail = "normal";
    int    q = 0, electrons = 0;                         // q 0 => resolve by policy below (valence / --semicore)
    bool   seed = false;                                 // also generate the seed valence density
    bool   semicore = false;                             // pick the first semicore variant instead of valence
    bool   iterations = false;                           // print the SCF per-iteration convergence trace
    bool   floor = false;                                // also compute the complete-basis (pool) energy floor
    int    ngrid = 400; double rmin = 1e-4, rmax = 20.0; // seed-density radial log grid
    std::vector<std::pair<int,std::vector<double>>> shells;

    // ---- help text ----
    auto usage = [&](std::ostream& os){
        os <<
        "Usage: valgen --element <sym> --shell <l:N:emin:emax> [--shell ...] [options]\n"
        "  Generate a GPW valence Gaussian basis for one element under its GTH pseudopotential.\n"
        "  The spherical pseudo-atom SCF validates the window, then the element block is printed.\n"
        "  Each option takes a value, written as \"--flag value\" or \"--flag=value\".\n"
        "\n"
        " Recipe:\n"
        "  --element <sym>    element symbol, e.g. F | Na | Al           (required)\n"
        "  --q <int>          which GTH pseudopotential variant (its 'q', the CP2K/GTH convention =\n"
        "                       the # valence electrons of the NEUTRAL atom / the core/valence split).\n"
        "                       Default = VALENCE (smallest q); most elements offer only one.  This is\n"
        "                       NOT the ion charge -- for an ion use --electrons.  (see --semicore)\n"
        "  --semicore         use the first SEMICORE variant (next q up) instead of valence, for the\n"
        "                       elements that offer one (e.g. Na, Ag); overridden by an explicit --q\n"
        "  --electrons <int>  how many valence electrons to OCCUPY = the CHARGE STATE; default = q\n"
        "                       (neutral).  THIS is where an ion goes: F- => --electrons 8 (q stays 7),\n"
        "                       a cation => fewer.  Does not change which PP is used (--q/--semicore do).\n"
        "  --functional <fn>  GTH parameter set                          (default LDA)\n"
        "  --shell <spec>     an even-tempered window l:N:emin:emax (REPEATABLE, >=1 required);\n"
        "                       N exponents geometric from emin to emax, e.g. --shell 0:8:0.12:40\n"
        "                       KEEP EXPONENTS DISJOINT ACROSS l (the reader merges shared exponents)\n"
        "\n"
        " Seed density (optional):\n"
        "  --seed             also generate the seed valence density (reports charge + <r> diffuseness)\n"
        "  --ngrid <int>      seed log-grid points                       (default 400)\n"
        "  --rmin <float>     seed log-grid inner radius (bohr)          (default 1e-4)\n"
        "  --rmax <float>     seed log-grid outer radius (bohr)          (default 20)\n"
        "\n"
        " Output:\n"
        "  --out <file>       write the assembled single-element .bsd to <file> (also printed to stdout)\n"
        "  --seed-out <file>  write the seed-density JSON library entry to <file> (implies --seed)\n"
        "  --name <title>     .bsd file title / BASIS= tag               (default: auto)\n"
        "  --detail <lvl>     console verbosity: terse | normal | verbose (default normal;\n"
        "                       verbose also dumps the raw basis.usage/basis.exponents blocks)\n"
        "  --iterations       print the pseudo-atom SCF's per-iteration convergence trace\n"
        "  --floor            also run the complete-basis (accuracy-pool) reference and report the\n"
        "                       energy floor + this window's gap above it (gap_mHa; ~how incomplete)\n"
        "\n"
        "  -h, --help         show this help and exit\n";
    };

    // ---- parse "--flag value" / "--flag=value" pairs (mirrors scfrun) ----
    for (int i = 1; i < argc; i++)
    {
        string a = argv[i], inlineVal; bool hasInline = false;
        auto eq = a.find('=');
        if (eq != string::npos) { inlineVal = a.substr(eq+1); a.erase(eq); hasInline = true; }
        auto need = [&](int& i){
            if (hasInline) return inlineVal;
            if (i+1 >= argc) { cout << "missing value for " << argv[i] << endl; exit(1); }
            return string(argv[++i]);
        };
        if      (a=="-h" || a=="--help") { usage(cout); return 0; }
        else if (a=="--element")    element    = need(i);
        else if (a=="--q")          q          = std::stoi(need(i));
        else if (a=="--semicore")   semicore   = true;
        else if (a=="--iterations") iterations = true;
        else if (a=="--floor")      floor      = true;
        else if (a=="--electrons")  electrons  = std::stoi(need(i));
        else if (a=="--functional") functional = need(i);
        else if (a=="--shell")      shells.push_back(ParseShell(need(i)));
        else if (a=="--seed")       seed       = true;
        else if (a=="--ngrid")      ngrid      = std::stoi(need(i));
        else if (a=="--rmin")       rmin       = std::stod(need(i));
        else if (a=="--rmax")       rmax       = std::stod(need(i));
        else if (a=="--out")        out        = need(i);
        else if (a=="--seed-out") { seedOut    = need(i); seed = true; }
        else if (a=="--name")       name       = need(i);
        else if (a=="--detail")     detail     = need(i);
        else { cout << "unknown option " << argv[i] << endl; usage(cout); return 1; }
    }

    if (element.empty()) { cout << "valgen: --element is required\n"; usage(cout); return 1; }
    if (shells.empty())  { cout << "valgen: at least one --shell l:N:emin:emax is required\n"; usage(cout); return 1; }

    // Resolve --element: accept a SYMBOL ("Al") or a numeric atomic number ("13").  Validate up front and
    // canonicalise to the table's symbol -- an unknown symbol / out-of-range Z otherwise resolves to Z=0 and
    // aborts deep in the basis build (GetSymbol(0) is out of bounds) with an opaque std::bad_alloc.
    const auto& pt = thePeriodicTable();
    int Z = 0;
    if (element.find_first_not_of("0123456789") == string::npos)      // all digits => atomic number
    {
        Z = std::stoi(element);
        if (Z < 1 || (size_t)Z > pt.GetNumElements())
        { cout << "valgen: --element " << Z << " is out of range (1.." << pt.GetNumElements() << ")\n"; return 1; }
        element = pt.GetSymbol(Z);                                     // canonical symbol for the recipe/report
    }
    else
    {
        Z = (int)pt.GetZ(element);
        if (Z == 0)
        { cout << "valgen: unknown element \"" << element << "\" (want a symbol like Al, or a Z like 13)\n"; return 1; }
    }

    // Resolve the PP variant q: an explicit --q wins; otherwise pick by POLICY from the tabulated variants --
    // default = VALENCE (smallest q), or the first SEMICORE (next q up) with --semicore.  (The GTH database's
    // own "default" is often the semicore, so a valence-first tool must select explicitly.)  If the element is
    // absent from the DB, leave q=0 and let GenerateValenceBasis report the clean "not in database" error.
    const std::vector<int> variants = qchem::Pseudopotential::GetGTHValences(element, functional);
    if (q <= 0 && !variants.empty())
        q = (semicore && variants.size() > 1) ? variants[1] : variants.front();

    ValenceBasisRecipe r;
    r.element    = element;
    r.Zion       = q;
    r.electrons  = electrons;
    r.functional = functional;
    r.shells     = shells;

    // Dogfood the reporting framework.  valgen runs at Normal: the raw basis.usage/basis.exponents blocks are
    // Verbose-only in the generic renderer, so at Normal they stay in the json but off the console -- valgen
    // instead renders its OWN joined heat map (exponent + population + bar) from those recorded blocks below.
    // --detail verbose additionally shows the raw generic blocks.
    const rpt::Detail lvl = detail=="terse"  ? rpt::Detail::Terse
                          : detail=="verbose"? rpt::Detail::Verbose
                                             : rpt::Detail::Normal;
    rpt::SetConsole(cout, lvl);

    GeneratedBasis g{};
    std::optional<GeneratedSeedDensity> s;
    try
    {
        // The recipe header is its own SHORT run so it renders at depth 1; then GenerateValenceBasis's
        // pseudo-atom SCF (an AtomCalculation) opens its OWN depth-1 run and renders scf/basis/usage/cache --
        // including the Verbose basis-usage heat map, the whole point of the tool.  (If valgen wrapped a run
        // AROUND the generate call, the SCF would nest at depth 2 and the heat map would stay console-quiet.)
        {
            rpt::Begin("valgen " + element + " recipe");
            struct RunScope { ~RunScope() { rpt::End(); } } runScope;
            rpt::Section recipe("recipe");
            rpt::Set("element",    element);
            rpt::Set("Z",          (long)Z);
            rpt::Set("q",          (long)q);               // the chosen PP variant (valence unless --semicore/--q)
            { std::ostringstream os; for (size_t i=0;i<variants.size();++i){ if(i) os<<","; os<<variants[i]; }
              rpt::Set("qAvailable", os.str()); }          // all tabulated variants (valence..semicore)
            rpt::Set("electrons",  (long)(electrons>0?electrons:q));   // occupancy; 0 => neutral (== q)
            rpt::Set("functional", functional);
            for (const auto& [l, es] : shells)
            {
                rpt::Row shell("shells");
                rpt::Set("l",     (long)l);
                rpt::Set("N",     (long)es.size());
                rpt::Set("emin",  es.empty() ? 0.0 : es.front());
                rpt::Set("emax",  es.empty() ? 0.0 : es.back());
            }
        }

        // -- validate: run the pseudo-atom SCF in exactly these shells (records basis.exponents + basis.usage) --
        rpt::Log("valgen: running the pseudo-atom SCF to validate the window...");
        g = GenerateValenceBasis(r, iterations);   // --iterations prints the per-iteration convergence trace
        RenderUsageHeatMap();      // valgen's joined exponent+population heat map (before seed adds a 2nd run)

        // -- result (+seed) as a closing run --
        {
            rpt::Begin("valgen " + element + " result");
            struct RunScope { ~RunScope() { rpt::End(); } } runScope;
            rpt::json result;
            result["energy"]     = g.energy;
            result["converged"]  = g.converged;
            long nExp = 0; for (const auto& [l, es] : shells) nExp += (long)es.size();
            result["nExponents"] = nExp;
            if (floor)   // the complete-basis (accuracy-pool) reference + how far this window sits above it
            {
                rpt::Log("valgen: computing the complete-basis energy floor (pool reference)...");
                const double ef = GenerateFloorEnergy(r);
                result["Efloor"]  = ef;                          // ~complete-basis PP energy
                result["gap_mHa"] = 1000.0 * (g.energy - ef);    // this window's incompleteness (positive = above floor)
            }
            rpt::EmitSection("result", result);

            // seed density (optional): the SAME pseudo-atom SCF, sampled to a spherical rho(r).  Its own SCF
            // nests here (depth 2, quiet) -- we only want the charge/⟨r⟩ metrics, not a second SCF report.
            if (seed)
            {
                rpt::Log("valgen: sampling the seed valence density...");
                s = GenerateSeedDensity(r, ngrid, rmin, rmax);
                rpt::json sj;
                sj["charge"]    = s->charge;
                sj["meanR"]     = s->meanR;     // <r> diffuseness (an anion's exceeds the neutral atom's)
                sj["converged"] = s->converged;
                sj["ngrid"]     = ngrid;
                sj["rmin"]      = rmin;
                sj["rmax"]      = rmax;
                rpt::EmitSection("seed", sj);
            }
        }
    }
    catch (const std::exception& e)
    {
        const string w = e.what();
        cout << "\nvalgen: " << w << "\n";
        if (w.find("angular momentum") != string::npos)         // the BasisSet_HF missing-l throw
            cout << "  Every angular momentum the element's valence configuration occupies needs its own --shell.\n"
                 << "  (the q valence electrons fill in aufbau order, e.g. Al q3 = 3s^2 3p^1 -> l=0 AND l=1.)\n";
        else                                                     // typically a too-diffuse exponent (Fock -> NaN)
            cout << "  The pseudo-atom SCF failed.  The usual cause is an exponent that is TOO DIFFUSE (a very\n"
                 << "  small emin, ~<=0.02): the density gets tiny over a large radius and the atomic XC grid\n"
                 << "  produces a non-finite value.  Try raising emin in the offending --shell.\n";
        return 1;
    }

    // ---- the PRODUCT: the Gaussian94 element block (print for capture; write files on request) ----
    cout << "\n===== BEGIN " << element << " valence block (Gaussian94 .bsd) =====\n"
         << g.block
         << "===== END " << element << " valence block =====" << endl;

    if (!out.empty())
    {
        const string title = name.empty()
            ? ("Low-q GTH-pseudopotential valence basis for " + element + " (valgen, doc/GPWPlan1 sec 1)")
            : name;
        const string file = AssembleBasisFile(title, { g.block });
        std::ofstream(out) << file;
        cout << "valgen: wrote " << out << " (" << element << " single-element .bsd)" << endl;
    }
    if (s && !seedOut.empty())
    {
        std::ofstream(seedOut) << s->jsonEntry << endl;
        cout << "valgen: wrote " << seedOut << " (seed valence-density library entry)" << endl;
    }
    // Always succeed: the ENERGY is the validation metric (reported in `result`), not the density gate --
    // a degenerate open-shell valence density need not settle (see qchem.ValenceBasisGen / doc/GPWPlan.md),
    // so `converged=false` is expected and NOT a failure.
    return 0;
}
