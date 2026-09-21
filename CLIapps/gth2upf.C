// File: CLIapps/gth2upf.C  Write a Quantum ESPRESSO UPF (v2.0.1) for one of OUR GTH pseudopotentials, with
// PP_CHI = OUR pseudo-atom's orbitals -- so pw.x / hp.x run the same pseudopotential AND the same Hubbard
// projector (QE's `atomic` = the PP_CHI 3d) as a qchem +U run (doc/OpenWork.md step 5 increment 3: hp.x is
// the U oracle, and an oracle is only as good as its match).
//
//   gth2upf --element Mn [--q 7] [--functional LDA] [--out Mn.pz-gth-q7.UPF] [--npool 16 --emin 0.05 --emax 200]
//
// What goes in:  PP_LOCAL   = 2 (V_long + V_short)(r)                       [Ry]    the analytic GTH local part
//                PP_BETA.p  = r * p_p(r), the KB-diagonalised GTH projectors   (normalised radials)
//                PP_DIJ     = diag(2 D_p)                                     [Ry]
//                PP_CHI.i   = r * chi_i(r), the pseudo-atom's occupied orbitals (LDA, unpolarized, a 16-exponent
//                             even-tempered pool per l -- the SAME atom that reproduces CP2K's ATOM code to mHa)
//                PP_RHOATOM = Sum_i f_i (r chi_i)^2
// on QE's logarithmic mesh r_i = exp(xmin + i dx)/Z (xmin = -7, dx = 0.0125, the HGH files' own).
// The functional label is LDA (SLA PZ NOGX NOGC): our atom is VWN, the GTH-PADE PP is Teter-Pade -- a few mHa on
// a total energy, nothing on U.

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

import qchem.AtomCalculation;                  // the pseudo-atom (PP_CHI, PP_RHOATOM)
import qchem.SCFParams;
import qchem.Pseudopotential.GTH_Potentials;   // GetGTH -> HGH_LocalPotential / HGH_SeparablePotential
import qchem.PeriodicTable;
import qchem.Symmetry.Atom.Spherical;          // AtomicSymmetry::Getl
import qchem.Symmetry.Irrep;
import qchem.Orbitals;
import qchem.Types;

using namespace qchem;

namespace
{
struct Wfc { int l; int n; double occ; double eps; std::vector<double> c; };   // c over the pool exponents

std::vector<double> Pool(int n, double emin, double emax)
{ std::vector<double> e; for (int i=0;i<n;i++) e.push_back(emin*std::pow(emax/emin, double(i)/(n-1))); return e; }

//! Unit-normalised r^l e^{-a r^2}: N^2 = 1 / int r^{2l+2} e^{-2 a r^2} dr = 2 (2a)^{l+3/2} / Gamma(l+3/2).
double RadialNorm(double a, int l) { return std::sqrt(2.0*std::pow(2.0*a, l+1.5)/std::tgamma(l+1.5)); }

int Period(int Z) { return Z<=2?1 : Z<=10?2 : Z<=18?3 : Z<=36?4 : Z<=54?5 : Z<=86?6 : 7; }

void WriteArray(std::ostream& os, const std::vector<double>& v)
{
    os<<std::scientific<<std::setprecision(15);
    for (size_t i=0;i<v.size();i++) { os<<std::setw(23)<<v[i]; if (i%4==3 || i+1==v.size()) os<<'\n'; else os<<' '; }
}
std::string E(double x) { std::ostringstream o; o<<std::scientific<<std::setprecision(15)<<x; return o.str(); }
}

int main(int argc, char** argv)
{
    std::string element, functional="LDA", out;
    int q=0, npool=16, lmaxPool=-1; double emin=0.05, emax=200.0;
    for (int i=1;i<argc;i++)
    {
        std::string a=argv[i]; auto val=[&]{ if (i+1>=argc) { std::cerr<<"missing value for "<<a<<"\n"; std::exit(1); } return std::string(argv[++i]); };
        if      (a=="--element")    element=val();
        else if (a=="--q")          q=std::stoi(val());
        else if (a=="--functional") functional=val();
        else if (a=="--out")        out=val();
        else if (a=="--npool")      npool=std::stoi(val());
        else if (a=="--emin")       emin=std::stod(val());
        else if (a=="--emax")       emax=std::stod(val());
        else if (a=="--lmax")       lmaxPool=std::stoi(val());
        else { std::cerr<<"usage: gth2upf --element <sym> [--q n] [--functional LDA] [--out file] [--npool 16 --emin 0.05 --emax 200 --lmax l]\n"; return a=="-h"||a=="--help" ? 0 : 1; }
    }
    if (element.empty()) { std::cerr<<"gth2upf: --element is required\n"; return 1; }
    const int Z=int(thePeriodicTable().GetZ(element));
    const Pseudopotential::GTH_PP pp=Pseudopotential::GetGTH(element, functional, q);
    const int Zion=pp.zion;
    if (out.empty()) out=element+".pz-gth-q"+std::to_string(Zion)+".UPF";

    // --- the pseudo-atom in the pool (every l up to 2; the EC occupies what the charge state fills) ---
    AtomCalcOptions o;
    o.type=AtomType::Gaussian; o.pseudopotential=true; o.valence=Zion;
    const std::vector<double> pool=Pool(npool, emin, emax);
    // The pool's l range: the occupied l's of the valence configuration (an empty extra channel only costs
    // the SCF its convergence): d-block from Sc on carries l=2, the p-block l=1, H/He/alkali l=0.
    if (lmaxPool<0) lmaxPool = (Z>=21 && Z<=30) || (Z>=39 && Z<=48) || (Z>=57 && Z<=80) ? 2 : (Z<=2 || Z==3 || Z==11 || Z==19 || Z==37 || Z==55) ? 0 : 1;
    for (int l=0;l<=lmaxPool;l++) o.exponentsByL.push_back({l, pool});
    SCFParams p; p.MinVirial=1e30; p.NMaxIter=200; p.Verbose=std::getenv("GTH2UPF_VERBOSE")!=nullptr;
    AtomCalculation atom(Z, Z-Zion, o, p);
    if (!atom.IsConverged()) { std::cerr<<"gth2upf: the pseudo-atom did not converge\n"; return 2; }
    std::vector<Wfc> wfcs;
    for (const Irrep& ir : atom.GetIrreps(Spin::None))
    {
        const auto* as=dynamic_cast<const Symmetry::Atom::AtomicSymmetry*>(ir.sym.get());
        const auto* os=dynamic_cast<const qchem::Orbitals::TOrbitals<double>*>(atom.Orbitals(ir));
        if (!as || !os) continue;
        const int l=int(as->Getl());
        // An open shell in an UNPOLARIZED spherical atom may be filled unevenly across m (O 2p^4 comes as a
        // 2/2 and a 2/4 entry at the same eigenvalue): one RADIAL per (l, eigenvalue), occupations summed.
        for (const auto* orb : os->template Iterate<qchem::Orbitals::TOrbital<double>>())
        {
            if (!orb->IsOccupied()) continue;
            bool merged=false;
            for (Wfc& w : wfcs) if (w.l==l && std::abs(w.eps-orb->GetEigenEnergy())<1e-6) { w.occ+=orb->GetOccupation(); merged=true; break; }
            if (merged) continue;
            int k=0; for (const Wfc& w : wfcs) if (w.l==l) k++;
            Wfc w; w.l=l; w.n=Period(Z)-(l==2?1:l==3?2:0)+k; w.occ=orb->GetOccupation(); w.eps=orb->GetEigenEnergy();
            const vec_t<double>& c=orb->GetCoeff(); w.c.assign(c.size(),0.0); for (size_t i=0;i<c.size();i++) w.c[i]=c[i];
            wfcs.push_back(w);
        }
    }
    // --- the mesh: QE's HGH files' own ---
    const double xmin=-7.0, dx=0.0125, zmesh=double(Z);
    std::vector<double> r, rab;
    for (int i=0;; i++) { const double x=xmin+i*dx, ri=std::exp(x)/zmesh; if (ri>100.0) break; r.push_back(ri); rab.push_back(ri*dx); }
    const size_t N=r.size();
    // --- local, projectors, wavefunctions, density on the mesh ---
    std::vector<double> vloc(N);
    for (size_t i=0;i<N;i++) vloc[i]=2.0*(pp.local.VlocLong(Z, r[i])+pp.local.VlocShort(Z, r[i]));
    const size_t nproj=pp.nonlocal.Count(Z);
    std::vector<std::vector<double>> beta(nproj, std::vector<double>(N)); std::vector<size_t> kkbeta(nproj); std::vector<int> lbeta(nproj);
    int lmax=0;
    for (size_t pj=0;pj<nproj;pj++)
    {
        lbeta[pj]=pp.nonlocal.L(Z,pj); lmax=std::max(lmax,lbeta[pj]);
        size_t last=0;
        for (size_t i=0;i<N;i++) { beta[pj][i]=r[i]*pp.nonlocal.RadialR(Z,pj,r[i]); if (std::abs(beta[pj][i])>1e-12) last=i; }
        kkbeta[pj]=std::min(N-1, last+1);
    }
    std::vector<std::vector<double>> chi(wfcs.size(), std::vector<double>(N)); std::vector<double> rho(N,0.0);
    for (size_t w=0;w<wfcs.size();w++)
    {
        for (size_t i=0;i<N;i++)
        {
            double s=0.0; for (size_t k=0;k<pool.size();k++) s+=wfcs[w].c[k]*RadialNorm(pool[k],wfcs[w].l)*std::pow(r[i],wfcs[w].l)*std::exp(-pool[k]*r[i]*r[i]);
            chi[w][i]=r[i]*s;
            rho[i]+=wfcs[w].occ*chi[w][i]*chi[w][i];
        }
    }
    // --- write ---
    std::ofstream f(out);
    f<<"<UPF version=\"2.0.1\">\n<PP_INFO>\n"
     <<"Generated by qchem gth2upf from the GTH/HGH parameter database (CP2K GTH_POTENTIALS transcoded):\n"
     <<"analytic local + separable-nonlocal parts on QE's mesh; PP_CHI are the GTH pseudo-atom's LDA orbitals\n"
     <<"(spherical, unpolarized, "<<npool<<"-exponent even-tempered Gaussian pool per l) -- the +U projector\n"
     <<"qchem uses (HubbardU_Atomic).  Element "<<element<<", q"<<Zion<<", functional "<<functional<<".\n"
     <<"Pseudo-atom E = "<<std::setprecision(8)<<atom.Energy()<<" Ha; orbitals:";
    for (const Wfc& w : wfcs) f<<"  "<<w.n<<"spdf"[w.l]<<" occ "<<w.occ<<" eps "<<std::setprecision(6)<<w.eps<<" Ha";
    f<<"\n</PP_INFO>\n<!--                               -->\n<!-- END OF HUMAN READABLE SECTION -->\n<!--                               -->\n";
    f<<"<PP_HEADER generated=\"qchem gth2upf\"\nauthor=\"Goedecker/Teter/Hutter parameters via qchem\"\ndate=\"2026\"\n"
     <<"comment=\"PP_CHI from the qchem GTH pseudo-atom\"\nelement=\""<<std::setw(2)<<element<<"\"\npseudo_type=\"NC\"\nrelativistic=\"scalar\"\n"
     <<"is_ultrasoft=\"F\"\nis_paw=\"F\"\nis_coulomb=\"F\"\nhas_so=\"F\"\nhas_wfc=\"F\"\nhas_gipaw=\"F\"\npaw_as_gipaw=\"F\"\ncore_correction=\"F\"\n"
     <<"functional=\"SLA PZ NOGX NOGC\"\nz_valence=\""<<E(double(Zion))<<"\"\ntotal_psenergy=\""<<E(2.0*atom.Energy())<<"\"\n"
     <<"wfc_cutoff=\""<<E(0.0)<<"\"\nrho_cutoff=\""<<E(0.0)<<"\"\nl_max=\""<<lmax<<"\"\nl_max_rho=\""<<2*lmax<<"\"\nl_local=\"-3\"\n"
     <<"mesh_size=\""<<N<<"\"\nnumber_of_wfc=\""<<wfcs.size()<<"\"\nnumber_of_proj=\""<<nproj<<"\"/>\n";
    f<<"<PP_MESH dx=\""<<E(dx)<<"\" mesh=\""<<N<<"\" xmin=\""<<E(xmin)<<"\" rmax=\""<<E(r.back())<<"\"\nzmesh=\""<<E(zmesh)<<"\">\n";
    f<<"<PP_R type=\"real\" size=\""<<N<<"\" columns=\"4\">\n";   WriteArray(f, r);   f<<"</PP_R>\n";
    f<<"<PP_RAB type=\"real\" size=\""<<N<<"\" columns=\"4\">\n"; WriteArray(f, rab); f<<"</PP_RAB>\n</PP_MESH>\n";
    f<<"<PP_LOCAL type=\"real\" size=\""<<N<<"\" columns=\"4\">\n"; WriteArray(f, vloc); f<<"</PP_LOCAL>\n";
    f<<"<PP_NONLOCAL>\n";
    for (size_t pj=0;pj<nproj;pj++)
    {
        f<<"<PP_BETA."<<pj+1<<" type=\"real\" size=\""<<N<<"\" columns=\"4\" index=\""<<pj+1<<"\" label=\""<<"spdf"[lbeta[pj]]<<pj+1<<"\" angular_momentum=\""<<lbeta[pj]
         <<"\" cutoff_radius_index=\""<<kkbeta[pj]+1<<"\"\ncutoff_radius=\""<<E(r[kkbeta[pj]])<<"\" ultrasoft_cutoff_radius=\""<<E(0.0)<<"\">\n";
        WriteArray(f, beta[pj]); f<<"</PP_BETA."<<pj+1<<">\n";
    }
    std::vector<double> dij(nproj*nproj, 0.0);
    for (size_t pj=0;pj<nproj;pj++) dij[pj*nproj+pj]=2.0*pp.nonlocal.Weight(Z,pj);
    f<<"<PP_DIJ type=\"real\" size=\""<<nproj*nproj<<"\" columns=\"4\">\n"; WriteArray(f, dij); f<<"</PP_DIJ>\n</PP_NONLOCAL>\n";
    f<<"<PP_PSWFC>\n";
    for (size_t w=0;w<wfcs.size();w++)
    {
        f<<"<PP_CHI."<<w+1<<" type=\"real\" size=\""<<N<<"\" columns=\"4\" index=\""<<w+1<<"\" label=\""<<wfcs[w].n<<"SPDF"[wfcs[w].l]<<"\" l=\""<<wfcs[w].l
         <<"\" occupation=\""<<E(wfcs[w].occ)<<"\" n=\""<<wfcs[w].n<<"\"\npseudo_energy=\""<<E(2.0*wfcs[w].eps)<<"\" cutoff_radius=\""<<E(0.0)<<"\" ultrasoft_cutoff_radius=\""<<E(0.0)<<"\">\n";
        WriteArray(f, chi[w]); f<<"</PP_CHI."<<w+1<<">\n";
    }
    f<<"</PP_PSWFC>\n";
    f<<"<PP_RHOATOM type=\"real\" size=\""<<N<<"\" columns=\"4\">\n"; WriteArray(f, rho); f<<"</PP_RHOATOM>\n</UPF>\n";
    // --- say what was written ---
    double nrho=0.0; for (size_t i=0;i<N;i++) nrho+=rho[i]*rab[i];
    std::cout<<"gth2upf: "<<out<<"  "<<element<<" q"<<Zion<<" "<<functional<<"  mesh "<<N<<"  proj "<<nproj<<" (l_max "<<lmax<<")  wfc";
    for (const Wfc& w : wfcs) std::cout<<" "<<w.n<<"spdf"[w.l]<<"("<<w.occ<<", eps "<<std::setprecision(5)<<w.eps<<" Ha)";
    std::cout<<"  int rhoatom = "<<std::setprecision(6)<<nrho<<" (expect "<<Zion<<")  E_atom = "<<std::setprecision(8)<<atom.Energy()<<" Ha"<<std::endl;
    return 0;
}
