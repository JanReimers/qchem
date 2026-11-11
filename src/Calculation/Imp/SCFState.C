// File: Calculation/Imp/SCFState.C  The saved SCF state every kind shares: the base fingerprint, the comparison, the framing,
// the writer and reader, the Löwdin carry and the restore checks.  (Moved out of Imp/SolidState.C, 2026-11-09, V-FACADE-STATE.)
module;
#include <algorithm>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.SCFState;

import qchem.HDF5;
import qchem.BasisSet.AoShellSource;    // the orbital shells (the basis half of the fingerprint)
import qchem.Blaze;                     // blazem::eigen (the Löwdin carry's S^{±1/2})

namespace qchem
{

static const char* kFormat          = "qchem-state 2.0";         //!< 2.0 (2026-11-09): every kind; + kind, irrep, S
static const char* kFormatMajor     = "qchem-state 2";           //!< every 2.x is read: a minor version only ADDS optional fields
static const char* kSolidFormatMajor= "qchem-solid-state 1";     //!< the solid's 1.x (1.0, 1.1): read as kind "solid"

size_t StateBlock::n() const
{
    if (real) return real->GetNumFunctions();
    if (cplx) return cplx->GetNumFunctions();
    throw std::logic_error("StateBlock: a block with no basis");
}
std::string StateBlock::Label() const {return irrep.sym ? irrep.sym->GetLabel() : std::string();}

//---------------------------------------------------------------------------------------------------
//  THE BASE FINGERPRINT
//---------------------------------------------------------------------------------------------------
void AddStructure(StateFingerprint& fp, const Structure& st)
{
    std::vector<double> Z, R;
    st.ForEachSite([&](int z, const rvec3_t& r, bool){ Z.push_back(z); R.push_back(r.x); R.push_back(r.y); R.push_back(r.z); });
    fp.refuseNum["structure.Z"]=Z;
    fp.warmNum["structure.positions"]=R;      // WARM: a geometry step carries D (the basis centres follow)
}

static std::string SpinGroupName(SpinGroup g) {return g==SpinGroup::Polarized ? "Polarized" : "UnPolarized";}

void AddElectrons(StateFingerprint& fp, int Nelec, SpinGroup g, int multiplicity, const std::string& functional)
{
    fp.refuseNum["electrons.Nelec"]={double(Nelec)};
    fp.refuseStr["spinGroup"]=SpinGroupName(g);
    fp.warmStr["functional"]=functional;
    fp.warmNum["spin.multiplicity"]={double(multiplicity)};
}

void AddBlocks(StateFingerprint& fp, const std::vector<StateBlock>& blocks)
{
    std::vector<double> ms, n, real;
    std::string labels;
    for (const StateBlock& b : blocks)
    {
        ms.push_back(double(int(b.irrep.ms)));
        n.push_back(double(b.n()));
        real.push_back(b.real ? 1.0 : 0.0);
        labels += (labels.empty() ? "" : ",") + b.Label();
    }
    fp.refuseNum["blocks.ms"]=ms;
    fp.refuseNum["blocks.n"]=n;
    fp.refuseStr["blocks.irrep"]=labels;
    fp.warmNum  ["blocks.real"]=real;          // a working-type change promotes/demotes D
}

namespace
{
void AddShellList(StateFingerprint& fp, const std::vector<std::vector<Symmetry::Molecule::AoShell>>& lists)
{
    std::vector<double> L, nc, np, cen, ex, co;
    for (const auto& shells : lists)
        for (const auto& sh : shells)
        {
            L.push_back(sh.rep->L());
            nc.push_back(sh.nComponents());
            np.push_back(sh.exponents.size());
            cen.push_back(sh.center.x); cen.push_back(sh.center.y); cen.push_back(sh.center.z);
            for (size_t p=0;p<sh.exponents.size();p++)    ex.push_back(sh.exponents[p]);
            for (size_t p=0;p<sh.coefficients.size();p++) co.push_back(sh.coefficients[p]);
        }
    fp.refuseNum["basis.shellL"]=L;
    fp.refuseNum["basis.shellComponents"]=nc;
    fp.refuseNum["basis.shellPrimitives"]=np;
    fp.warmNum["basis.shellCenters"]=cen;
    fp.refuseNum["basis.exponents"]=ex;
    fp.refuseNum["basis.coefficients"]=co;
}
template <class T> std::vector<Symmetry::Molecule::AoShell> ShellsOf(const BasisSet::Orbital_1E_IBS<T>& b)
{
    const auto* src=dynamic_cast<const BasisSet::AoShellSource*>(&b);
    if (!src) throw std::logic_error("SCFState: a basis block does not report its AO shells (AoShellSource) -- a saved state "
                                     "cannot identify this basis");
    return src->GetAoShells();
}
} //anonymous namespace

void AddShells(StateFingerprint& fp, const std::vector<const BasisSet::Orbital_1E_IBS<double>*>& blocks)
{
    if (blocks.empty()) throw std::logic_error("SCFState::AddShells: no blocks");
    std::vector<std::vector<Symmetry::Molecule::AoShell>> lists;
    for (const auto* b : blocks) lists.push_back(ShellsOf(*b));
    AddShellList(fp, lists);
}
void AddShells(StateFingerprint& fp, const BasisSet::Orbital_1E_IBS<dcmplx>& block) {AddShellList(fp, {ShellsOf(block)});}

//---------------------------------------------------------------------------------------------------
//  THE COMPARISON
//---------------------------------------------------------------------------------------------------
namespace
{
std::string Show(const std::vector<double>& v)
{
    std::ostringstream os;
    os<<std::setprecision(10)<<"[";
    const size_t shown=std::min<size_t>(v.size(), 6);
    for (size_t k=0;k<shown;k++) os<<(k?",":"")<<v[k];
    if (v.size()>shown) os<<",... ("<<v.size()<<" values)";
    os<<"]";
    return os.str();
}
//! Equal to 1e-10 RELATIVE, with a 1e-14 absolute floor for entries at zero.  (Was 1e-10*max(1,|a|): an ABSOLUTE 1e-10 below 1,
//! so a tolerance recorded as a value -- screenEps 1e-10 vs 1e-12 -- compared as the same run.  Found 2026-10-30 by the gpw.* keys.)
bool Same(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size()!=b.size()) return false;
    for (size_t k=0;k<a.size();k++)
        if (std::fabs(a[k]-b[k]) > std::max(1e-14, 1e-10*std::max(std::fabs(a[k]), std::fabs(b[k])))) return false;
    return true;
}

//! WHAT a warm difference MEANS for the restart, keyed on the fingerprint name (exact, else the longest matching prefix ending
//! in '.'): the user's ask (2026-10-29) was the consequence, not only the two values.  An unknown key prints raw.
std::string Consequence(const std::string& key)
{
    static const std::string moves=": the converged ground state moves -- expect iterations to reach the new one";
    static const std::map<std::string, std::string> table=
    {
        {"functional",            "a different XC functional"+moves},
        {"hubbard.",              "the +U term changed"+moves},
        {"spin.multiplicity",     "a different spin state: the occupations are re-filled, expect iterations"},
        {"grids.densityEcut",     "the density grid changed"+moves},
        {"grids.cutoffFactor",    "the density grid's floor changed"+moves},
        {"grids.ladderFactor",    "the multigrid ladder changed"+moves},
        {"grids.raster",          "the density raster changed"+moves},
        {"grids.xcMesh",          "the XC quadrature changed"+moves},
        {"grids.beckeEps",        "the Becke partition tolerance changed"+moves},
        {"xc.fit",                "the v_xc fit basis changed"+moves},
        {"xc.",                   "the v_xc route changed"+moves},
        {"coulomb.",              "the Hartree route changed"+moves},
        {"gpw.",                  "a GPW numerical tolerance changed"+moves},
        {"policy.",               "a declared CP2K-deviation route changed the Hamiltonian"+moves},
        {"ortho.",                "the basis orthogonalisation changed: functions may be dropped or kept differently"+moves},
        {"scf.kT",                "the smearing changed: the occupations are re-filled at the new temperature"+moves},
        {"accuracy.eps",          "a different accuracy target: the meshes it sizes may differ"},
        {"blocks.real",           "a block's working type changed: D is promoted/demoted exactly, the physics does not move"},
        {"basis.shellCenters",    "the basis centres moved with the atoms"+moves},
    };
    if (auto it=table.find(key); it!=table.end()) return it->second;
    std::string best;
    for (const auto& [k,v] : table)
        if (k.back()=='.' && key.rfind(k,0)==0 && k.size()>best.size()) best=k;
    return best.empty() ? std::string() : table.at(best);
}

//! Every differing entry of one class, one line each.
std::vector<std::string> Differences(const std::map<std::string,std::vector<double>>& sn,
                                     const std::map<std::string,std::vector<double>>& nn,
                                     const std::map<std::string,std::string>& ss,
                                     const std::map<std::string,std::string>& ns)
{
    std::vector<std::string> out;
    std::set<std::string> keys;
    for (const auto& [k,v] : sn) keys.insert(k);
    for (const auto& [k,v] : nn) keys.insert(k);
    for (const std::string& k : keys)
    {
        auto a=sn.find(k), b=nn.find(k);
        if      (a==sn.end()) out.push_back(k+": not in the saved state (now "+Show(b->second)+")");
        else if (b==nn.end()) out.push_back(k+": saved "+Show(a->second)+", absent from this run");
        else if (!Same(a->second, b->second)) out.push_back(k+": saved "+Show(a->second)+", now "+Show(b->second));
    }
    keys.clear();
    for (const auto& [k,v] : ss) keys.insert(k);
    for (const auto& [k,v] : ns) keys.insert(k);
    for (const std::string& k : keys)
    {
        auto a=ss.find(k), b=ns.find(k);
        const std::string sa = a==ss.end() ? "(absent)" : "'"+a->second+"'";
        const std::string sb = b==ns.end() ? "(absent)" : "'"+b->second+"'";
        if (sa!=sb) out.push_back(k+": saved "+sa+", now "+sb);
    }
    for (std::string& l : out)
        if (const std::string c=Consequence(l.substr(0, l.find(':'))); !c.empty()) l+=" -- "+c;
    return out;
}

double Displacement(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size()!=b.size() || a.size()%3!=0) return -1.0;
    double m=0.0;
    for (size_t i=0;i<a.size();i+=3)
        m=std::max(m, std::sqrt((a[i]-b[i])*(a[i]-b[i]) + (a[i+1]-b[i+1])*(a[i+1]-b[i+1]) + (a[i+2]-b[i+2])*(a[i+2]-b[i+2])));
    return m;
}
} //anonymous namespace

double MaxDisplacement(const StateFingerprint& saved, const StateFingerprint& now)
{
    const auto sp=saved.warmNum.find("structure.positions"), np=now.warmNum.find("structure.positions");
    if (sp==saved.warmNum.end() || np==now.warmNum.end()) return 0.0;
    return std::max(0.0, Displacement(sp->second, np->second));
}

Outcome<std::vector<std::string>, RestartRefusal> CompareFingerprints(const StateFingerprint& savedAsWritten,
                                                                     const StateFingerprint& now, const FingerprintMigration& migrate)
{
    using O=Outcome<std::vector<std::string>, RestartRefusal>;
    const StateFingerprint saved=migrate ? migrate(savedAsWritten, now) : savedAsWritten;
    const std::vector<std::string> refuse=Differences(saved.refuseNum, now.refuseNum, saved.refuseStr, now.refuseStr);
    if (!refuse.empty())
    {
        std::string d="the saved state was made for a DIFFERENT run:";
        for (const std::string& r : refuse) d+="\n    "+r;
        return O::Fail({RestartRefusal::Why::Mismatch, d});
    }
    std::vector<std::string> warm=Differences(saved.warmNum, now.warmNum, saved.warmStr, now.warmStr);
    // A moved atom is ONE fact: say how far, and drop the basis-centre line that is only its consequence.
    if (const double d=MaxDisplacement(saved, now); d>1e-10)
    {
        std::erase_if(warm, [](const std::string& l){ return l.rfind("structure.positions:",0)==0 || l.rfind("basis.shellCenters:",0)==0; });
        std::ostringstream os;
        os<<"structure.positions: atoms moved, max displacement "<<std::setprecision(6)<<d<<" bohr (the basis follows)";
        warm.push_back(os.str());
    }
    return O::Ok(std::move(warm));
}

//---------------------------------------------------------------------------------------------------
//  THE BLOCK MATRICES (shared by the writer, the snapshot and the carry)
//---------------------------------------------------------------------------------------------------
namespace
{
//! A block's overlap, promoted to complex (a real block's exactly).
mat_t<dcmplx> BlockOverlap(const StateBlock& r)
{
    const size_t n=r.n();
    mat_t<dcmplx> S(n, n);
    if (r.real) { const auto& s=r.real->Overlap(); for (size_t a=0;a<n;a++) for (size_t b=0;b<n;b++) S(a,b)=dcmplx(s(a,b)); }
    else        { const auto& s=r.cplx->Overlap(); for (size_t a=0;a<n;a++) for (size_t b=0;b<n;b++) S(a,b)=s(a,b); }
    return S;
}

//! The orbitals of a block, of its basis's scalar.
template <class U> const Orbitals::TOrbitals<U>* Typed(const StateBlock& r)
{
    if (!r.orbitals) throw std::logic_error("SCFState: a block handed over without its orbitals");
    return dynamic_cast<const Orbitals::TOrbitals<U>*>(r.orbitals);
}

// One block's physical D = Σ_i f_i c_i c_i^† exactly as TOrbital::AddDensityMatrix, held complex (a real block's promotes exactly);
// with \a energyWeighted, W = Σ_i f_i ε_i c_i c_i^† (the force's and the stress's Pulay term).
template <class U> mat_t<dcmplx> BlockD(const Orbitals::TOrbitals<U>& os, size_t n, bool energyWeighted)
{
    mat_t<dcmplx> D(n,n);
    for (size_t a=0;a<n;a++) for (size_t b=0;b<n;b++) D(a,b)=dcmplx(0.0);
    for (const auto* o : os.template Iterate<Orbitals::TOrbital<U>>())
    {
        const double f=o->GetOccupation()*(energyWeighted ? o->GetEigenEnergy() : 1.0);
        if (f==0.0) continue;
        const vec_t<U>& c=o->GetCoeff();
        if (c.size()!=n) throw std::logic_error("SCFState: an orbital's coefficient count differs from its block's size");
        for (size_t a=0;a<n;a++)
            for (size_t b=0;b<n;b++)
            {
                if constexpr (std::is_same_v<U,double>) D(a,b)+=f*c[a]*c[b];
                else                                    D(a,b)+=f*c[a]*std::conj(c[b]);
            }
    }
    return D;
}
mat_t<dcmplx> BlockD(const StateBlock& r, bool energyWeighted=false)
{
    if (const auto* ro=Typed<double>(r); ro && r.real) return BlockD(*ro, r.n(), energyWeighted);
    if (const auto* co=Typed<dcmplx>(r); co && r.cplx) return BlockD(*co, r.n(), energyWeighted);
    throw std::logic_error("SCFState: a block's orbitals are not of its basis block's scalar");
}

//! \f$S^p\f$ of a Hermitian positive (semi)definite S, through its eigenpairs.  An eigenvalue below 1e-12 of the largest is a null
//! direction: it contributes NOTHING (the inverse root would amplify noise there; the ortho step drops those functions anyway).
mat_t<dcmplx> HermitianPower(const mat_t<dcmplx>& S, double p)
{
    const size_t n=S.rows();
    hmat_t<dcmplx> H(n);
    for (size_t a=0;a<n;a++) for (size_t b=a;b<n;b++) H(a,b)=0.5*(S(a,b)+std::conj(S(b,a)));
    rvec_t w; mat_t<dcmplx> U;
    blazem::eigen(H, w, U);
    const double floor=1e-12*w[n-1];
    mat_t<dcmplx> Uf(U);
    for (size_t k=0;k<n;k++)
    {
        const double f = w[k]>floor ? std::pow(w[k], p) : 0.0;
        for (size_t a=0;a<n;a++) Uf(a,k)*=f;
    }
    return Uf*blazem::ctrans(U);
}
} //anonymous namespace

//---------------------------------------------------------------------------------------------------
//  WRITE
//---------------------------------------------------------------------------------------------------
namespace
{
void WriteMaps(H5::Group g, const std::map<std::string,std::vector<double>>& num, const std::map<std::string,std::string>& str)
{
    for (const auto& [k,v] : num) g.Write(k, v);
    for (const auto& [k,v] : str) g.SetAttr(k, v);
}

// One block's D, C, ε, f.  C and D in C order; D = Σ_i f_i c_i c_i^† exactly as TOrbital::AddDensityMatrix.
template <class U> void WriteOrbitals(H5::Group& g, const Orbitals::TOrbitals<U>& os, size_t n)
{
    std::vector<const Orbitals::TOrbital<U>*> orbs;
    for (const auto* o : os.template Iterate<Orbitals::TOrbital<U>>()) orbs.push_back(o);
    const size_t nmo=orbs.size();
    std::vector<U> C(n*nmo), D(n*n, U(0));
    std::vector<double> eps(nmo), f(nmo);
    for (size_t i=0;i<nmo;i++)
    {
        const vec_t<U>& c=orbs[i]->GetCoeff();
        if (c.size()!=n) throw std::logic_error("WriteState: an orbital's coefficient count differs from its block's size");
        eps[i]=orbs[i]->GetEigenEnergy();
        f[i]  =orbs[i]->GetOccupation();
        for (size_t a=0;a<n;a++) C[a*nmo+i]=c[a];
        if (f[i]==0.0) continue;
        for (size_t a=0;a<n;a++)
            for (size_t b=0;b<n;b++)
            {
                if constexpr (std::is_same_v<U,double>) D[a*n+b]+=f[i]*c[a]*c[b];
                else                                    D[a*n+b]+=f[i]*c[a]*std::conj(c[b]);
            }
    }
    g.SetAttr("n",   std::int64_t(n));
    g.SetAttr("nmo", std::int64_t(nmo));
    g.Write("D",   D, {n,n});
    g.Write("C",   C, {n,nmo});
    g.Write("eps", eps);
    g.Write("f",   f);
}

// The block's overlap, in the block's own scalar (a real block's stays real on disk).
void WriteOverlap(H5::Group& g, const StateBlock& r)
{
    const size_t n=r.n();
    if (r.real)
    {
        const auto& s=r.real->Overlap();
        std::vector<double> S(n*n);
        for (size_t a=0;a<n;a++) for (size_t b=0;b<n;b++) S[a*n+b]=s(a,b);
        g.Write("S", S, {n,n});
    }
    else
    {
        const auto& s=r.cplx->Overlap();
        std::vector<dcmplx> S(n*n);
        for (size_t a=0;a<n;a++) for (size_t b=0;b<n;b++) S[a*n+b]=s(a,b);
        g.Write("S", S, {n,n});
    }
}

std::string UTCNow()
{
    const std::time_t t=std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&t));
    return buf;
}
} //anonymous namespace

void WriteState(const std::string& path, const std::string& kind, const StateFingerprint& fp, const StateSummary& sum,
                const std::vector<StateBlock>& blocks)
{
    const std::string tmp=path+".tmp";
    {
        H5::File f=H5::File::Create(tmp);
        f.SetAttr("format",     kFormat);
        f.SetAttr("kind",       kind);
        f.SetAttr("created",    UTCNow());
        f.SetAttr("label",      sum.label);
        f.SetAttr("converged",  std::int64_t(sum.converged ? 1 : 0));
        f.SetAttr("iterations", std::int64_t(sum.iterations));
        f.SetAttr("energy",     sum.energy);
        f.SetAttr("charge",     sum.charge);
        f.SetAttr("commutator", sum.commutator);
        f.SetAttr("recipe",     sum.recipe);
        f.SetAttr("densityChange", sum.densityChange);
        {
            H5::Group g=f.CreateGroup("fingerprint");
            WriteMaps(g.CreateGroup("refuse"), fp.refuseNum, fp.refuseStr);
            WriteMaps(g.CreateGroup("warm"),   fp.warmNum,   fp.warmStr);
        }
        H5::Group bg=f.CreateGroup("blocks");
        for (size_t b=0;b<blocks.size();b++)
        {
            const StateBlock& r=blocks[b];
            H5::Group g=bg.CreateGroup(std::to_string(b));
            g.SetAttr("ms",     std::int64_t(int(r.irrep.ms)));
            g.SetAttr("irrep",  r.Label());
            g.SetAttr("weight", r.weight);
            g.SetAttr("real",   std::int64_t(r.real ? 1 : 0));
            for (const auto& [k,v] : r.extras) g.Write(k, v);
            if      (const auto* ro=Typed<double>(r); ro && r.real) WriteOrbitals(g, *ro, r.n());
            else if (const auto* co=Typed<dcmplx>(r); co && r.cplx) WriteOrbitals(g, *co, r.n());
            else throw std::logic_error("WriteState: block "+std::to_string(b)+"'s orbitals are not of its basis block's scalar");
            WriteOverlap(g, r);
        }
        f.Flush();
    }
    std::filesystem::rename(tmp, path);   // atomic on one filesystem: the previous state survives a crash above
}

//---------------------------------------------------------------------------------------------------
//  READ
//---------------------------------------------------------------------------------------------------
namespace
{
void ReadMaps(const H5::Group& g, std::map<std::string,std::vector<double>>& num, std::map<std::string,std::string>& str)
{
    for (const std::string& k : g.Children())  num[k]=g.ReadReal(k);
    for (const std::string& k : g.AttrNames()) str[k]=g.AttrString(k);
}
bool IsVersionOf(const std::string& fmt, const std::string& major) {return fmt==major || fmt.rfind(major+".",0)==0;}
} //anonymous namespace

Outcome<SavedState, RestartRefusal> ReadState(const std::string& path)
{
    using O=Outcome<SavedState, RestartRefusal>;
    auto opened=H5::File::Open(path);
    if (!opened) return O::Fail({RestartRefusal::Why::Unreadable, opened.Error()});
    H5::File f=opened.TakeValue();
    if (!f.HasAttr("format") || !f.AttrIsString("format"))
        return O::Fail({RestartRefusal::Why::Unreadable, "'"+path+"' is HDF5 but not a qchem state (no format attribute)"});
    // Any 2.x, and the solid's 1.x: the major version is the one that may break a reader, a minor version only adds OPTIONAL
    // fields (read with HasAttr first, below).
    const std::string fmt=f.AttrString("format");
    const bool solid1=IsVersionOf(fmt, kSolidFormatMajor);
    if (!solid1 && !IsVersionOf(fmt, kFormatMajor))
        return O::Fail({RestartRefusal::Why::Format, "'"+path+"' is format '"+fmt+"'; this build reads '"+kFormatMajor+".x' and '"
                        +kSolidFormatMajor+".x'"});

    // Past the format check the file CLAIMS to be ours, so a missing piece is a broken file, not a legitimate "no state" -- the H5
    // layer's throws are allowed through.
    SavedState st;
    st.path=path;
    st.kind=solid1 ? "solid" : f.AttrString("kind");
    st.summary.label      =f.AttrString("label");
    st.summary.converged  =f.AttrInt("converged")!=0;
    st.summary.iterations =size_t(f.AttrInt("iterations"));
    st.summary.energy     =f.AttrReal("energy");
    st.summary.charge     =f.AttrReal("charge");
    st.summary.commutator =f.AttrReal("commutator");
    st.summary.recipe     =f.AttrString("recipe");
    if (f.HasAttr("densityChange")) st.summary.densityChange=f.AttrReal("densityChange");   // solid 1.1 on
    {
        H5::Group g=f.OpenGroup("fingerprint");
        ReadMaps(g.OpenGroup("refuse"), st.fingerprint.refuseNum, st.fingerprint.refuseStr);
        ReadMaps(g.OpenGroup("warm"),   st.fingerprint.warmNum,   st.fingerprint.warmStr);
    }
    static const std::set<std::string> payload={"D","C","eps","f","S"};
    H5::Group blocks=f.OpenGroup("blocks");
    for (size_t b=0; blocks.Has(std::to_string(b)); ++b)
    {
        H5::Group g=blocks.OpenGroup(std::to_string(b));
        SavedBlock sb;
        sb.ms    =g.AttrInt("ms");
        sb.weight=g.AttrReal("weight");
        sb.real  =g.AttrInt("real")!=0;
        if (g.HasAttr("irrep")) sb.irrep=g.AttrString("irrep");
        for (const std::string& k : g.Children())
            if (!payload.count(k)) sb.extras[k]=g.ReadReal(k);    // a solid's k (1.x and 2.x alike)
        const size_t n=size_t(g.AttrInt("n"));
        auto square=[&](const char* name)
        {
            if (g.Shape(name)!=std::vector<size_t>{n,n})
                throw std::runtime_error("ReadState ("+path+"): block "+std::to_string(b)+" "+name+" is not n x n");
            const std::vector<dcmplx> v=g.ReadComplex(name);
            mat_t<dcmplx> M(n,n);
            for (size_t a=0;a<n;a++) for (size_t c=0;c<n;c++) M(a,c)=v[a*n+c];
            return M;
        };
        sb.D=square("D");
        if (g.Has("S")) sb.S=square("S");
        st.blocks.push_back(std::move(sb));
    }
    if (st.fingerprint.refuseNum.count("blocks.n") && st.fingerprint.refuseNum["blocks.n"].size()!=st.blocks.size())
        throw std::runtime_error("ReadState ("+path+"): the fingerprint lists "
                                 +std::to_string(st.fingerprint.refuseNum["blocks.n"].size())+" blocks, the file holds "
                                 +std::to_string(st.blocks.size()));
    return O::Ok(std::move(st));
}

//---------------------------------------------------------------------------------------------------
//  SNAPSHOT, CARRY, RESTORE CHECKS
//---------------------------------------------------------------------------------------------------
SavedState SnapshotState(const std::vector<StateBlock>& blocks, const std::string& kind, bool energyWeighted)
{
    SavedState st;
    st.path="(in memory)";
    st.kind=kind;
    for (const StateBlock& r : blocks)
    {
        SavedBlock sb;
        sb.ms    =int(r.irrep.ms);
        sb.irrep =r.Label();
        sb.weight=r.weight;
        sb.real  =r.real!=nullptr;
        sb.extras=r.extras;
        sb.D     =BlockD(r, energyWeighted);
        sb.S     =BlockOverlap(r);
        st.blocks.push_back(std::move(sb));
    }
    return st;
}

void CarryToMovedAtoms(SavedState& carry, const std::vector<StateBlock>& blocks)
{
    if (blocks.size()!=carry.blocks.size())
        throw std::logic_error("CarryToMovedAtoms: the carry has "+std::to_string(carry.blocks.size())+" blocks, the basis "+std::to_string(blocks.size()));
    for (size_t b=0;b<blocks.size();b++)
    {
        SavedBlock& sb=carry.blocks[b];
        if (sb.D.rows()!=blocks[b].n())
            throw std::logic_error("CarryToMovedAtoms: block "+std::to_string(b)+" has "+std::to_string(sb.D.rows())+" functions, the basis "
                                   +std::to_string(blocks[b].n()));
        sb.S=BlockOverlap(blocks[b]);   // D stays: its AO coefficients follow their atoms
    }
}

void LowdinCarry(SavedState& carry, const std::vector<StateBlock>& blocks)
{
    if (blocks.size()!=carry.blocks.size())
        throw std::logic_error("LowdinCarry: the carry has "+std::to_string(carry.blocks.size())+" blocks, the basis "+std::to_string(blocks.size()));
    for (size_t b=0;b<blocks.size();b++)
    {
        SavedBlock& sb=carry.blocks[b];
        const size_t n=blocks[b].n();
        if (sb.S.rows()!=n || sb.D.rows()!=n)
            throw std::logic_error("LowdinCarry: block "+std::to_string(b)+" has no saved overlap of its size (a solid 1.x file has none)");
        const mat_t<dcmplx> Snew=BlockOverlap(blocks[b]);
        const mat_t<dcmplx> A=HermitianPower(Snew, -0.5)*HermitianPower(sb.S, 0.5);
        sb.D=A*sb.D*blazem::ctrans(A);
        sb.S=Snew;   // the carry now lives on the new basis
    }
}

double Electrons(const SavedState& st, const std::vector<StateBlock>& blocks)
{
    if (blocks.size()!=st.blocks.size()) throw std::logic_error("SCFState::Electrons: the block counts differ");
    double N=0.0;
    for (size_t b=0;b<blocks.size();b++)
    {
        const mat_t<dcmplx> S=BlockOverlap(blocks[b]);
        const mat_t<dcmplx>& D=st.blocks[b].D;
        if (D.rows()!=S.rows()) throw std::logic_error("SCFState::Electrons: block "+std::to_string(b)+"'s sizes differ");
        double tr=0.0;
        for (size_t a=0;a<D.rows();a++) for (size_t c=0;c<D.rows();c++) tr+=(D(a,c)*S(c,a)).real();
        N+=blocks[b].weight*tr;
    }
    return N;
}

Outcome<std::vector<mat_t<dcmplx>>, RestartRefusal> CheckedBlockDensities(const SavedState& st, const std::vector<StateBlock>& blocks)
{
    using O=Outcome<std::vector<mat_t<dcmplx>>, RestartRefusal>;
    // The fingerprint comparison has already judged these; a disagreement HERE means the fingerprint is missing an entry that
    // matters -- a defect in the code, not the caller's input.
    if (blocks.size()!=st.blocks.size())
        throw std::logic_error("SCFState restore: "+std::to_string(st.blocks.size())+" saved blocks for "+std::to_string(blocks.size())
                               +" in this run -- the fingerprint should have refused");
    std::vector<mat_t<dcmplx>> out;
    for (size_t b=0;b<blocks.size();b++)
    {
        const StateBlock& r=blocks[b];
        const SavedBlock& sb=st.blocks[b];
        const size_t n=r.n();
        if (sb.D.rows()!=n || sb.ms!=int(r.irrep.ms))
            throw std::logic_error("SCFState restore: block "+std::to_string(b)+" does not match this run's -- the fingerprint should have refused");
        if (r.real)
        {
            // A complex-saved TRIM block restarting REAL: its D is real up to eigenvector-phase noise, or the two runs are not the
            // same state.
            double imax=0.0;
            for (size_t a=0;a<n;a++) for (size_t c=0;c<n;c++) imax=std::max(imax, std::fabs(sb.D(a,c).imag()));
            if (imax>1e-8)
                return O::Fail({RestartRefusal::Why::Realness, "block "+std::to_string(b)+" restarts REAL but its saved D has |Im D| up to "
                                +std::to_string(imax)+" (run with forceComplex, or resave)"});
        }
        mat_t<dcmplx> D(n,n);
        for (size_t a=0;a<n;a++)
        {
            D(a,a)=dcmplx(sb.D(a,a).real(), 0.0);
            for (size_t c=a+1;c<n;c++) { D(a,c)=0.5*(sb.D(a,c)+std::conj(sb.D(c,a))); D(c,a)=std::conj(D(a,c)); }
        }
        if (r.real) for (size_t a=0;a<n;a++) for (size_t c=0;c<n;c++) D(a,c)=dcmplx(D(a,c).real(), 0.0);
        out.push_back(std::move(D));
    }
    return O::Ok(std::move(out));
}

} //namespace qchem
