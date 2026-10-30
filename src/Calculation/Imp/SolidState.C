// File: Calculation/Imp/SolidState.C  The saved periodic SCF state: fingerprint, writer, reader, restore.
module;
#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdint>
#include <ctime>
#include <filesystem>
#include <iomanip>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.SolidState;

import qchem.HDF5;
import qchem.Symmetry.Irrep;                    // Irrep, Spin, SpinIrreps
import qchem.Symmetry.Lattice_3D.BlochQN;       // Getk (the block's crystal momentum)
import qchem.BasisSet.Orbital_1E_IBS;           // the per-block bases (real and complex)
import qchem.BasisSet.AoShellSource;            // the orbital shells (the basis half of the fingerprint)
import qchem.Orbitals;                          // TOrbitals/TOrbital (C, ε, f)
import qchem.CompositeCD;                       // tComposite_CD (the restored density's shape)
import qchem.ChargeDensity.Factory;             // IrrepCD_Factory

namespace qchem
{

static const char* kFormat = "qchem-solid-state 1.1";   //!< 1.1 (2026-10-30): + densityChange
static const char* kFormatMajor = "qchem-solid-state 1"; //!< every 1.x is read: a minor version only ADDS optional fields

//---------------------------------------------------------------------------------------------------
//  The canonical BLOCK WALK -- spin irreps outer, basis blocks inner: the order tCompositeWF builds its
//  children in.  Every writer, reader and restorer walks it, so "block b" means one thing in this file.
//---------------------------------------------------------------------------------------------------
namespace
{
struct BlockRef
{
    Spin                                     s;
    size_t                                   i;
    const BasisSet::Orbital_1E_IBS<double>*  real;   //!< set on a REAL (TRIM) block, else null
    const BasisSet::Orbital_1E_IBS<dcmplx>*  cplx;   //!< set on a complex block, else null
    Irrep                                    irrep;
    size_t                                   n;
};
std::vector<BlockRef> Blocks(const BasisSet::Complex_BS& bs, SpinGroup g)
{
    std::vector<BlockRef> out;
    for (Spin s : SpinIrreps(g))
        for (size_t i=0;i<bs.GetNumIBS();++i)
        {
            if (const auto* rb=bs.GetRealIBS(i))
                out.push_back({s, i, rb, nullptr, rb->GetIrrep(s), rb->GetNumFunctions()});
            else
            {
                const auto* cb=bs[i];
                out.push_back({s, i, nullptr, cb, cb->GetIrrep(s), cb->GetNumFunctions()});
            }
        }
    return out;
}

std::string SpinGroupName(SpinGroup g) {return g==SpinGroup::Polarized ? "Polarized" : "UnPolarized";}

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
        {"gpw.",                  "a GPW numerical tolerance changed"+moves},
        {"ortho.",                "the basis orthogonalisation changed: functions may be dropped or kept differently"+moves},
        {"scf.kT",                "the smearing changed: the occupations are re-filled at the new temperature"+moves},
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
} //anonymous namespace

//---------------------------------------------------------------------------------------------------
StateFingerprint MakeStateFingerprint(const rmat3d_t& M, const Structure& st, const BasisSet::Complex_BS& bs,
                                      SpinGroup g, const RunIdentity& id)
{
    StateFingerprint fp;
    // ---- REFUSE: the saved D is a matrix over these functions, at these k, for this many electrons ----
    {
        std::vector<double> A;                                   // rows = a_1, a_2, a_3 (Matrix3D is 1-based)
        for (size_t col=1;col<=3;col++) for (size_t row=1;row<=3;row++) A.push_back(M(row,col));
        fp.refuseNum["structure.cell"]=A;
        std::vector<double> Z, R;
        st.ForEachSite([&](int z, const rvec3_t& r, bool){ Z.push_back(z); R.push_back(r.x); R.push_back(r.y); R.push_back(r.z); });
        fp.refuseNum["structure.Z"]=Z;
        fp.warmNum["structure.positions"]=R;                     // WARM: a geometry step carries D (the basis centres follow)
    }
    {
        // Canonical (sorted) so listing the same species in another order is not a "different run".
        std::vector<std::pair<std::string,int>> sp=id.species;
        std::sort(sp.begin(), sp.end());
        std::string s;
        for (const auto& [el,val] : sp) s += (s.empty() ? "" : ",") + el + ":" + std::to_string(val);
        fp.refuseStr["species"]=s;
        fp.refuseStr["pp.set"]=id.ppSet;
    }
    fp.refuseStr["basis.family"]=id.orbitalBasis;
    if (id.orbitalBasis=="planeWave")
        // A PLANE-WAVE basis has no shells: its functions are the G-vectors with |k+G|^2/2 < Ecut, so the cutoff -- with the cell
        // and each block's k and size, all REFUSE -- names them.  (OpenWork §3 "PW states cannot be saved", 2026-10-30.)
        fp.refuseNum["basis.planeWaveEcut"]={id.planeWaveEcut};
    else
    {
        // The ORBITAL basis, shell by shell in block order: every Bloch block is the Bloch sum of these shells
        // (AoShellSource's own contract), so block 0 speaks for all of them.
        const BasisSet::AoShellSource* src = bs.GetRealIBS(0)
            ? dynamic_cast<const BasisSet::AoShellSource*>(bs.GetRealIBS(0))
            : dynamic_cast<const BasisSet::AoShellSource*>(bs[0]);
        if (!src) throw std::logic_error("MakeStateFingerprint: the Bloch block does not report its AO shells "
                                         "(AoShellSource) -- a saved state cannot identify this basis");
        std::vector<double> L, nc, np, cen, ex, co;
        for (const auto& sh : src->GetAoShells())
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
    {
        std::vector<double> ms, k, w, n, real;
        for (const BlockRef& b : Blocks(bs, g))
        {
            ms.push_back(double(int(b.irrep.ms)));
            const rvec3_t kk=Symmetry::Lattice_3D::Getk(*b.irrep.sym);
            k.push_back(kk.x); k.push_back(kk.y); k.push_back(kk.z);
            w.push_back(b.irrep.sym->GetWeight());
            n.push_back(double(b.n));
            real.push_back(b.real ? 1.0 : 0.0);
        }
        fp.refuseNum["blocks.ms"]=ms;
        fp.refuseNum["blocks.k"]=k;
        fp.refuseNum["blocks.weight"]=w;
        fp.refuseNum["blocks.n"]=n;
        fp.warmNum  ["blocks.real"]=real;                         // a working-type change promotes/demotes D
    }
    fp.refuseNum["electrons.Nelec"]={double(id.Nelec)};
    fp.refuseStr["spinGroup"]=SpinGroupName(g);
    fp.warmStr["functional"]=id.functional;

    // ---- WARM: a restart across these is what the mechanism is FOR ----
    {
        std::vector<double> site, l, U, alpha, alphaUp, alphaDown, projector, ortho, nUirrep, Uirrep;
        for (const auto& M : id.hubbard)
        {
            site.push_back(M.site); l.push_back(M.l); U.push_back(M.U); alpha.push_back(M.alpha); alphaUp.push_back(M.alphaUp); alphaDown.push_back(M.alphaDown);
            projector.push_back(M.projector==Hamiltonian::HubbardProjectorKind::PseudoAtomic ? 1 : 0); ortho.push_back(M.orthoAtomic ? 1 : 0);
            nUirrep.push_back(M.Uirrep.size());
            for (double u : M.Uirrep) Uirrep.push_back(u);
        }
        fp.warmNum["hubbard.site"]=site;   fp.warmNum["hubbard.l"]=l;
        fp.warmNum["hubbard.U"]=U;         fp.warmNum["hubbard.alpha"]=alpha;
        bool spinAlpha=false; for (const auto& M : id.hubbard) spinAlpha|=(M.alphaUp!=0.0 || M.alphaDown!=0.0);
        if (spinAlpha) { fp.warmNum["hubbard.alphaUp"]=alphaUp; fp.warmNum["hubbard.alphaDown"]=alphaDown; }   // only when set: older states stay EXACT
        std::vector<double> nAlphaSlot, alphaSlot; bool slotAlpha=false;
        for (const auto& M : id.hubbard) { nAlphaSlot.push_back(M.alphaSlot.size()); for (double a : M.alphaSlot) alphaSlot.push_back(a); slotAlpha|=!M.alphaSlot.empty(); }
        if (slotAlpha) { fp.warmNum["hubbard.nAlphaSlot"]=nAlphaSlot; fp.warmNum["hubbard.alphaSlot"]=alphaSlot; }        // likewise
        fp.warmNum["hubbard.projector"]=projector; fp.warmNum["hubbard.orthoAtomic"]=ortho;   // 1 = pseudo-atomic (was hubbard.atomicRadial until 2026-10-25)
        fp.warmNum["hubbard.nUirrep"]=nUirrep;     fp.warmNum["hubbard.Uirrep"]=Uirrep;
    }
    fp.warmNum["spin.multiplicity"]={double(id.multiplicity)};
    fp.warmNum["grids.densityEcut"]={id.densityEcut};
    fp.warmNum["grids.cutoffFactor"]={id.cutoffFactor};
    fp.warmStr["grids.xcMesh"]=id.xcMesh;
    fp.warmNum["ortho.kind"]={double(id.ortho)};
    fp.warmNum["ortho.tol"]={id.orthoTol};
    fp.warmNum["scf.kT"]={id.kT};
    // The numerics that move the ground state (gap 2, 2026-10-30).  Each key exists only where it means something, so a
    // plane-wave state carries no GPW key and a uniform-mesh state no Becke eps.
    if (!id.vxcFit.empty())     fp.warmStr["xc.fit"]=id.vxcFit;
    if (id.beckeEps>0.0)        fp.warmNum["grids.beckeEps"]={id.beckeEps};
    if (!id.raster.empty())     fp.warmStr["grids.raster"]=id.raster;
    if (id.ladderFactor>0.0)    fp.warmNum["grids.ladderFactor"]={id.ladderFactor};
    for (const auto& [k,v] : id.gpwTolerances) fp.warmNum["gpw."+k]=v;
    return fp;
}

//! A state written before C4 step 5 filed the functional, the atom positions and the shell centres under REFUSE and had no
//! pp.set (nor, before 2026-10-30, basis.family): move them to the class they belong to now, and take the run's pp.set for a state that never recorded one (its
//! species string, still REFUSE, carries the q per element).  Additive: an old file keeps loading.
static StateFingerprint Migrated(StateFingerprint s, const StateFingerprint& now)
{
    for (const char* k : {"structure.positions", "basis.shellCenters"})
        if (auto it=s.refuseNum.find(k); it!=s.refuseNum.end()) { s.warmNum[k]=std::move(it->second); s.refuseNum.erase(it); }
    if (auto it=s.refuseStr.find("functional"); it!=s.refuseStr.end()) { s.warmStr["functional"]=std::move(it->second); s.refuseStr.erase(it); }
    if (!s.refuseStr.count("pp.set") && now.refuseStr.count("pp.set")) s.refuseStr["pp.set"]=now.refuseStr.at("pp.set");
    // basis.family (2026-10-30): every state written before it was Gaussian -- a plane-wave state could not be written at all.
    if (!s.refuseStr.count("basis.family")) s.refuseStr["basis.family"]="gaussian";
    // The numerics keys (2026-10-30): a state written before them recorded none, so none of them is compared -- the run's value
    // stands in.  Only when the state lacks the WHOLE set: a newer state missing one of them really was made without it.
    const bool predatesNumerics=!s.warmStr.count("xc.fit");
    if (predatesNumerics)
    {
        for (const auto& [k,v] : now.warmStr) if (k=="xc.fit" || k=="grids.raster") s.warmStr[k]=v;
        for (const auto& [k,v] : now.warmNum) if (k=="grids.beckeEps" || k=="grids.ladderFactor" || k.rfind("gpw.",0)==0) s.warmNum[k]=v;
    }
    return s;
}

//! Largest atom displacement between two flattened position lists, or -1 if they are not comparable.
static double MaxDisplacement(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size()!=b.size() || a.size()%3!=0) return -1.0;
    double m=0.0;
    for (size_t i=0;i<a.size();i+=3)
        m=std::max(m, std::sqrt((a[i]-b[i])*(a[i]-b[i]) + (a[i+1]-b[i+1])*(a[i+1]-b[i+1]) + (a[i+2]-b[i+2])*(a[i+2]-b[i+2])));
    return m;
}

Outcome<std::vector<std::string>, RestartRefusal> CompareFingerprints(const StateFingerprint& savedAsWritten,
                                                                     const StateFingerprint& now)
{
    using O=Outcome<std::vector<std::string>, RestartRefusal>;
    const StateFingerprint saved=Migrated(savedAsWritten, now);
    const std::vector<std::string> refuse=Differences(saved.refuseNum, now.refuseNum, saved.refuseStr, now.refuseStr);
    if (!refuse.empty())
    {
        std::string d="the saved state was made for a DIFFERENT run:";
        for (const std::string& r : refuse) d+="\n    "+r;
        return O::Fail({RestartRefusal::Why::Mismatch, d});
    }
    std::vector<std::string> warm=Differences(saved.warmNum, now.warmNum, saved.warmStr, now.warmStr);
    // A moved atom is ONE fact: say how far, and drop the basis-centre line that is only its consequence.
    const auto sp=saved.warmNum.find("structure.positions"), np=now.warmNum.find("structure.positions");
    if (sp!=saved.warmNum.end() && np!=now.warmNum.end())
        if (const double d=MaxDisplacement(sp->second, np->second); d>1e-10)
        {
            std::erase_if(warm, [](const std::string& l){ return l.rfind("structure.positions:",0)==0 || l.rfind("basis.shellCenters:",0)==0; });
            std::ostringstream os;
            os<<"structure.positions: atoms moved, max displacement "<<std::setprecision(6)<<d<<" bohr (the basis follows)";
            warm.push_back(os.str());
        }
    return O::Ok(std::move(warm));
}

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

// One block's physical D = Σ_i f_i c_i c_i^† exactly as TOrbital::AddDensityMatrix, held complex (a real block's promotes exactly).
template <class U> mat_t<dcmplx> BlockD(const Orbitals::TOrbitals<U>& os, size_t n)
{
    mat_t<dcmplx> D(n,n);
    for (size_t a=0;a<n;a++) for (size_t b=0;b<n;b++) D(a,b)=dcmplx(0.0);
    for (const auto* o : os.template Iterate<Orbitals::TOrbital<U>>())
    {
        const double f=o->GetOccupation();
        if (f==0.0) continue;
        const vec_t<U>& c=o->GetCoeff();
        if (c.size()!=n) throw std::logic_error("SnapshotState: an orbital's coefficient count differs from its block's size");
        for (size_t a=0;a<n;a++)
            for (size_t b=0;b<n;b++)
            {
                if constexpr (std::is_same_v<U,double>) D(a,b)+=f*c[a]*c[b];
                else                                    D(a,b)+=f*c[a]*std::conj(c[b]);
            }
    }
    return D;
}

// One block's D, C, ε, f.  C and D in C order; D = Σ_i f_i c_i c_i^† exactly as TOrbital::AddDensityMatrix.
template <class U> void WriteBlock(H5::Group& g, const Orbitals::TOrbitals<U>& os, size_t n)
{
    std::vector<const Orbitals::TOrbital<U>*> orbs;
    for (const auto* o : os.template Iterate<Orbitals::TOrbital<U>>()) orbs.push_back(o);
    const size_t nmo=orbs.size();
    std::vector<U> C(n*nmo), D(n*n, U(0));
    std::vector<double> eps(nmo), f(nmo);
    for (size_t i=0;i<nmo;i++)
    {
        const vec_t<U>& c=orbs[i]->GetCoeff();
        if (c.size()!=n) throw std::logic_error("WriteSolidState: an orbital's coefficient count differs from its block's size");
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

std::string UTCNow()
{
    const std::time_t t=std::chrono::system_clock::to_time_t(std::chrono::system_clock::now());
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%dT%H:%M:%SZ", std::gmtime(&t));
    return buf;
}
} //anonymous namespace

void WriteSolidState(const std::string& path, const StateFingerprint& fp, const StateSummary& sum,
                     const BasisSet::Complex_BS& bs, const WaveFunction::cWaveFunction& wf)
{
    const std::string tmp=path+".tmp";
    {
        H5::File f=H5::File::Create(tmp);
        f.SetAttr("format",     kFormat);
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
        H5::Group blocks=f.CreateGroup("blocks");
        const std::vector<BlockRef> refs=Blocks(bs, wf.GetSpinGroup());
        for (size_t b=0;b<refs.size();b++)
        {
            const BlockRef& r=refs[b];
            H5::Group g=blocks.CreateGroup(std::to_string(b));
            g.SetAttr("ms",     std::int64_t(int(r.irrep.ms)));
            g.SetAttr("weight", r.irrep.sym->GetWeight());
            g.SetAttr("real",   std::int64_t(r.real ? 1 : 0));
            const rvec3_t k=Symmetry::Lattice_3D::Getk(*r.irrep.sym);
            g.Write("k", std::vector<double>{k.x, k.y, k.z});
            const Orbitals::Orbitals* os=wf.GetOrbitals(r.irrep);
            if      (const auto* ro=dynamic_cast<const Orbitals::TOrbitals<double>*>(os); ro && r.real) WriteBlock(g, *ro, r.n);
            else if (const auto* co=dynamic_cast<const Orbitals::TOrbitals<dcmplx>*>(os); co && r.cplx) WriteBlock(g, *co, r.n);
            else throw std::logic_error("WriteSolidState: block "+std::to_string(b)+"'s orbitals are not of its basis block's scalar");
        }
        f.Flush();
    }
    std::filesystem::rename(tmp, path);   // atomic on one filesystem: the previous state survives a crash above
}

SavedState SnapshotState(const BasisSet::Complex_BS& bs, const WaveFunction::cWaveFunction& wf)
{
    SavedState st;
    st.path="(in memory)";
    for (const BlockRef& r : Blocks(bs, wf.GetSpinGroup()))
    {
        SavedBlock sb;
        sb.ms    =int(r.irrep.ms);
        sb.k     =Symmetry::Lattice_3D::Getk(*r.irrep.sym);
        sb.weight=r.irrep.sym->GetWeight();
        sb.real  =r.real!=nullptr;
        const Orbitals::Orbitals* os=wf.GetOrbitals(r.irrep);
        if      (const auto* ro=dynamic_cast<const Orbitals::TOrbitals<double>*>(os); ro && r.real) sb.D=BlockD(*ro, r.n);
        else if (const auto* co=dynamic_cast<const Orbitals::TOrbitals<dcmplx>*>(os); co && r.cplx) sb.D=BlockD(*co, r.n);
        else throw std::logic_error("SnapshotState: a block's orbitals are not of its basis block's scalar");
        st.blocks.push_back(std::move(sb));
    }
    return st;
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
} //anonymous namespace

Outcome<SavedState, RestartRefusal> ReadSolidState(const std::string& path)
{
    using O=Outcome<SavedState, RestartRefusal>;
    auto opened=H5::File::Open(path);
    if (!opened) return O::Fail({RestartRefusal::Why::Unreadable, opened.Error()});
    H5::File f=opened.TakeValue();
    if (!f.HasAttr("format") || !f.AttrIsString("format"))
        return O::Fail({RestartRefusal::Why::Unreadable, "'"+path+"' is HDF5 but not a qchem solid state (no format attribute)"});
    // Any 1.x: the major version is the one that may break a reader, a minor version only adds OPTIONAL fields (read with
    // HasAttr first, below).  "qchem-solid-state 1" itself is 1.0.
    if (const std::string fmt=f.AttrString("format"); fmt!=kFormatMajor && fmt.rfind(std::string(kFormatMajor)+".",0)!=0)
        return O::Fail({RestartRefusal::Why::Format, "'"+path+"' is format '"+fmt+"'; this build reads '"+kFormatMajor+".x'"});

    // Past the format check the file CLAIMS to be ours, so a missing piece is a broken file, not a
    // legitimate "no state" -- the H5 layer's throws are allowed through.
    SavedState st;
    st.path=path;
    st.summary.label      =f.AttrString("label");
    st.summary.converged  =f.AttrInt("converged")!=0;
    st.summary.iterations =size_t(f.AttrInt("iterations"));
    st.summary.energy     =f.AttrReal("energy");
    st.summary.charge     =f.AttrReal("charge");
    st.summary.commutator =f.AttrReal("commutator");
    st.summary.recipe     =f.AttrString("recipe");
    // OPTIONAL fields: absent = a state written before they existed (the default stands, and says so where it is used).
    if (f.HasAttr("densityChange")) st.summary.densityChange=f.AttrReal("densityChange");   // 1.1
    {
        H5::Group g=f.OpenGroup("fingerprint");
        ReadMaps(g.OpenGroup("refuse"), st.fingerprint.refuseNum, st.fingerprint.refuseStr);
        ReadMaps(g.OpenGroup("warm"),   st.fingerprint.warmNum,   st.fingerprint.warmStr);
    }
    H5::Group blocks=f.OpenGroup("blocks");
    for (size_t b=0; blocks.Has(std::to_string(b)); ++b)
    {
        H5::Group g=blocks.OpenGroup(std::to_string(b));
        SavedBlock sb;
        sb.ms    =g.AttrInt("ms");
        sb.weight=g.AttrReal("weight");
        sb.real  =g.AttrInt("real")!=0;
        const std::vector<double> k=g.ReadReal("k");
        if (k.size()!=3) throw std::runtime_error("ReadSolidState ("+path+"): block "+std::to_string(b)+" k is not a 3-vector");
        sb.k=rvec3_t(k[0], k[1], k[2]);
        const size_t n=size_t(g.AttrInt("n"));
        if (g.Shape("D")!=std::vector<size_t>{n,n})
            throw std::runtime_error("ReadSolidState ("+path+"): block "+std::to_string(b)+" D is not n x n");
        const std::vector<dcmplx> D=g.ReadComplex("D");
        sb.D.resize(n,n);
        for (size_t a=0;a<n;a++) for (size_t c=0;c<n;c++) sb.D(a,c)=D[a*n+c];
        st.blocks.push_back(std::move(sb));
    }
    if (st.fingerprint.refuseNum.count("blocks.n") && st.fingerprint.refuseNum["blocks.n"].size()!=st.blocks.size())
        throw std::runtime_error("ReadSolidState ("+path+"): the fingerprint lists "
                                 +std::to_string(st.fingerprint.refuseNum["blocks.n"].size())+" blocks, the file holds "
                                 +std::to_string(st.blocks.size()));
    return O::Ok(std::move(st));
}

//---------------------------------------------------------------------------------------------------
//  RESTORE
//---------------------------------------------------------------------------------------------------
Outcome<std::unique_ptr<ChargeDensity::cDM_CD>, RestartRefusal>
RestoreDensity(const SavedState& st, const BasisSet::Complex_BS& bs, SpinGroup g)
{
    using O=Outcome<std::unique_ptr<ChargeDensity::cDM_CD>, RestartRefusal>;
    using namespace ChargeDensity;
    const std::vector<BlockRef> refs=Blocks(bs, g);
    // The fingerprint comparison has already judged these; a disagreement HERE means the fingerprint is
    // missing an entry that matters -- a defect in this file, not the caller's input.
    if (refs.size()!=st.blocks.size())
        throw std::logic_error("RestoreDensity: "+std::to_string(st.blocks.size())+" saved blocks for "
                               +std::to_string(refs.size())+" in this run -- the fingerprint should have refused");
    // Pass 1: every block's weighted D, as saved.
    std::vector<hmat_t<double>> Dr(refs.size());
    std::vector<hmat_t<dcmplx>> Dc(refs.size());
    for (size_t b=0;b<refs.size();b++)
    {
        const BlockRef& r=refs[b];
        const SavedBlock& sb=st.blocks[b];
        const rvec3_t k=Symmetry::Lattice_3D::Getk(*r.irrep.sym);
        if (sb.D.rows()!=r.n || sb.ms!=int(r.irrep.ms) || norm(sb.k-k)>1e-10)
            throw std::logic_error("RestoreDensity: block "+std::to_string(b)+" does not match this run's -- the "
                                   "fingerprint should have refused");
        const double w=r.irrep.sym->GetWeight();
        const size_t n=r.n;
        if (r.real)
        {
            // A complex-saved TRIM block restarting REAL: its D is real up to eigenvector-phase noise, or
            // the two runs are not the same state.
            double imax=0.0;
            for (size_t a=0;a<n;a++) for (size_t c=0;c<n;c++) imax=std::max(imax, std::fabs(sb.D(a,c).imag()));
            if (imax>1e-8)
                return O::Fail({RestartRefusal::Why::Realness, "block "+std::to_string(b)+" restarts REAL but its saved D "
                                "has |Im D| up to "+std::to_string(imax)+" (run with forceComplex, or resave)"});
            hmat_t<double> Dw(n);
            for (size_t a=0;a<n;a++) for (size_t c=a;c<n;c++) Dw(a,c)=0.5*w*(sb.D(a,c).real()+sb.D(c,a).real());
            Dr[b]=std::move(Dw);
        }
        else
        {
            hmat_t<dcmplx> Dw(n);
            for (size_t a=0;a<n;a++)
            {
                Dw(a,a)=dcmplx(w*sb.D(a,a).real(), 0.0);           // a HermitianMatrix diagonal must be real
                for (size_t c=a+1;c<n;c++) Dw(a,c)=0.5*w*(sb.D(a,c)+std::conj(sb.D(c,a)));
            }
            Dc[b]=std::move(Dw);
        }
    }
    // Pass 2: THE DENSITY THE RUN USED.  A polarized run's density is the Shubnikov PAIR projection of each block's
    // (D_up, D_dn) (tCompositeWF::GetChargeDensity, 3b2a163a), but the file holds the ORBITALS' D (Sigma f c c^dag).  At a
    // converged imposed state the two agree (beta-MnO2 A1AFM r015: E[D_saved] reproduced the saved E(internal) to 1e-10 with
    // and without this), so this is the contract, not a measured repair: the restored density is the projected one, as every
    // later density is.  Blocks are spin-major (Blocks), so block i's twin is i + N.  A no-op without a magnetic group.
    if (g==SpinGroup::Polarized)
    {
        const size_t N=refs.size()/2;
        assert(2*N==refs.size());
        for (size_t b=0;b<N;b++)
            if (refs[b].real) bs.ProjectSpinPair(Dr[b], Dr[b+N]);
            else              bs.ProjectSpinPair(Dc[b], Dc[b+N]);
    }
    auto cd=std::make_unique<tComposite_CD<dcmplx>>(bs.GetReciprocalPointOps());
    for (size_t b=0;b<refs.size();b++)
    {
        const BlockRef& r=refs[b];
        if (r.real) cd->Insert(std::unique_ptr<tDM_CD<double>>(IrrepCD_Factory<double>(Dr[b], r.real, r.irrep)), r.irrep);
        else        cd->Insert(std::unique_ptr<tDM_CD<dcmplx>>(IrrepCD_Factory<dcmplx>(Dc[b], r.cplx, r.irrep)), r.irrep);
    }
    return O::Ok(std::unique_ptr<cDM_CD>(std::move(cd)));
}

} //namespace qchem
