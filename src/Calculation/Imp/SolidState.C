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

static const char* kFormat = "qchem-solid-state 1";

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
bool Same(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size()!=b.size()) return false;
    for (size_t k=0;k<a.size();k++)
        if (std::fabs(a[k]-b[k]) > 1e-10*std::max(1.0, std::fabs(a[k]))) return false;
    return true;
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
        fp.refuseNum["structure.positions"]=R;
    }
    {
        // Canonical (sorted) so listing the same species in another order is not a "different run".
        std::vector<std::pair<std::string,int>> sp=id.species;
        std::sort(sp.begin(), sp.end());
        std::string s;
        for (const auto& [el,val] : sp) s += (s.empty() ? "" : ",") + el + ":" + std::to_string(val);
        fp.refuseStr["species"]=s;
    }
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
        fp.refuseNum["basis.shellCenters"]=cen;
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
    fp.refuseStr["functional"]=id.functional;

    // ---- WARM: a restart across these is what the mechanism is FOR ----
    {
        std::vector<double> site, l, U, alpha, atomic, ortho, nUirrep, Uirrep;
        for (const auto& M : id.hubbard)
        {
            site.push_back(M.site); l.push_back(M.l); U.push_back(M.U); alpha.push_back(M.alpha);
            atomic.push_back(M.atomicRadial ? 1 : 0); ortho.push_back(M.orthoAtomic ? 1 : 0);
            nUirrep.push_back(M.Uirrep.size());
            for (double u : M.Uirrep) Uirrep.push_back(u);
        }
        fp.warmNum["hubbard.site"]=site;   fp.warmNum["hubbard.l"]=l;
        fp.warmNum["hubbard.U"]=U;         fp.warmNum["hubbard.alpha"]=alpha;
        fp.warmNum["hubbard.atomicRadial"]=atomic; fp.warmNum["hubbard.orthoAtomic"]=ortho;
        fp.warmNum["hubbard.nUirrep"]=nUirrep;     fp.warmNum["hubbard.Uirrep"]=Uirrep;
    }
    fp.warmNum["spin.multiplicity"]={double(id.multiplicity)};
    fp.warmNum["grids.densityEcut"]={id.densityEcut};
    fp.warmNum["grids.cutoffFactor"]={id.cutoffFactor};
    fp.warmStr["grids.xcMesh"]=id.xcMesh;
    fp.warmNum["ortho.kind"]={double(id.ortho)};
    fp.warmNum["ortho.tol"]={id.orthoTol};
    fp.warmNum["scf.kT"]={id.kT};
    return fp;
}

Outcome<std::vector<std::string>, RestartRefusal> CompareFingerprints(const StateFingerprint& saved,
                                                                     const StateFingerprint& now)
{
    using O=Outcome<std::vector<std::string>, RestartRefusal>;
    const std::vector<std::string> refuse=Differences(saved.refuseNum, now.refuseNum, saved.refuseStr, now.refuseStr);
    if (!refuse.empty())
    {
        std::string d="the saved state was made for a DIFFERENT run:";
        for (const std::string& r : refuse) d+="\n    "+r;
        return O::Fail({RestartRefusal::Why::Mismatch, d});
    }
    return O::Ok(Differences(saved.warmNum, now.warmNum, saved.warmStr, now.warmStr));
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
    if (const std::string fmt=f.AttrString("format"); fmt!=kFormat)
        return O::Fail({RestartRefusal::Why::Format, "'"+path+"' is format '"+fmt+"'; this build reads '"+kFormat+"'"});

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
    auto cd=std::make_unique<tComposite_CD<dcmplx>>(bs.GetReciprocalPointOps());
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
            cd->Insert(std::unique_ptr<tDM_CD<double>>(IrrepCD_Factory<double>(Dw, r.real, r.irrep)), r.irrep);
        }
        else
        {
            hmat_t<dcmplx> Dw(n);
            for (size_t a=0;a<n;a++)
            {
                Dw(a,a)=dcmplx(w*sb.D(a,a).real(), 0.0);           // a HermitianMatrix diagonal must be real
                for (size_t c=a+1;c<n;c++) Dw(a,c)=0.5*w*(sb.D(a,c)+std::conj(sb.D(c,a)));
            }
            cd->Insert(std::unique_ptr<tDM_CD<dcmplx>>(IrrepCD_Factory<dcmplx>(Dw, r.cplx, r.irrep)), r.irrep);
        }
    }
    return O::Ok(std::unique_ptr<cDM_CD>(std::move(cd)));
}

} //namespace qchem
