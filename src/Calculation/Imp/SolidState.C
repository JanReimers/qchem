// File: Calculation/Imp/SolidState.C  The saved periodic SCF state: the solid's fingerprint, block walk and density restore.
// (The framing, the writer and reader, the comparison and the Löwdin carry are qchem.SCFState's since 2026-11-09.)
module;
#include <algorithm>
#include <cassert>
#include <cmath>
#include <complex>
#include <cstdint>
#include <iomanip>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.SolidState;

import qchem.BasisSet.Gaussian.Lattice.GPWTolerances;   // kRungGateC: what a state with no gpw.rungGateC was built on
import qchem.Symmetry.Irrep;                    // Irrep, Spin, SpinIrreps
import qchem.Symmetry.Lattice_3D.BlochQN;       // Getk (the block's crystal momentum)
import qchem.BasisSet.Orbital_1E_IBS;           // the per-block bases (real and complex)
import qchem.CompositeCD;                       // tComposite_CD (the restored density's shape)
import qchem.ChargeDensity.Factory;             // IrrepCD_Factory

namespace qchem
{

//---------------------------------------------------------------------------------------------------
//  The canonical BLOCK WALK -- spin irreps outer, basis blocks inner: the order tCompositeWF builds its
//  children in.  Every writer, reader and restorer walks it, so "block b" means one thing in a state.
//---------------------------------------------------------------------------------------------------
std::vector<StateBlock> SolidStateBlocks(const BasisSet::Complex_BS& bs, SpinGroup g, const WaveFunction::cWaveFunction* wf)
{
    std::vector<StateBlock> out;
    for (Spin s : SpinIrreps(g))
        for (size_t i=0;i<bs.GetNumIBS();++i)
        {
            StateBlock b;
            if (const auto* rb=bs.GetRealIBS(i)) { b.real=rb;    b.irrep=rb->GetIrrep(s); }
            else                                 { b.cplx=bs[i]; b.irrep=b.cplx->GetIrrep(s); }
            const rvec3_t k=Symmetry::Lattice_3D::Getk(*b.irrep.sym);
            b.extras["k"]={k.x, k.y, k.z};
            b.weight=b.irrep.sym->GetWeight();
            if (wf) b.orbitals=wf->GetOrbitals(b.irrep);
            out.push_back(std::move(b));
        }
    return out;
}

namespace
{
bool Same3(const std::vector<double>& a, const std::vector<double>& b)
{
    if (a.size()!=3 || b.size()!=3) return false;
    return std::fabs(a[0]-b[0])+std::fabs(a[1]-b[1])+std::fabs(a[2]-b[2]) <= 1e-10;
}
} //anonymous namespace

//---------------------------------------------------------------------------------------------------
void AddSolidKeys(StateFingerprint& fp, const rmat3d_t& M, const BasisSet::Complex_BS& bs, SpinGroup g, const RunIdentity& id)
{
    const std::vector<StateBlock> blocks=SolidStateBlocks(bs, g);
    // ---- REFUSE: the saved D is a matrix over these functions, at these k (the structure, electron and block keys are the base's) ----
    {
        std::vector<double> A;                                   // rows = a_1, a_2, a_3 (Matrix3D is 1-based)
        for (size_t col=1;col<=3;col++) for (size_t row=1;row<=3;row++) A.push_back(M(row,col));
        fp.refuseNum["structure.cell"]=A;
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
        // The ORBITAL basis, shell by shell: every Bloch block is the Bloch sum of these shells (AoShellSource's own contract), so
        // block 0 speaks for all of them.
        if (const auto* rb=bs.GetRealIBS(0)) AddShells(fp, std::vector<const BasisSet::Orbital_1E_IBS<double>*>{rb});
        else                                 AddShells(fp, *bs[0]);
    }
    {   // the solid's own block identity: each block's k and BZ weight (ms, n and the labels are the shared keys)
        std::vector<double> k, w;
        for (const StateBlock& b : blocks)
        {
            const std::vector<double>& kk=b.extras.at("k");
            k.insert(k.end(), kk.begin(), kk.end());
            w.push_back(b.weight);
        }
        fp.refuseNum["blocks.k"]=k;
        fp.refuseNum["blocks.weight"]=w;
    }

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
    fp.warmNum["grids.densityEcut"]={id.densityEcut};
    fp.warmNum["grids.cutoffFactor"]={id.cutoffFactor};
    fp.warmStr["grids.xcMesh"]=id.xcMesh;
    fp.warmNum["ortho.kind"]={double(id.ortho)};
    fp.warmNum["ortho.tol"]={id.orthoTol};
    // The numerics that move the ground state (gap 2, 2026-10-30).  Each key exists only where it means something, so a
    // plane-wave state carries no GPW key and a uniform-mesh state no Becke eps.
    if (!id.vxcFit.empty())     fp.warmStr["xc.fit"]=id.vxcFit;
    if (id.beckeEps>0.0)        fp.warmNum["grids.beckeEps"]={id.beckeEps};
    if (!id.raster.empty())     fp.warmStr["grids.raster"]=id.raster;
    if (id.ladderFactor>0.0)    fp.warmNum["grids.ladderFactor"]={id.ladderFactor};
    for (const auto& [k,v] : id.gpwTolerances) fp.warmNum["gpw."+k]=v;
    for (const auto& [knob,on] : id.policy)    fp.warmStr[knob]=on ? "on" : "off";
}

//! A state written before C4 step 5 filed the functional, the atom positions and the shell centres under REFUSE and had no
//! pp.set (nor, before 2026-10-30, basis.family): move them to the class they belong to now, and take the run's pp.set for a state that never recorded one (its
//! species string, still REFUSE, carries the q per element).  Additive: an old file keeps loading.
StateFingerprint MigrateSolidFingerprint(StateFingerprint s, const StateFingerprint& now)
{
    // The block irrep labels (qchem-state 2.0, 2026-11-09): a 1.x state has none; its blocks are already identified by k, ms and n.
    if (!s.refuseStr.count("blocks.irrep") && now.refuseStr.count("blocks.irrep")) s.refuseStr["blocks.irrep"]=now.refuseStr.at("blocks.irrep");
    for (const char* k : {"structure.positions", "basis.shellCenters"})
        if (auto it=s.refuseNum.find(k); it!=s.refuseNum.end()) { s.warmNum[k]=std::move(it->second); s.refuseNum.erase(it); }
    if (auto it=s.refuseStr.find("functional"); it!=s.refuseStr.end()) { s.warmStr["functional"]=std::move(it->second); s.refuseStr.erase(it); }
    if (!s.refuseStr.count("pp.set") && now.refuseStr.count("pp.set")) s.refuseStr["pp.set"]=now.refuseStr.at("pp.set");
    // grids.raster's old name for the ball's own rasters (2026-11-04): the same run, so the same value.
    if (auto it=s.warmStr.find("grids.raster"); it!=s.warmStr.end() && it->second=="BallOnly") it->second=BallRastersLabel;
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
    // The rung gate's calibration (2026-11-05): a state that records none was built on the default, whatever this run states.
    if (!s.warmNum.count("gpw.rungGateC") && now.warmNum.count("gpw.rungGateC")) s.warmNum["gpw.rungGateC"]={BasisSet::Gaussian::kRungGateC};
    // The pair rule's stiffness (V-BS-RASTER step 2, 2026-11-07): before it was a tolerance it was the evaluator's constant, the default.
    if (!s.warmNum.count("gpw.relCutoffSafety") && now.warmNum.count("gpw.relCutoffSafety"))
        s.warmNum["gpw.relCutoffSafety"]={BasisSet::Gaussian::kRelCutoffSafety};
    // The raster demands (C2 step 3, 2026-10-31): a state written before them sampled on a NAMED raster.  The ball's own raster has
    // every demand 0; the retired AliasFree has no demand equivalent, so its keys stay absent and differ (warm).
    // The accuracy target itself moves no state (the demands it derives do): an older state takes the run's.
    if (!s.warmNum.count("gpw.rasterX") && now.warmNum.count("gpw.rasterX"))
    {
        auto r=s.warmStr.find("grids.raster");
        if (r!=s.warmStr.end() && r->second==BallRastersLabel)
            for (const char* k : {"gpw.rasterX","gpw.vlocRasterX","gpw.xcRasterX"}) s.warmNum[k]={0.0};
        if (now.warmNum.count("gpw.accuracyEps")) s.warmNum["gpw.accuracyEps"]=now.warmNum.at("gpw.accuracyEps");
    }
    // The policy routes (2026-10-30): likewise, all or none -- a state that records any route was written with the whole table.
    auto isPolicy=[](const auto& e){ return e.first.rfind("policy.",0)==0; };
    if (std::none_of(s.warmStr.begin(), s.warmStr.end(), isPolicy))
        for (const auto& e : now.warmStr) if (isPolicy(e)) s.warmStr.insert(e);
    return s;
}


Outcome<std::vector<std::string>, RestartRefusal> CompareSolidFingerprints(const StateFingerprint& saved, const StateFingerprint& now)
{
    return CompareFingerprints(saved, now, MigrateSolidFingerprint);
}

SavedState SnapshotState(const BasisSet::Complex_BS& bs, const WaveFunction::cWaveFunction& wf, bool energyWeighted)
{
    return SnapshotState(SolidStateBlocks(bs, wf.GetSpinGroup(), &wf), "solid", energyWeighted);
}

void LowdinCarry(SavedState& carry, const BasisSet::Complex_BS& bs, SpinGroup g) {LowdinCarry(carry, SolidStateBlocks(bs, g));}

Outcome<SavedState, RestartRefusal> ReadSolidState(const std::string& path)
{
    using O=Outcome<SavedState, RestartRefusal>;
    auto st=ReadState(path);
    if (!st) return st;
    if (st->kind!="solid")
        return O::Fail({RestartRefusal::Why::Kind, "'"+path+"' holds a "+st->kind+" state; a solid run restarts only from a solid one"});
    return st;
}

//---------------------------------------------------------------------------------------------------
//  RESTORE
//---------------------------------------------------------------------------------------------------
Outcome<std::unique_ptr<ChargeDensity::cDM_CD>, RestartRefusal>
RestoreDensity(const SavedState& st, const BasisSet::Complex_BS& bs, SpinGroup g)
{
    using O=Outcome<std::unique_ptr<ChargeDensity::cDM_CD>, RestartRefusal>;
    using namespace ChargeDensity;
    const std::vector<StateBlock> refs=SolidStateBlocks(bs, g);
    // Pass 1: every block's weighted D, as saved -- the shared checks (counts, sizes, spins, realness), then the solid's k.
    auto checked=CheckedBlockDensities(st, refs);
    if (!checked) return O::Fail(checked.Error());
    std::vector<hmat_t<double>> Dr(refs.size());
    std::vector<hmat_t<dcmplx>> Dc(refs.size());
    for (size_t b=0;b<refs.size();b++)
    {
        const StateBlock& r=refs[b];
        const SavedBlock& sb=st.blocks[b];
        const auto sk=sb.extras.find("k");
        if (sk==sb.extras.end() || !Same3(sk->second, r.extras.at("k")))
            throw std::logic_error("RestoreDensity: block "+std::to_string(b)+" is not at this run's k -- the fingerprint should have refused");
        const double w=r.weight;
        const size_t n=r.n();
        const mat_t<dcmplx>& D=(*checked)[b];
        if (r.real)
        {
            hmat_t<double> Dw(n);
            for (size_t a=0;a<n;a++) for (size_t c=a;c<n;c++) Dw(a,c)=w*D(a,c).real();
            Dr[b]=std::move(Dw);
        }
        else
        {
            hmat_t<dcmplx> Dw(n);
            for (size_t a=0;a<n;a++) for (size_t c=a;c<n;c++) Dw(a,c)=w*D(a,c);
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
        {
            const auto* imp=bs.GetImposition();
            if (!imp) break;
            if (refs[b].real) Symmetry::ProjectSpinPair(*imp, *refs[b].real, Dr[b], Dr[b+N]);
            else              Symmetry::ProjectSpinPair(*imp, *refs[b].cplx, Dc[b], Dc[b+N]);
        }
    }
    auto cd=std::make_unique<tComposite_CD<dcmplx>>(bs.GetImposition() ? bs.GetImposition()->DensityGroup() : nullptr);
    for (size_t b=0;b<refs.size();b++)
    {
        const StateBlock& r=refs[b];
        if (r.real) cd->Insert(std::unique_ptr<tDM_CD<double>>(IrrepCD_Factory<double>(Dr[b], r.real, r.irrep)), r.irrep);
        else        cd->Insert(std::unique_ptr<tDM_CD<dcmplx>>(IrrepCD_Factory<dcmplx>(Dc[b], r.cplx, r.irrep)), r.irrep);
    }
    return O::Ok(std::unique_ptr<cDM_CD>(std::move(cd)));
}

} //namespace qchem
