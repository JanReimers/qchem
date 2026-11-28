// File: src/Hamiltonian/Internal/Imp/Mixed.C  The linear mixer term's bodies (V-HAM-COMPOSE step 3, lane S 2026-11-26); see the interface.
module;
#include <algorithm>
#include <cmath>
#include <memory>
#include <optional>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.Hamiltonian.Internal.Mixed;
import qchem.ChargeDensity;   // SpinIrrepsOf
import qchem.ChargeDensity.DensitySampler;   // the mesh a child record carries
import qchem.Streamable;
import qchem.Blaze;

namespace qchem::Hamiltonian
{

template <class T> tTerm<T> Mix(std::vector<tWeightedTerm<T>> parts)
{
    if (parts.empty()) throw std::invalid_argument("Hamiltonian::Mix: a mix of no terms is not a term");
    std::shared_ptr<const ChargeDensity::tDensitySampler<T>> mesh;
    bool anyHF=false, anyDynamic=false;
    for (auto& p : parts)
    {
        if (!std::isfinite(p.w)) throw std::invalid_argument("Hamiltonian::Mix: a weight is not finite");
        const int filled=int(bool(p.term.statik))+int(bool(p.term.dynamic))+int(bool(p.term.hf));
        if (filled!=1) throw std::logic_error("Hamiltonian::Mix: a child record holds "+std::to_string(filled)+" terms (exactly one expected)");
        anyHF      = anyHF      || bool(p.term.hf);
        anyDynamic = anyDynamic || bool(p.term.dynamic);
        if (p.term.mesh)
        {
            if (mesh) throw std::logic_error("Hamiltonian::Mix: two children bring a nonlinear integration mesh -- one Hamiltonian, one mesh "
                                             "(V-NONLINEAR-MESH): the XC parts on one mesh belong in ONE XC term (CompositeExFunctional)");
            mesh=std::move(p.term.mesh);
        }
    }
    tMixedChildren<T> children(std::move(parts));
    tTerm<T> out;
    if      (anyHF)      out=tTerm<T>::Of(new tMixed_HF_HT<T>     (std::move(children)));
    else if (anyDynamic) out=tTerm<T>::Of(new tMixed_Dynamic_HT<T>(std::move(children)));
    else                 out=tTerm<T>::Of(new tMixed_Static_HT<T> (std::move(children)));
    out.mesh=std::move(mesh);
    return out;
}

template <class T> tMixedChildren<T>::tMixedChildren(std::vector<tWeightedTerm<T>>&& parts)
{
    for (auto& p : parts)
    {
        if (p.term.statik)  itsStatic .push_back({p.w, std::move(p.term.statik)});
        if (p.term.dynamic) itsDynamic.push_back({p.w, std::move(p.term.dynamic)});
        if (p.term.hf)      itsHF     .push_back({p.w, std::move(p.term.hf)});
    }
}

template <class T> void tMixedChildren<T>::RefreshState(const tbs_t<T>* wholeBasis, const tChargeDensity<T>* cd) const
{
    if (!cd) return;
    for (const auto& c : itsDynamic) c.t->RefreshForDensity(cd);
    for (const auto& c : itsHF)      c.t->RefreshForDensity(wholeBasis, cd);
}
template <class T> void tMixedChildren<T>::BuildBlocks(const tbs_t<T>* bs, const tChargeDensity<T>* cd) const
{
    for (const auto& c : itsStatic)  c.t->BuildBlocks(bs, cd);
    for (const auto& c : itsDynamic) c.t->BuildBlocks(bs, cd);
}

template <class T> hmat_t<T> tMixedChildren<T>::Sum(const tobs_t<T>* bs, const Spin& s, DensityID id, const tbs_t<T>* wholeBasis, bool energy) const
{
    if (!itsHF.empty() && !wholeBasis) throw std::logic_error("Hamiltonian mixer: an HF child needs the whole (composite) basis");
    hmat_t<T> H=blazem::zeroH<T>(bs->GetNumFunctions());
    for (const auto& c : itsStatic)  H+=c.w*c.t->GetMatrix(bs,s);
    for (const auto& c : itsDynamic) H+=c.w*(energy ? c.t->GetEMatrix(bs,s,id) : c.t->GetMatrix(bs,s,id));
    for (const auto& c : itsHF)      H+=c.w*c.t->GetMatrix(bs,s,id,wholeBasis);
    return H;
}

namespace
{
//! \a c scaled by \a w into \a e: every contribution, with its role and expectation; "Exc" through the one XC door with the child's pin-28
//! certificate; the diagnostics; and the charge accounting the child wrote (a property of the density, not of the term: never scaled).
void AddScaled(EnergyBreakdown& e, const EnergyBreakdown& c, double w)
{
    for (const auto& [n,t] : c.Terms())
    {
        const std::optional<double> trDV = t.TrDV ? std::optional<double>(w**t.TrDV) : std::nullopt;
        if (n=="Exc") e.AddXC(w*t.E, trDV, c.charge.xcRho, c.charge.xcRhoCertified);
        else          e.Add  (n, w*t.E, t.role, trDV);
    }
    for (const auto& [n,d] : c.Diagnostics()) e.AddDiagnostic(n, w*d);
    e.charge.lost+=c.charge.lost;
    if (e.charge.siteMoments.size()==0) e.charge.siteMoments=c.charge.siteMoments;
    if (!e.charge.xcGroupResidual) e.charge.xcGroupResidual=c.charge.xcGroupResidual;
    if (!e.charge.samplingOverTau) e.charge.samplingOverTau=c.charge.samplingOverTau;
}
template <class P> std::string NameOf(const P& t) {std::ostringstream os; os<<t; std::string l=os.str(); return l.substr(0,l.find('\n'));}
} //anon

template <class T> void tMixedChildren<T>::Energy(EnergyBreakdown& e, const tDM_CD<T>* cd) const
{
    auto one=[&](const auto& c) {EnergyBreakdown ec; c.t->GetEnergy(ec,cd); AddScaled(e,ec,c.w);};
    for (const auto& c : itsStatic)  one(c);
    for (const auto& c : itsDynamic) one(c);
    for (const auto& c : itsHF)      one(c);
}

// The forces and virials: every child must carry the face (abstract -> abstract cross-cast); one that does not is NAMED, never a silent zero.
namespace
{
template <class Face, class C, class F> void ForEachChild(const C& c, F&& f)
{
    for (const auto& x : c)
    {
        const auto* g=dynamic_cast<const Face*>(x.t.get());
        f(x.w, g, *x.t);
    }
}
} //anon

template <class T> std::string tMixedChildren<T>::MissingGradient(const tbs_t<T>* bs, const Structure* st) const
{
    std::string why;
    auto ask=[&](double, const tGeometryGradient_HT<T>* g, const auto& t)
    {
        if (!why.empty()) return;
        why = g ? g->MissingGradient(bs,st) : "Hamiltonian mixer: its child '"+NameOf(t)+"' has no force (tGeometryGradient_HT)";
    };
    ForEachChild<tGeometryGradient_HT<T>>(itsStatic,ask); ForEachChild<tGeometryGradient_HT<T>>(itsDynamic,ask); ForEachChild<tGeometryGradient_HT<T>>(itsHF,ask);
    return why;
}
template <class T> void tMixedChildren<T>::AddForces(rvec3vec_t& F, const tDM_CD<T>* D, const Structure* st) const
{
    auto add=[&](double w, const tGeometryGradient_HT<T>* g, const auto& t)
    {
        if (!g) throw std::logic_error("Hamiltonian mixer: its child '"+NameOf(t)+"' has no force (MissingGradient says so first)");
        rvec3vec_t f(F.size(), rvec3_t(0.0,0.0,0.0));
        g->AddForces(f,D,st);
        for (size_t a=0;a<F.size();a++) F[a]+=w*f[a];
    };
    ForEachChild<tGeometryGradient_HT<T>>(itsStatic,add); ForEachChild<tGeometryGradient_HT<T>>(itsDynamic,add); ForEachChild<tGeometryGradient_HT<T>>(itsHF,add);
}
template <class T> std::string tMixedChildren<T>::MissingVirial(const tbs_t<T>* bs, const Structure* st) const
{
    std::string why;
    auto ask=[&](double, const tStrainGradient_HT<T>* g, const auto& t)
    {
        if (!why.empty()) return;
        why = g ? g->MissingVirial(bs,st) : "Hamiltonian mixer: its child '"+NameOf(t)+"' has no virial (tStrainGradient_HT)";
    };
    ForEachChild<tStrainGradient_HT<T>>(itsStatic,ask); ForEachChild<tStrainGradient_HT<T>>(itsDynamic,ask); ForEachChild<tStrainGradient_HT<T>>(itsHF,ask);
    return why;
}
template <class T> void tMixedChildren<T>::AddVirial(rmat3d_t& V, const tDM_CD<T>* D, const Structure* st) const
{
    auto add=[&](double w, const tStrainGradient_HT<T>* g, const auto& t)
    {
        if (!g) throw std::logic_error("Hamiltonian mixer: its child '"+NameOf(t)+"' has no virial (MissingVirial says so first)");
        rmat3d_t v(0,0,0, 0,0,0, 0,0,0);   // NOT rmat3d_t(): Matrix3D default-constructs the IDENTITY
        g->AddVirial(v,D,st);
        V+=w*v;
    };
    ForEachChild<tStrainGradient_HT<T>>(itsStatic,add); ForEachChild<tStrainGradient_HT<T>>(itsDynamic,add); ForEachChild<tStrainGradient_HT<T>>(itsHF,add);
}

template <class T> bool tMixedChildren<T>::RequiresDensityMatrix() const
{
    auto needs=[](const auto& t, bool byDefault)
    {
        const auto* r=dynamic_cast<const DensityMatrixRequirement*>(t.get());   // abstract -> abstract
        return r ? r->RequiresDensityMatrix() : byDefault;
    };
    for (const auto& c : itsDynamic) if (needs(c.t,false)) return true;
    for (const auto& c : itsHF)      if (needs(c.t,true))  return true;   // tHamiltonianImp::Add(HF)'s default
    return false;
}
template <class T> bool tMixedChildren<T>::IsVirialValid() const
{
    for (const auto& c : itsStatic)  if (!c.t->IsVirialValid()) return false;
    for (const auto& c : itsDynamic) if (!c.t->IsVirialValid()) return false;
    for (const auto& c : itsHF)      if (!c.t->IsVirialValid()) return false;
    return true;
}
template <class T> bool tMixedChildren<T>::PreservesReal() const
{
    for (const auto& c : itsStatic)  if (!c.t->PreservesReal()) return false;
    for (const auto& c : itsDynamic) if (!c.t->PreservesReal()) return false;
    for (const auto& c : itsHF)      if (!c.t->PreservesReal()) return false;
    return true;
}
template <class T> bool tMixedChildren<T>::IsRelativistic() const
{
    for (const auto& c : itsStatic)  if (c.t->IsRelativistic()) return true;
    for (const auto& c : itsDynamic) if (c.t->IsRelativistic()) return true;
    for (const auto& c : itsHF)      if (c.t->IsRelativistic()) return true;
    return false;
}

template <class T> std::ostream& tMixedChildren<T>::Write(std::ostream& os) const
{
    os << "Linear mix of " << itsStatic.size()+itsDynamic.size()+itsHF.size() << " terms:";
    auto one=[&](const auto& c) {os << "\n      " << c.w << " x " << NameOf(*c.t);};
    for (const auto& c : itsStatic)  one(c);
    for (const auto& c : itsDynamic) one(c);
    for (const auto& c : itsHF)      one(c);
    return os << std::endl;
}

// ---- the dynamic mixer: the children's blocks, then the two sums, stamped ------------------------------------------------------------
template <class T> void tMixed_Dynamic_HT<T>::BuildBlocks(const tbs_t<T>* bs, const tChargeDensity<T>* cd) const
{
    itsChildren.BuildBlocks(bs, cd);                      // the children's BLOCKS (their state was refreshed by RefreshForDensity(cd))
    if (!cd || !bs) return;
    const auto refreshV=this->itsCache.BeginRefresh();    // the potential sum, into the base's cache ...
    const auto refreshE=itsECache.BeginRefresh();         // ... and the energy sum
    this->itsCacheID=DensityID(cd);
    const DensityID id(cd);
    for (size_t i=0;i<bs->GetNumIBS();++i)
    {
        if (bs->GetRealIBS(i)) continue;                  // the real cache's (a mixer has no real-block face: see the interface)
        if (const tobs_t<T>* b=(*bs)[i])
            for (const Spin& s : ChargeDensity::SpinIrrepsOf(cd))
            {
                refreshV.Fill(b->GetIrrep(s), [&]{return itsChildren.Sum(b,s,id,nullptr,false);});
                refreshE.Fill(b->GetIrrep(s), [&]{return itsChildren.Sum(b,s,id,nullptr,true);});
            }
    }
}

template <class T> const hmat_t<T>& tMixed_Dynamic_HT<T>::GetEMatrix(const tobs_t<T>* bs, const Spin& s, DensityID id) const
{
    if (id.IsNone() || !(id==this->itsCacheID))
        throw std::logic_error("Hamiltonian mixer GetEMatrix for a density it was not refreshed for (doc/KAxisPlan.md P2)");
    const hmat_t<T>* block=itsECache.Read(bs->GetIrrep(s));
    if (!block) throw std::logic_error("Hamiltonian mixer GetEMatrix: a block the refresh did not build");
    return *block;
}

// ---- the HF mixer: latch the basis, refresh the children, sum every block of every spin irrep cd resolves --------------------------------
template <class T> void tMixed_HF_HT<T>::RefreshForDensity(const tbs_t<T>* wholeBasis, const tChargeDensity<T>* cd) const
{
    itsChildren.RefreshState(wholeBasis, cd);             // each child latches the whole basis itself (a null one: keep the latched)
    itsChildren.BuildBlocks (wholeBasis, cd);
    if (!cd || !wholeBasis) return;                       // an energy evaluation: the children are refreshed, no Fock blocks to build
    const auto refresh=itsCache.BeginRefresh();
    itsCacheID=DensityID(cd);
    const DensityID id(cd);
    for (size_t i=0;i<wholeBasis->GetNumIBS();++i)
    {
        if (wholeBasis->GetRealIBS(i)) continue;
        if (const tobs_t<T>* b=(*wholeBasis)[i])
            for (const Spin& s : ChargeDensity::SpinIrrepsOf(cd))
                refresh.Fill(b->GetIrrep(s), [&]{return itsChildren.Sum(b,s,id,wholeBasis,false);});
    }
}

template <class T> const hmat_t<T>& tMixed_HF_HT<T>::GetMatrix(const tobs_t<T>* bs, const Spin& s, DensityID id, const tbs_t<T>*) const
{
    if (id.IsNone() || !(id==itsCacheID))
        throw std::logic_error("Hamiltonian mixer GetMatrix for a density it was not refreshed for -- call "
                               "tHamiltonian::RefreshForDensity(bs, cd) first (doc/KAxisPlan.md P2)");
    const hmat_t<T>* block=itsCache.Read(bs->GetIrrep(s));
    if (!block) throw std::logic_error("Hamiltonian mixer GetMatrix: a block the refresh did not build (a spin the density does not resolve)");
    return *block;
}

template tTerm<double> Mix(std::vector<tWeightedTerm<double>>);
template tTerm<dcmplx> Mix(std::vector<tWeightedTerm<dcmplx>>);
template class tMixedChildren<double>;    template class tMixedChildren<dcmplx>;
template class tMixed_Static_HT<double>;  template class tMixed_Static_HT<dcmplx>;
template class tMixed_Dynamic_HT<double>; template class tMixed_Dynamic_HT<dcmplx>;
template class tMixed_HF_HT<double>;      template class tMixed_HF_HT<dcmplx>;

} //namespace
