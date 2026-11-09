// File: Hamiltonian/Internal/Imp/PWTerms_PP_Forces.C  The plane-wave PP terms' forces: the STRUCTURE-FACTOR Hellmann-Feynman term
// (doc/GeometryOptimizationPlan.md §13.2; lane B 2026-11-08).  Its own TU so the force work stays out of Imp/PWTerms_PP.C and the term
// files lane S's sprint edits.
//
// Each of Ven_PP_Short / Ven_PP_Long / Ven_PP_NonLocal is STATIC and LINEAR in D, E = Tr(D V) (+ the G = 0 alignment N·α, which a fixed
// D and a fixed cell leave constant), so ∂E/∂τ_a = Tr(D ∂V/∂τ_a): each block's from its Orbital_PP_IBS pullback (the basis owns the
// structure factor), summed by the density walk.  A real TRIM block inside a complex run is not walked yet (the composite's mixed-run
// pullback throws); a basis whose pullback refuses (GPW, a Bloch Gaussian) is a MissingGradients() entry.
module;
#include <stdexcept>
#include <string>
module qchem.Hamiltonian.Internal.PWTerms;
import qchem.ChargeDensity;
import qchem.BasisSet.Orbital_PP_IBS;
import qchem.Structure;

namespace qchem::Hamiltonian
{

namespace
{
const BasisSet::Orbital_PP_IBS<dcmplx>* PP(const cobs_t* b) {return dynamic_cast<const BasisSet::Orbital_PP_IBS<dcmplx>*>(b);}
//! Empty when every block can give \a can; else why not, prefixed by \a who.
template <class Can> std::string Missing(const tbs_t<dcmplx>* bs, const char* who, const char* what, Can&& can)
{
    if (!bs) return std::string(who)+": no basis to differentiate";
    for (const auto* b : bs->template Iterate<cobs_t>())
    {
        const auto* pp=PP(b);
        if (!pp || !can(*pp)) return std::string(who)+": the basis block "+b->Name()+" gives no "+what;
    }
    return "";
}
rvec3vec_t Value(const Outcome<rvec3vec_t,AD::NoGradient>& r, const char* who)
{
    if (!r) throw std::logic_error(std::string(who)+"::Pullback: "+r.Error().why+" (MissingGradients() names it before a force is asked)");
    return *r;
}
void Subtract(rvec3vec_t& F, const rvec3vec_t& g, const char* who)
{
    if (g.size()!=F.size()) throw std::logic_error(std::string(who)+"::AddForces: the pullback and the force array differ in atom count");
    for (size_t A=0;A<F.size();A++) F[A]-=g[A];
}
const BasisSet::Orbital_PP_IBS<dcmplx>& Need(const cobs_t* b, const char* who)
{
    const auto* pp=PP(b);
    if (!pp) throw std::logic_error(std::string(who)+"::Pullback: the block has no species-field face");
    return *pp;
}
} // namespace

std::string Ven_PP_Short::MissingGradient(const tbs_t<dcmplx>* bs, const Structure*) const
{return Missing(bs, "Ven_PP_Short", "species-field pullback", [&](const auto& pp){return pp.CanPullbackSpeciesField(*itsLocal);});}
std::string Ven_PP_Long::MissingGradient(const tbs_t<dcmplx>* bs, const Structure*) const
{return Missing(bs, "Ven_PP_Long", "species-field pullback", [&](const auto& pp){return pp.CanPullbackSpeciesField(*itsLocal);});}
std::string Ven_PP_NonLocal::MissingGradient(const tbs_t<dcmplx>* bs, const Structure*) const
{return Missing(bs, "Ven_PP_NonLocal", "projector pullback", [&](const auto& pp){return pp.CanPullbackProjectors(*itsSep);});}

rvec3vec_t Ven_PP_Short::Pullback(const cobs_t* bs, const Spin&, const chmat_t& D, const Structure*) const
{return Value(Need(bs,"Ven_PP_Short").MakeSpeciesFieldPullback(&*theStructure, *itsLocal, BasisSet::FieldRange::Short, D), "Ven_PP_Short");}
rvec3vec_t Ven_PP_Long::Pullback(const cobs_t* bs, const Spin&, const chmat_t& D, const Structure*) const
{return Value(Need(bs,"Ven_PP_Long").MakeSpeciesFieldPullback(&*theStructure, *itsLocal, BasisSet::FieldRange::Long, D), "Ven_PP_Long");}
rvec3vec_t Ven_PP_NonLocal::Pullback(const cobs_t* bs, const Spin&, const chmat_t& D, const Structure*) const
{return Value(Need(bs,"Ven_PP_NonLocal").MakeProjectorPullback(&*theStructure, *itsSep, D), "Ven_PP_NonLocal");}

void Ven_PP_Short::AddForces   (rvec3vec_t& F, const cDM_CD* D, const Structure* st) const {Subtract(F, D->DM_Pullback(this, st), "Ven_PP_Short");}
void Ven_PP_Long::AddForces    (rvec3vec_t& F, const cDM_CD* D, const Structure* st) const {Subtract(F, D->DM_Pullback(this, st), "Ven_PP_Long");}
void Ven_PP_NonLocal::AddForces(rvec3vec_t& F, const cDM_CD* D, const Structure* st) const {Subtract(F, D->DM_Pullback(this, st), "Ven_PP_NonLocal");}

} // namespace qchem::Hamiltonian
