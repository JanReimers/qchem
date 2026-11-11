// File: Hamiltonian/Internal/Imp/PWTerms_PP_Virials.C  The plane-wave PP terms' VIRIALS (doc/Records/CellHalf_Design.md §5; user ruling
// §6.2, 2026-11-10; lane H).  Its own TU, beside Imp/PWTerms_PP_Forces.C, so the strain work stays out of the term files lane S edits.
//
// Each of Ven_PP_Short / Ven_PP_Long / Ven_PP_NonLocal is STATIC and LINEAR in D, E = Re Tr(D† V) (+ a local half's G = 0 alignment
// E_α = N·α/Ω), so at fixed D, ∂E/∂ε_ab = Σ_blocks MakeSpeciesFieldVirial / MakeProjectorVirial (the basis owns Ω, G and the structure
// factor) − δ_ab E_α (α/Ω scales as 1/Ω; N is fixed).  A basis whose virial refuses is a MissingVirials() entry, never a silent zero.
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
template <class U> const BasisSet::Orbital_PP_IBS<U>* PP(const tobs_t<U>* b) {return dynamic_cast<const BasisSet::Orbital_PP_IBS<U>*>(b);}
//! Empty when every block, a real TRIM one included, can give \a can; else why not, prefixed by \a who.
template <class Can> std::string Missing(const tbs_t<dcmplx>* bs, const char* who, const char* what, Can&& can)
{
    if (!bs) return std::string(who)+": no basis to differentiate";
    std::string why;
    ForEachBlock(bs, [&](const auto* b)
    {
        const auto* pp=PP(b);
        if (why.empty() && (!pp || !can(*pp))) why=std::string(who)+": the basis block "+b->Name()+" gives no "+what;
    });
    return why;
}
rmat3d_t Value(const Outcome<rmat3d_t,AD::NoGradient>& r, const char* who)
{
    if (!r) throw std::logic_error(std::string(who)+"::Virial: "+r.Error().why+" (MissingVirials() names it before a stress is asked)");
    return *r;
}
template <class U> const BasisSet::Orbital_PP_IBS<U>& Need(const tobs_t<U>* b, const char* who)
{
    const auto* pp=PP(b);
    if (!pp) throw std::logic_error(std::string(who)+"::Virial: the block has no species-field face");
    return *pp;
}
//! V += v − δ E_α: a block sum plus a local half's alignment (0 for the KB term and for a finite structure).
void Add(rmat3d_t& V, const rmat3d_t& v, double Ealpha)
{
    for (size_t a=1;a<=3;a++) for (size_t b=1;b<=3;b++) V(a,b)+=v(a,b);
    for (size_t a=1;a<=3;a++) V(a,a)-=Ealpha;
}
} // namespace

std::string Ven_PP_Short::MissingVirial(const tbs_t<dcmplx>* bs, const Structure*) const
{return Missing(bs, "Ven_PP_Short", "species-field virial", [&](const auto& pp){return pp.CanVirialSpeciesField(*itsLocal);});}
std::string Ven_PP_Long::MissingVirial(const tbs_t<dcmplx>* bs, const Structure*) const
{return Missing(bs, "Ven_PP_Long", "species-field virial", [&](const auto& pp){return pp.CanVirialSpeciesField(*itsLocal);});}
std::string Ven_PP_NonLocal::MissingVirial(const tbs_t<dcmplx>* bs, const Structure*) const
{return Missing(bs, "Ven_PP_NonLocal", "projector virial", [&](const auto& pp){return pp.CanVirialProjectors(*itsSep);});}

rmat3d_t Ven_PP_Short::Virial(const cobs_t* bs, const Spin&, const chmat_t& D, const Structure*) const
{return Value(Need(bs,"Ven_PP_Short").MakeSpeciesFieldVirial(&*theStructure, *itsLocal, BasisSet::FieldRange::Short, D), "Ven_PP_Short");}
rmat3d_t Ven_PP_Short::Virial(const robs_t* bs, const Spin&, const hmat_t<double>& D, const Structure*) const
{return Value(Need(bs,"Ven_PP_Short").MakeSpeciesFieldVirial(&*theStructure, *itsLocal, BasisSet::FieldRange::Short, D), "Ven_PP_Short");}
rmat3d_t Ven_PP_Long::Virial(const cobs_t* bs, const Spin&, const chmat_t& D, const Structure*) const
{return Value(Need(bs,"Ven_PP_Long").MakeSpeciesFieldVirial(&*theStructure, *itsLocal, BasisSet::FieldRange::Long, D), "Ven_PP_Long");}
rmat3d_t Ven_PP_Long::Virial(const robs_t* bs, const Spin&, const hmat_t<double>& D, const Structure*) const
{return Value(Need(bs,"Ven_PP_Long").MakeSpeciesFieldVirial(&*theStructure, *itsLocal, BasisSet::FieldRange::Long, D), "Ven_PP_Long");}
rmat3d_t Ven_PP_NonLocal::Virial(const cobs_t* bs, const Spin&, const chmat_t& D, const Structure*) const
{return Value(Need(bs,"Ven_PP_NonLocal").MakeProjectorVirial(&*theStructure, *itsSep, D), "Ven_PP_NonLocal");}
rmat3d_t Ven_PP_NonLocal::Virial(const robs_t* bs, const Spin&, const hmat_t<double>& D, const Structure*) const
{return Value(Need(bs,"Ven_PP_NonLocal").MakeProjectorVirial(&*theStructure, *itsSep, D), "Ven_PP_NonLocal");}

void Ven_PP_Short::AddVirial   (rmat3d_t& V, const cDM_CD* D, const Structure* st) const {Add(V, D->DM_Virial(this, st), D->GetTotalCharge()*itsAlphaZ);}
void Ven_PP_Long::AddVirial    (rmat3d_t& V, const cDM_CD* D, const Structure* st) const {Add(V, D->DM_Virial(this, st), D->GetTotalCharge()*itsAlphaZ);}
void Ven_PP_NonLocal::AddVirial(rmat3d_t& V, const cDM_CD* D, const Structure* st) const {Add(V, D->DM_Virial(this, st), 0.0);}

} // namespace qchem::Hamiltonian
