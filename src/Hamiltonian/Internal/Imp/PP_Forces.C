// File: Hamiltonian/Internal/Imp/PP_Forces.C  The two pseudopotential terms' forces (doc/GeometryOptimizationPlan.md §13.2; lane B
// 2026-11-08).  Its own TU so the force work stays out of Imp/PP_Local.C, Imp/PP_NonLocal.C and the term files lane S's sprint edits.
//
// Both terms are STATIC and LINEAR in D, E = Tr(D V) with V = the closed-form species-field (local) or projector (KB) matrix, so
// ∂E/∂R_A = Tr(D ∂V/∂R_A), the basis centres AND the species' sites moving: each block's from its Orbital_PP_IBS pullback, summed by the
// density walk.  Only the CLOSED-FORM route has one.  The mesh-quadrature route (a block that cannot assemble the field) would need the
// mesh's own motion, and the radial (atomic) route has one centre and no geometry: both are MissingGradients() entries.
module;
#include <stdexcept>
#include <string>
module qchem.Hamiltonian.Internal.Terms;
import qchem.ChargeDensity;
import qchem.Structure;

namespace qchem::Hamiltonian
{

namespace
{
const BasisSet::Orbital_PP_IBS<double>* ClosedForm(const robs_t* b) {return dynamic_cast<const BasisSet::Orbital_PP_IBS<double>*>(b);}

void Subtract(rvec3vec_t& F, const rvec3vec_t& g, const char* who)
{
    if (g.size()!=F.size()) throw std::logic_error(std::string(who)+"::AddForces: the pullback and the force array differ in atom count");
    for (size_t A=0;A<F.size();A++) F[A]-=g[A];
}
} // namespace

// --- the local part: Tr(D V_loc), V_loc = MakeSpeciesFieldMatrix(st, v_loc, Full) on the closed-form route -----------------------------
std::string PP_Local::MissingGradient(const tbs_t<double>* bs, const Structure*) const
{
    if (!bs) return "PP_Local: no basis to differentiate";
    for (const auto* b : bs->template Iterate<robs_t>())
    {
        const auto* pp=ClosedForm(b);
        if (!pp || !pp->CanAssembleField(*itsVloc))
            return "PP_Local: the basis block "+b->Name()+" takes the mesh-quadrature route, which has no geometry gradient (the mesh moves)";
        if (!pp->CanPullbackSpeciesField(*itsVloc)) return "PP_Local: the basis block "+b->Name()+" gives no species-field pullback";
    }
    return "";
}

rvec3vec_t PP_Local::Pullback(const robs_t* bs, const Spin&, const hmat_t<double>& D, const Structure*) const
{
    const auto* pp=ClosedForm(bs);
    if (!pp) throw std::logic_error("PP_Local::Pullback: the block has no closed-form PP face (MissingGradients() names it before a force is asked)");
    const auto r=pp->MakeSpeciesFieldPullback(theStructure.get(), *itsVloc, BasisSet::FieldRange::Full, D);
    if (!r) throw std::logic_error("PP_Local::Pullback: "+r.Error().why+" (MissingGradients() names it before a force is asked)");
    return *r;
}

void PP_Local::AddForces(rvec3vec_t& F, const rDM_CD* D, const Structure* st) const {Subtract(F, D->DM_Pullback(this, st), "PP_Local");}

// --- the KB part: Tr(D V_NL), V_NL = MakeProjectorMatrix(st, projectors) on the closed-form route -------------------------------------
std::string PP_NonLocal::MissingGradient(const tbs_t<double>* bs, const Structure*) const
{
    if (!bs) return "PP_NonLocal: no basis to differentiate";
    for (const auto* b : bs->template Iterate<robs_t>())
    {
        if (dynamic_cast<const BasisSet::ImplicitAngular_IBS*>(b))
            return "PP_NonLocal: the basis block "+b->Name()+" takes the radial (atomic) route: one centre, no geometry gradient";
        const auto* pp=ClosedForm(b);
        if (!pp || !pp->CanAssembleProjectors(*itsSep))
            return "PP_NonLocal: the basis block "+b->Name()+" takes the mesh-quadrature route, which has no geometry gradient (the mesh moves)";
        if (!pp->CanPullbackProjectors(*itsSep)) return "PP_NonLocal: the basis block "+b->Name()+" gives no projector pullback";
    }
    return "";
}

rvec3vec_t PP_NonLocal::Pullback(const robs_t* bs, const Spin&, const hmat_t<double>& D, const Structure*) const
{
    const auto* pp=ClosedForm(bs);
    if (!pp || dynamic_cast<const BasisSet::ImplicitAngular_IBS*>(bs))
        throw std::logic_error("PP_NonLocal::Pullback: the block has no closed-form PP route (MissingGradients() names it before a force is asked)");
    const auto r=pp->MakeProjectorPullback(theStructure.get(), *itsSep, D);
    if (!r) throw std::logic_error("PP_NonLocal::Pullback: "+r.Error().why+" (MissingGradients() names it before a force is asked)");
    return *r;
}

void PP_NonLocal::AddForces(rvec3vec_t& F, const rDM_CD* D, const Structure* st) const {Subtract(F, D->DM_Pullback(this, st), "PP_NonLocal");}

} // namespace qchem::Hamiltonian
