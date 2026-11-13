// File: Calculation/Imp/StateView.C  A saved state, drawn: the levels, the orbitals and m(r) from the saved blocks.
module;
#include <algorithm>
#include <cmath>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.StateView;
import qchem.Types;
import qchem.Reporting;
import qchem.Symmetry.Orbital;          // Orbital_QNs: a level's label, spelled as the live orbital spells it
import qchem.BasisSet.Orbital_1E_IBS;   // a real block's functions at a point

namespace qchem
{
namespace
{
//! m(r) = rho_up(r) - rho_down(r), over two channel VIEWS the view's density owns.
struct ChannelDifference : public virtual ScalarFunction<double>
{
    ChannelDifference(const ScalarFunction<double>* up, const ScalarFunction<double>* down) : itsUp(up), itsDown(down) {}
    double  operator()(const rvec3_t& r) const override {return (*itsUp)(r)-(*itsDown)(r);}
    rvec3_t Gradient  (const rvec3_t& r) const override {return itsUp->Gradient(r)-itsDown->Gradient(r);}
    const ScalarFunction<double>* itsUp;
    const ScalarFunction<double>* itsDown;
};
//! m of an unpolarized state: a field, not an absence (the facades' ZeroField twin).
struct ZeroSpin : public virtual ScalarFunction<double>
{
    double  operator()(const rvec3_t&) const override {return 0.0;}
    rvec3_t Gradient  (const rvec3_t&) const override {return rvec3_t(0.0, 0.0, 0.0);}
};
//! One saved orbital of a REAL block: \f$\psi(r)=\sum_\mu c_\mu\chi_\mu(r)\f$ on the block's own functions.
struct SavedRealOrbital : public virtual ScalarFunction<double>
{
    SavedRealOrbital(const BasisSet::Orbital_1E_IBS<double>* bs, rvec_t c) : itsBS(bs), itsC(std::move(c)) {}
    double operator()(const rvec3_t& r) const override
    {
        const vec_t<double> chi=(*itsBS)(r);
        double psi=0.0;
        for (size_t a=0;a<itsC.size();a++) psi+=itsC[a]*chi[a];
        return psi;
    }
    rvec3_t Gradient(const rvec3_t& r) const override
    {
        const vec3vec_t<double> g=itsBS->Gradient(r);
        rvec3_t out(0.0, 0.0, 0.0);
        for (size_t a=0;a<itsC.size();a++) out+=itsC[a]*g[a];
        return out;
    }
    const BasisSet::Orbital_1E_IBS<double>* itsBS;
    rvec_t itsC;
};
} //anonymous namespace

struct StateView::Imp
{
    // Declaration order IS the lifetime order: the owner (the basis) outlives every field built over it, and rho outlives m's views.
    std::unique_ptr<const Calculation>               owner;
    SavedState                                       saved;
    RestartSummary                                   judged;
    std::vector<StateBlock>                          blocks;
    std::unique_ptr<const ScalarFunction<double>>    rho, m;
    std::vector<Calculation::Level>                  levels;
    //! Level i's orbital, or null with the reason in \c whyNot[i] (a complex block).
    std::vector<std::unique_ptr<const ScalarFunction<double>>> orbitals;
    std::vector<std::string>                         whyNot;
};

StateView::StateView(std::unique_ptr<const Calculation> owner, SavedState saved, RestartSummary judged, std::vector<StateBlock> blocks,
                     std::unique_ptr<const sf_t> rho, const sf_t* up, const sf_t* down)
    : itsImp(std::make_unique<Imp>())
{
    if (!owner || !rho) throw std::logic_error("StateView: a view needs the facade that owns its basis and the density it draws");
    if ((up==nullptr)!=(down==nullptr)) throw std::logic_error("StateView: one spin channel without the other");
    if (blocks.size()!=saved.blocks.size())
        throw std::logic_error("StateView: "+std::to_string(blocks.size())+" basis blocks for "+std::to_string(saved.blocks.size())
                               +" saved blocks -- the facade's judgement should have refused");
    Imp& I=*itsImp;
    I.owner =std::move(owner);
    I.saved =std::move(saved);
    I.judged=std::move(judged);
    I.blocks=std::move(blocks);
    I.rho   =std::move(rho);
    if (up) I.m=std::make_unique<ChannelDifference>(up, down);
    else    I.m=std::make_unique<ZeroSpin>();

    // THE LEVELS, built exactly as Calculation::WalkLevels builds the live list: block by block in the saved (canonical) order, each
    // block's orbitals in their saved order (WriteOrbitals walked the live block's Iterate()), then ONE stable sort on e.  The label
    // is the live orbital's: Orbital_QNs(principal offset + its 1-based index, the block's irrep), streamed.
    struct Entry { Calculation::Level level; std::unique_ptr<const ScalarFunction<double>> phi; std::string whyNot; };
    std::vector<Entry> entries;
    bool any=false;
    for (size_t b=0;b<I.blocks.size();b++)
    {
        const StateBlock& blk=I.blocks[b];
        const SavedBlock& sb =I.saved.blocks[b];
        if (sb.C.columns()==0) continue;   // a state that saved no orbitals: no levels from this block
        any=true;
        if (sb.C.rows()!=blk.n())
            throw std::logic_error("StateView: block "+std::to_string(b)+"'s saved orbitals have "+std::to_string(sb.C.rows())
                                   +" coefficients on a "+std::to_string(blk.n())+"-function block");
        rvec3_t k(0,0,0);
        if (auto it=blk.extras.find("k"); it!=blk.extras.end() && it->second.size()==3) k=rvec3_t(it->second[0], it->second[1], it->second[2]);
        const size_t offset=blk.irrep.sym ? blk.irrep.sym->GetPrincipleOffset() : 0;
        for (size_t i=0;i<sb.C.columns();i++)
        {
            std::ostringstream label;
            label<<Orbital_QNs(offset+i+1, blk.irrep);
            Entry e{.level={.e=sb.eps[i], .occupation=sb.f[i], .spin=blk.irrep.ms, .k=k, .weight=blk.weight, .label=label.str()}};
            if (blk.real)
            {
                rvec_t c(sb.C.rows());
                for (size_t a=0;a<c.size();a++) c[a]=sb.C(a,i).real();
                e.phi=std::make_unique<SavedRealOrbital>(blk.real, std::move(c));
            }
            else e.whyNot="level "+e.level.label+" is a complex Bloch orbital; drawing one needs a choice of part or modulus, not made yet";
            entries.push_back(std::move(e));
        }
    }
    std::stable_sort(entries.begin(), entries.end(), [](const Entry& a, const Entry& b){ return a.level.e < b.level.e; });
    for (Entry& e : entries)
    {
        I.levels.push_back(std::move(e.level));
        I.orbitals.push_back(std::move(e.phi));
        I.whyNot.push_back(std::move(e.whyNot));
    }

    // THE VIEW REPORTS ITSELF one level in, under the scope its caller opened (report::Section "view").
    namespace rpt = qchem::report;
    const StateSummary& s=I.saved.summary;
    rpt::Set("state",      I.saved.path);
    rpt::Set("savedBy",    s.label);
    rpt::Set("converged",  s.converged);
    rpt::Set("iterations", (long)s.iterations);
    rpt::Set("energy",     s.energy);
    rpt::Set("exact",      I.judged.exact);
    rpt::Set("warm",       I.judged.warm);
    rpt::Set("electrons",  I.judged.electrons);
    rpt::Set("levels",     (long)I.levels.size());
    if (!any) rpt::Set("orbitals", "none saved (a solid 1.x state): the density and the spin density draw, no level does");
}

StateView::~StateView() = default;

const StateView::sf_t& StateView::Density()     const {return *itsImp->rho;}
const StateView::sf_t& StateView::SpinDensity() const {return *itsImp->m;}
const std::vector<Calculation::Level>& StateView::GetLevels() const {return itsImp->levels;}
const StateSummary&   StateView::Summary()   const {return itsImp->saved.summary;}
const RestartSummary& StateView::Judgement() const {return itsImp->judged;}

Outcome<const StateView::sf_t*, RestartRefusal> StateView::Orbital(size_t i) const
{
    using O=Outcome<const sf_t*, RestartRefusal>;
    const Imp& I=*itsImp;
    if (I.levels.empty())
        return O::Fail({RestartRefusal::Why::Format, "'"+I.saved.path+"' saved no orbitals (a solid 1.x state): only the density draws"});
    if (i>=I.levels.size())
        return O::Fail({RestartRefusal::Why::Mismatch, "level "+std::to_string(i)+" asked of a state with "+std::to_string(I.levels.size())+" levels"});
    if (!I.orbitals[i]) return O::Fail({RestartRefusal::Why::Realness, I.whyNot[i]});
    return O::Ok(I.orbitals[i].get());
}

} // namespace qchem
