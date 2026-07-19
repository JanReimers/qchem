// File: Hamiltonian/Internal/Imp/PWTerms.C  Plane-wave Kohn-Sham term implementations.
module;
#include <cassert>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <stdexcept>
module qchem.Hamiltonian.Internal.PWTerms;
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.ChargeDensity.FourierDensity;   // cast cd UP to its reciprocal-space coefficients rho-tilde
import qchem.BasisSet.Band_FT_IBS;         // cast bs UP to the reciprocal-space DFT capability (Hartree/XC)
import qchem.BasisSet.G_FieldEvaluator;    // the fit basis's FFT grid engine (RhoOnGrid/Integral for XC)
import qchem.Pseudopotential.Integrals_Pseudo;   // cast bs ACROSS to the external-PP operator-assembly mixin (PW_Pseudo)
import qchem.Fitting.FunctionFitter;        // Fitting::Factory (both PW fitters) + ProjectedDensity_G / ProjectedScalar_R
import qchem.Structure;                       // Structure::isFinite()/SumFormFactors() -- the G=0 alignment (term-side)
import qchem.Blaze;                            // blazem::zeroH<dcmplx> (the null-PP V_long block)

namespace qchem::Hamiltonian
{

bool& ReportGridCharge() { static bool on = false; return on; }

PW_Pseudo::PW_Pseudo(const st_t& st, const Pseudopotential::LocalPotential* loc,
                         const Pseudopotential::SeparablePotential* nl)
    : cStatic_HT_Imp()
    , theStructure(st)
    , itsLocal(loc)
    , itsSep(nl)
{
    assert(st->GetNumAtoms()>0);
    assert(loc && "PW_Pseudo: the term owns the local pseudopotential model (must be non-null)");
}

// Assemble the external matrix from the MODELS the term owns: hand the basis the abstract local +
// optional KB nonlocal models and let it assemble <i|V_ext|j>.  The dynamic_cast is the sanctioned
// abstract->abstract move (cobs_t = Orbital_1E_IBS<dcmplx> ACROSS to the Integrals_Pseudo capability); only a
// basis that supports reciprocal-space PP assembly answers it.  V = V_loc + V_NL.
chmat_t PW_Pseudo::CalculateMatrix(const cobs_t* bs, const Spin&) const
{
    auto pw=dynamic_cast<const Pseudopotential::Integrals_Pseudo<dcmplx>*>(bs);
    assert(pw && "PW_Pseudo requires an Integrals_Pseudo<dcmplx> (e.g. plane-wave) basis");
    // SHORT-range local only: the LONG (softened-Coulomb) part folds into the Hartree Poisson via PW_Hartree
    // (the CP2K local-PP split, doc/GPWPlan.md 0e-PP).  V_ext(short) = V_loc(short) + V_NL.
    chmat_t V=pw->MakeLocalPotentialShort(&*theStructure, *itsLocal);
    if (itsSep) V += pw->MakeSeparablePotential(&*theStructure, *itsSep);
    return V;
}

void PW_Pseudo::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    // Een stays the band expectation over the (G!=0) external matrix (== the prototype's electron-ion
    // energy).  The dropped-G=0 alignment alpha is a separate constant in te.Ealign -- kept in the total
    // energy but NOT the matrix, and out of the band-structure cross-check.  It is E_alpha = (N/Omega)
    // Sum_a alpha_a with alpha_a = the model's finite G->0 limit (FormFactorG0 = integral[V_loc^a+Z/r]).
    // The alignment is a PERIODIC neutralising-background artifact: it exists only when the Structure is
    // periodic (!isFinite(), Omega finite); a finite/molecular Structure has no G=0 background, so no
    // alignment term.  The term owns the model (alpha) and reads Omega straight off the geometry -- no basis
    // (the old basis-side PseudoG0Energy was a scalar PP formula leaking into the neutral basis interface).
    te.Een=cd->DM_Contract(this);                                          // integral rho V_ext(short) (G!=0)
    // A finite/molecular structure has no G=0 neutralising background, so NO alignment (even though its atoms
    // DO have form factors -- the physics decision lives here, via isFinite(), not in the geometry).  For a
    // periodic cell the structure folds in 1/Omega: SumFormFactors returns (1/Omega) Sum_a alpha_a.  The SHORT
    // part's alignment only -- the LONG part's G=0 shift moves to PW_Hartree with V_long (doc/GPWPlan.md 0e-PP).
    if (!theStructure->isFinite())                                         // periodic only (Omega finite)
        te.Ealign = cd->GetTotalCharge() *
                    theStructure->SumFormFactors([this](int Z){return itsLocal->FormFactorG0Short(Z);});
}

std::ostream& PW_Pseudo::Write(std::ostream& os) const
{
    return os << "    PW external (pseudo)potential with " << theStructure->GetNumAtoms() << " atoms." << std::endl;
}

//----------------------------------------------------------------------------------- Kinetic
chmat_t PW_Kinetic::CalculateMatrix(const cobs_t* bs, const Spin&) const
{
    chmat_t T=bs->MakeKinetic();   // <p^2> block (no 1/2)
    T*=0.5;                        // T = 1/2 <p^2>
    return T;
}

void PW_Kinetic::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    te.Kinetic=cd->DM_Contract(this);   // <T> = integral rho (1/2 p^2)
}

std::ostream& PW_Kinetic::Write(std::ostream& os) const
{
    return os << "    PW kinetic energy 1/2<p^2>." << std::endl;
}

// Ion-ion (Ewald) is now the shared IonIon<dcmplx> term (qchem.Hamiltonian.Internal.IonIon).

//----------------------------------------------------------------------------------- Hartree
// Holds its CD fit basis (from Ham_PW_DFT::BuildTerms via the orbital basis's factory, never assuming
// orbital==fit) and hands it to the density's GetRepulsion3C each SCF cycle -- mirrors FittedVee holding its
// CD fit basis and calling IrrepCD::GetRepulsion3C(fbs).
PW_Hartree::PW_Hartree(fbs_t fb, st_t st, const Pseudopotential::LocalPotential* loc)
    : itsFitBasis(fb)
    , theStructure(st)
    , itsLocal(loc)
{}

// The fixed <i|V_long|j> block for this orbital basis: the long-range (softened-Coulomb) local-PP matrix,
// built ONCE per irrep basis (density-independent) and cached by BasisSetID.  Null model => zero (pure
// Hartree, no core charge to fold).  Assembled through the SAME Integrals_Pseudo cross-cast PW_Pseudo uses.
const chmat_t& PW_Hartree::LongBlock(const cobs_t* bs) const
{
    const std::string id=bs->BasisSetID();
    auto it=itsLongBlocks.find(id);
    if (it!=itsLongBlocks.end()) return it->second;
    chmat_t VL;
    if (itsLocal)
    {
        auto pp=dynamic_cast<const Pseudopotential::Integrals_Pseudo<dcmplx>*>(bs);
        assert(pp && "PW_Hartree: the V_long fold needs an Integrals_Pseudo<dcmplx> orbital basis");
        VL=pp->MakeLocalPotentialLong(&*theStructure, *itsLocal);
    }
    else
        VL=blazem::zeroH<dcmplx>(bs->GetNumFunctions());
    return itsLongBlocks.emplace(id, std::move(VL)).first->second;
}

chmat_t PW_Hartree::CalcMatrix(const cobs_t* bs, const Spin&, const cChargeDensity* cd) const
{
    newCD(cd);   // dirty the Irrep cache if cd is new (the cross-iteration freshness mechanism)
    auto fd=dynamic_cast<const qchem::ChargeDensity::FourierDensity*>(cd);
    assert(fd && "PW_Hartree requires a FourierDensity (periodic) charge density");
    auto bft=dynamic_cast<const BasisSet::Band_FT_IBS*>(bs);
    assert(bft && "PW_Hartree requires a Band_FT_IBS (reciprocal-space DFT) orbital basis");
    // The density contracts D against the basis's D-free Coulomb tensor Repulsion3C (kernel baked) to give
    // V_H(dm); the term assembles <i|V_H|j> = V_H(G_i-G_j) via the orbital basis's MakeOverlap bridge.  D
    // never crosses into the basis.  Then ADD the fixed long-range core-charge field V_long (the CP2K
    // local-PP split): the total electrostatic field is V_H[rho_elec] + V_long -- one Poisson-solved matrix.
    ΔG_Map VH=fd->GetRepulsion3C(*itsFitBasis);
    chmat_t H=bft->MakeOverlap([&VH](const ivec3_t& dm)->dcmplx
        { auto it=VH.find(dm); return it==VH.end()?dcmplx(0.0):it->second; });
    H+=LongBlock(bs);
    return H;
}

void PW_Hartree::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    newCD(cd);
    // The Fock matrix carries V_H + V_long, so DM_Contract(this) = Tr(D(V_H+V_long)).  Split into the two
    // physical energies (doc/GPWPlan.md 0e-PP): the electron-electron Hartree gets the 1/2 double-counting
    // factor, the electron-ion long-range gets none.  V_long is the FIXED per-basis block (DM_ContractBlocks).
    // (this->CalcMatrix has run for every irrep by now -- DM_Contract triggers GetMatrix -- so itsLongBlocks
    // is fully populated.)
    double total=cd->DM_Contract(this,cd);                 // Tr(D (V_H + V_long))
    double eLong=cd->DM_ContractBlocks(itsLongBlocks);     // Tr(D V_long) = E_een,long
    te.Eee += 0.5*(total-eLong);                           // E_H = 1/2 integral rho V_H[rho]
    te.Een += eLong;                                       // electron-ion long-range (no 1/2)
    // The dropped-G=0 alignment of the LONG part moves here with V_long (the SHORT part's stays in PW_Pseudo).
    if (itsLocal && !theStructure->isFinite())             // periodic only (Omega finite)
        te.Ealign += cd->GetTotalCharge() *
                     theStructure->SumFormFactors([this](int Z){return itsLocal->FormFactorG0Long(Z);});
}

std::ostream& PW_Hartree::Write(std::ostream& os) const
{
    return os << "    PW electrostatics: Hartree V_H[rho] + long-range core-charge V_long (G-space)." << std::endl;
}

//----------------------------------------------------------------------------------- XC
namespace
{
// v_xc = functional(rho) on the fit basis's OWN grid, presented as a ProjectedScalar_R the ortho scalar fitter
// samples in bulk (this is what fixed the old flaw of quadraturing on the ORBITAL basis's grid).  rho(r) is
// PRECOMPUTED by the term (one inverse FFT, shared with GetEnergy); this field just maps the functional over
// those grid values -- the functional stays Hamiltonian-side (no qcBasisSet->qcHamiltonian library cycle).
// It is GRID-BOUND: only the ortho fitter samples it, in bulk, on exactly the grid rho was transformed onto.
class PWVxcField
    : public virtual ScalarFunction<double>
    , public         Fitting::ProjectedScalar_R
{
public:
    PWVxcField(const ExFunctional* xc, const rvec_t& rhoGrid, const BasisSet::G_FieldEvaluator* grid)
        : itsXc(xc), itsRhoGrid(rhoGrid), itsGrid(grid) {}

    // Pointwise is NOT supported: this field carries only grid values, and nothing samples it pointwise (the
    // ortho fitter uses the bulk overload).  Make the grid-bound contract explicit rather than silently wrong.
    virtual double  operator()(const rvec3_t&) const override
        {throw std::logic_error("PWVxcField is grid-bound: sample it in bulk on the fit grid, not pointwise");}
    virtual rvec3_t Gradient  (const rvec3_t&) const override {return rvec3_t(0,0,0);}

    // Bulk: v_xc = functional(rho) at the precomputed grid values.  The points MUST be the fit grid rho was
    // transformed onto -- assert IDENTITY (not merely size), so a future diagnostic that samples a different
    // same-cardinality point set fails loudly instead of pairing values to the wrong points (review #2).
    virtual rvec_t  operator()(const rvec3vec_t& rs) const override
    {
        assert(SampledOnGrid(rs) && "PWVxcField: must be sampled on the fit basis's own grid (identity)");
        rvec_t v(itsRhoGrid.size());
        for (size_t q=0;q<v.size();q++) v[q]=itsXc->GetVxc(itsRhoGrid[q]);
        return v;
    }

    virtual const ScalarFunction<double>* GetScalarFunction() const override {return this;}
private:
    bool SampledOnGrid(const rvec3vec_t& rs) const   // debug-only identity check (GridPoints() is cached)
    {
        const rvec3vec_t& g=itsGrid->GridPoints();
        if (rs.size()!=g.size()) return false;
        for (size_t q=0;q<rs.size();q++)
        {
            double dx=rs[q].x-g[q].x, dy=rs[q].y-g[q].y, dz=rs[q].z-g[q].z;
            if (dx*dx+dy*dy+dz*dz > 1e-24) return false;
        }
        return true;
    }
    const ExFunctional*               itsXc;
    const rvec_t&                     itsRhoGrid;   // precomputed rho(r) on the fit grid (owned by PW_XC; field is transient)
    const BasisSet::G_FieldEvaluator* itsGrid;
};
} // anonymous

// Built once (in Ham_PW_DFT::BuildTerms) with its Vxc fit basis -- the overlap-metric sibling of PW_Hartree.
// The XC QUADRATURE GRID comes from the FIT basis (not the orbital), so relCutoff / the functional's
// GridCutoffFactor control the Vxc/E_xc grid.  The fitter OWNS that grid; this term borrows it via
// itsScalarFitter->Grid() -- one owner, no second cross-cast of the fit basis (#7).
PW_XC::PW_XC(const xc_t& xc, fbs_t fb)
    : itsXc(xc)
    , itsVxcFitBasis(fb)                       // hand it to the density's GetFourierDensity (its Overlap3C key)
    , itsScalarFitter(Fitting::Factory(fb))   // the ortho (G-space) scalar fitter -- owns the FFT quadrature grid
{}
PW_XC::~PW_XC() = default;   // itsScalarFitter's abstract type is complete here

// rho(r) on the fit grid for cd -- one inverse FFT, recomputed only on a new density serial (newCD), so
// CalcMatrix and GetEnergy share it (whichever runs first this iteration pays; the other reuses).
void PW_XC::RefreshRhoGrid(const cChargeDensity* cd) const
{
    if (!newCD(cd)) return;
    auto fd=dynamic_cast<const qchem::ChargeDensity::FourierDensity*>(cd);
    assert(fd && "PW_XC requires a FourierDensity (periodic) charge density");
    itsRhoGrid=itsScalarFitter->Grid().RhoOnGrid(fd->GetFourierDensity(*itsVxcFitBasis));   // rho-tilde via Overlap3C, onto the FIT grid
    if (ReportGridCharge())
    {
        // Grid charge vs analytic charge: the electrons LOST to grid truncation (high-G aliasing of rho).
        // == CP2K's "Electronic density on regular grids: <int rho> <error>" -- a controlled cutoff metric.
        const double qGrid=itsScalarFitter->Grid().Integral(itsRhoGrid);   // integral rho_grid d3r
        const double qDM  =cd->GetTotalCharge();                           // Tr(D S) (analytic, ~ N)
        std::cout << "[grid charge] integral rho_grid=" << std::fixed << std::setprecision(6) << qGrid
                  << "  Tr(DS)=" << qDM
                  << "  lost=" << std::scientific << std::setprecision(3) << (qGrid-qDM)
                  << std::defaultfloat << std::endl;
    }
}

// XC through the pre-built ortho scalar fitter, mirroring the molecular FittedVxc: the fitter batch-samples
// the v_xc(rho) field on the FIT basis's grid and forward-transforms it (the projection IS the fit on the
// orthonormal {G}); the ORBITAL basis then assembles <i|v_xc|j>.  No O(Npts*n^2) pointwise density sampling.
chmat_t PW_XC::CalcMatrix(const cobs_t* bs, const Spin&, const cChargeDensity* cd) const
{
    RefreshRhoGrid(cd);
    itsScalarFitter->DoFit(PWVxcField(itsXc.get(), itsRhoGrid, &itsScalarFitter->Grid()));
    return itsScalarFitter->Overlap(bs);                                        // <i|v_xc|j> (no kernel)
}

void PW_XC::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    RefreshRhoGrid(cd);   // reuses CalcMatrix's transform this iteration (same density serial)
    rvec_t exc(itsRhoGrid.size());
    for (size_t q=0;q<itsRhoGrid.size();q++) {double ro=itsRhoGrid[q]; exc[q]=itsXc->GetEpsXc(ro)*ro;}
    te.Exc += itsScalarFitter->Grid().Integral(exc);   // E_xc = integral eps_xc(rho) rho, on the fitter's grid
}

std::ostream& PW_XC::Write(std::ostream& os) const
{
    return os << "    PW exchange-correlation potential v_xc(rho(r))." << std::endl;
}

} //namespace
