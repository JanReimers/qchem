// File: Hamiltonian/Internal/Imp/PWTerms.C  Plane-wave Kohn-Sham term implementations.
module;
#include <algorithm>   // std::min (the threaded quadrature's output-column blocking)
#include <cassert>
#include <complex>
#include <cstdlib>
#include <exception>   // std::exception_ptr (throw containment across the threaded Phi build)
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
import qchem.Pseudopotential.Integrals_Pseudo;   // cast bs ACROSS to the external-PP operator-assembly mixin (Ven_PP_*)
import qchem.Fitting.FunctionFitter;        // Fitting::Factory (both PW fitters) + ProjectedDensity_G / ProjectedScalar_R
import qchem.Structure;                       // Structure::isFinite()/SumFormFactors() -- the G=0 alignment (term-side)
import qchem.Blaze;                            // blazem::zeroH<dcmplx> (the null-PP V_long block)
import qchem.Mesh.Quadrature;                 // qcMesh::Mesh (the DeltaFittedVxc engine's quadrature mesh)
import qchem.Reporting;                       // Timed (the setup/scf timing ledger)
import qchem.Parallel;                         // WorkerThreads (GPW_OMP_THREADS -- the XC-mesh table + quadrature loops)

namespace qchem::Hamiltonian
{

bool& ReportGridCharge() { static bool on = false; return on; }

// The dropped-G=0 alignment, PER ELECTRON, evaluated ONCE at construction (R2.16).
//
// E_alpha = N * (1/Omega) Sum_a alpha_a, with alpha_a the model's finite G->0 limit
// (FormFactorG0 = integral[V_loc^a + Z/r]).  It is kept in the total energy but NOT in the matrix, so it
// stays out of the band-structure cross-check.  The alignment is a PERIODIC neutralising-background
// artifact: a finite/molecular Structure has no G=0 background, so its coefficient is exactly ZERO (even
// though its atoms DO have form factors -- the physics decision lives here, not in the geometry).
//
// Both the isFinite() question and the SumFormFactors sum are answered HERE rather than per energy call:
// the structure and the model are fixed at construction, so the answer cannot change during a run.  Each
// GetEnergy then just scales by the current electron count.
namespace
{
double G0AlignmentPerElectron(const Structure& st, const std::function<double(int)>& formFactorG0)
{
    if (st.isFinite()) return 0.0;              // no neutralising background => no alignment, exactly 0
    return st.SumFormFactors(formFactorG0);     // periodic: folds in 1/Omega
}
} // namespace

Ven_PP_Short::Ven_PP_Short(const st_t& st, const Pseudopotential::LocalPotential* loc)
    : cStatic_HT_Imp()
    , theStructure(st)
    , itsLocal(loc)
{
    assert(st->GetNumAtoms()>0);
    assert(loc && "Ven_PP_Short: the term owns the local pseudopotential model (must be non-null)");
    itsAlphaZ=G0AlignmentPerElectron(*st, [loc](int Z){return loc->FormFactorG0Short(Z);});
}

// Assemble the external matrix from the MODEL the term owns: hand the basis the abstract local model and
// let it assemble <i|V_loc,short|j>.  The dynamic_cast is the sanctioned abstract->abstract move
// (cobs_t = Orbital_1E_IBS<dcmplx> ACROSS to the Integrals_Pseudo capability); only a basis that supports
// reciprocal-space PP assembly answers it.
chmat_t Ven_PP_Short::MakeMatrix(const cobs_t* bs, const Spin&) const
{
    auto pw=dynamic_cast<const Pseudopotential::Integrals_Pseudo<dcmplx>*>(bs);
    assert(pw && "Ven_PP_Short requires an Integrals_Pseudo<dcmplx> (e.g. plane-wave) basis");
    // SHORT-range local only.  The LONG (softened-Coulomb) half is Ven_PP_Long and the KB projectors are
    // Ven_PP_NonLocal (the CP2K local-PP split, doc/GPWPlan.md 0e-PP).
    return pw->MakeLocalPotentialShort(&*theStructure, *itsLocal);
}

void Ven_PP_Short::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    // Een is the band expectation over the (G!=0) matrix (== the prototype's electron-ion energy).
    te.Een     += cd->DM_Contract(this);                 // integral rho V_loc,short (G!=0)
    te.E_alphaZ+= cd->GetTotalCharge()*itsAlphaZ;        // SHORT G=0 alignment (0 for a finite structure)
}

std::ostream& Ven_PP_Short::Write(std::ostream& os) const
{
    return os << "    PW electron-ion: SHORT-range local PP, "
              << theStructure->GetNumAtoms() << " atoms." << std::endl;
}

//--------------------------------------------------------------- Ven_PP_NonLocal (the KB projectors)
Ven_PP_NonLocal::Ven_PP_NonLocal(const st_t& st, const Pseudopotential::SeparablePotential* nl)
    : cStatic_HT_Imp()
    , theStructure(st)
    , itsSep(nl)
{
    assert(st && st->GetNumAtoms()>0);
    // REQUIRED, not optional: a local-only PP omits this TERM rather than adding an all-zeros one.
    if (!nl)
        throw std::runtime_error("Ven_PP_NonLocal: the KB projector term requires a SeparablePotential "
                                 "model.  A local-only pseudopotential must not add this term.");
}

chmat_t Ven_PP_NonLocal::MakeMatrix(const cobs_t* bs, const Spin&) const
{
    auto pw=dynamic_cast<const Pseudopotential::Integrals_Pseudo<dcmplx>*>(bs);
    assert(pw && "Ven_PP_NonLocal requires an Integrals_Pseudo<dcmplx> (e.g. plane-wave) basis");
    if (std::getenv("GPW_NL_PER_L"))
    {   // I0 diagnostic (doc/SphericalLatticePlan.md): bank the per-l blocks once per irrep block.
        const std::string id=bs->BasisSetID();
        if (itsByLSeen.insert(id).second)
            for (auto& lH : pw->MakeSeparablePotentialByL(&*theStructure, *itsSep))
                itsByL[lH.first].emplace(id, std::move(lH.second));
    }
    return pw->MakeSeparablePotential(&*theStructure, *itsSep);
}

void Ven_PP_NonLocal::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    // Electron-ion, and short-ranged by construction, so no G=0 alignment of its own.
    const double eNL = cd->DM_Contract(this);            // Tr(D V_NL)
    te.Een   += eNL;
    te.EenNL += eNL;                                     // the diagnostic V_loc/V_NL split (Een keeps the total)
    if (!itsByL.empty())
    {   // GPW_NL_PER_L: the per-channel decomposition of eNL (l=-1 = a basis that answered lumped).
        std::cout << "[NL per-l]";
        double sum=0;
        for (const auto& lB : itsByL)
        {
            const double el=cd->DM_ContractBlocks(lB.second);
            sum+=el;
            std::cout << "  l=" << lB.first << ": " << el;
        }
        std::cout << "  (sum=" << sum << "  Tr(D V_NL)=" << eNL << ")" << std::endl;
    }
}

std::ostream& Ven_PP_NonLocal::Write(std::ostream& os) const
{
    return os << "    PW electron-ion: KB separable nonlocal projectors, "
              << theStructure->GetNumAtoms() << " atoms." << std::endl;
}

// Kinetic is now the shared Kinetic<dcmplx> term (qchem.Hamiltonian.Internal.Kinetic).
// Ion-ion (Ewald) is now the shared IonIon<dcmplx> term (qchem.Hamiltonian.Internal.IonIon).

//------------------------------------------------------------------- Ven_PP_Long (the LONG range half)
// The long-range (softened-Coulomb / Gaussian core-charge) local-PP matrix.  DENSITY-INDEPENDENT --
// MakeLocalPotentialLong takes only (structure, model) -- so this is an ordinary static term and rides the
// standard per-Irrep static cache, exactly as Ven_PP_Short does.  It was previously a side-block inside the
// Hartree term (summed into that term's matrix, then subtracted back out of its energy) purely because the
// two are solved through the same G-space Poisson machinery; that is a COMPUTATIONAL kinship, not a
// physical one, and it cost a nullable model, a second block cache, and a "Hartree" term that contributed
// to E_een.  Assembled through the SAME Integrals_Pseudo cross-cast Ven_PP_Short uses.
Ven_PP_Long::Ven_PP_Long(const st_t& st, const Pseudopotential::LocalPotential* loc)
    : cStatic_HT_Imp()
    , theStructure(st)
    , itsLocal(loc)
{
    assert(st && st->GetNumAtoms()>0);
    // REQUIRED, not optional: a run with no local PP omits this TERM rather than constructing an
    // all-zeros one.  (Absence of a capability belongs in the term LIST, not in a per-call branch.)
    if (!loc)
        throw std::runtime_error("Ven_PP_Long: the long-range local-PP term requires a LocalPotential "
                                 "model.  A run without a local pseudopotential must not add this term.");
    itsAlphaZ=G0AlignmentPerElectron(*st, [loc](int Z){return loc->FormFactorG0Long(Z);});
}

chmat_t Ven_PP_Long::MakeMatrix(const cobs_t* bs, const Spin&) const
{
    auto pp=dynamic_cast<const Pseudopotential::Integrals_Pseudo<dcmplx>*>(bs);
    assert(pp && "Ven_PP_Long requires an Integrals_Pseudo<dcmplx> (e.g. plane-wave) basis");
    return pp->MakeLocalPotentialLong(&*theStructure, *itsLocal);
}

void Ven_PP_Long::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    // Electron-ION, so NO 1/2: the double-counting factor belongs to the electron-electron Hartree alone.
    te.Een     += cd->DM_Contract(this);                 // Tr(D V_long) = E_een,long
    te.E_alphaZ+= cd->GetTotalCharge()*itsAlphaZ;        // LONG G=0 alignment (0 for a finite structure)
}

std::ostream& Ven_PP_Long::Write(std::ostream& os) const
{
    return os << "    PW electron-ion: LONG-range local PP (Gaussian core charge, G-space), "
              << theStructure->GetNumAtoms() << " atoms." << std::endl;
}

//----------------------------------------------------------------------------------- Hartree
// Holds its CD fit basis (from Ham_PW_DFT::BuildTerms via the orbital basis's factory, never assuming
// orbital==fit) and hands it to the density's GetRepulsion3C each SCF cycle -- mirrors FittedVee holding its
// CD fit basis and calling IrrepCD::GetRepulsion3C(fbs).  Pure V_H[rho_elec]: no structure, no PP model.
Vee_Hartree::Vee_Hartree(fbs_t fb)
    : itsFitBasis(fb)
{}

chmat_t Vee_Hartree::MakeMatrix(const cobs_t* bs, const Spin&, const cChargeDensity* cd) const
{
    newCD(cd);   // dirty the Irrep cache if cd is new (the cross-iteration freshness mechanism)
    auto fd=dynamic_cast<const qchem::ChargeDensity::FourierDensity*>(cd);
    assert(fd && "Vee_Hartree requires a FourierDensity (periodic) charge density");
    auto bft=dynamic_cast<const BasisSet::Band_FT_IBS*>(bs);
    assert(bft && "Vee_Hartree requires a Band_FT_IBS (reciprocal-space DFT) orbital basis");
    // The density contracts D against the basis's D-free Coulomb tensor Repulsion3C (kernel baked) to give
    // V_H(dm) [FORWARD]; the KS matrix <i|V_H|j> = Σ_k V_H(G_k) <i|e^{iG_k}|j> is the BACKWARD contraction of the
    // SAME Repulsion3C tensor over the CD fit basis (its applyAdjoint -- the overlap integrate-back on the fit
    // grid; the Coulomb kernel is forward-only, already in V_H).  So forward AND backward run on the one fit grid
    // (doc/GPWPlan §0e step 2).
    ΔG_Map VH=fd->GetRepulsion3C(*itsFitBasis);
    return ContractAdjointG_ERI3(bft->Repulsion3C(*itsFitBasis),
        [&VH](const ivec3_t& dm)->dcmplx { auto it=VH.find(dm); return it==VH.end()?dcmplx(0.0):it->second; });
}

void Vee_Hartree::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    newCD(cd);
    // E_H = 1/2 integral rho V_H[rho] -- the 1/2 is the electron-electron double-counting factor.  The Fock
    // block is now V_H ALONE, so this is a direct contraction: the old form computed Tr(D(V_H+V_long)) and
    // subtracted Tr(D V_long) back off, because V_long was folded into the same matrix.
    te.Eee += 0.5*cd->DM_Contract(this,cd);
}

std::ostream& Vee_Hartree::Write(std::ostream& os) const
{
    return os << "    PW electron-electron: Hartree V_H[rho] (G-space Poisson)." << std::endl;
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
    const rvec_t&                     itsRhoGrid;   // precomputed rho(r) on the fit grid (owned by PWFittedVxc; field is transient)
    const BasisSet::G_FieldEvaluator* itsGrid;
};
} // anonymous

// Built once (in Ham_PW_DFT::BuildTerms) with its Vxc fit basis -- the overlap-metric sibling of Vee_Hartree.
// The XC QUADRATURE GRID comes from the FIT basis (not the orbital), so relCutoff / the functional's
// GridCutoffFactor control the Vxc/E_xc grid.  The fitter OWNS that grid; this term borrows it via
// itsScalarFitter->Grid() -- one owner, no second cross-cast of the fit basis (#7).
PWFittedVxc::PWFittedVxc(const xc_t& xc, fbs_t fb)
    : itsXc(xc)
    , itsVxcFitBasis(fb)                       // hand it to the density's GetFourierDensity (its Overlap3C key)
    , itsScalarFitter(Fitting::Factory(fb))   // the ortho (G-space) scalar fitter -- owns the FFT quadrature grid
{
    // Announce the quadrature grid ONCE per distinct fit basis (the exchange+correlation pair shares one):
    // the uniform-route sibling of the [Becke grid] line, so the console/report always carries the XC
    // grid's N + spacing.  (UnitCell::CreateIntegrationMesh cannot report this grid -- the XC raster is
    // the FFT engine's own 5-smooth grid, never built through the Structure mesh factory.)
    static const void* lastAnnounced=nullptr;
    if (itsVxcFitBasis.get()!=lastAnnounced)
    {
        lastAnnounced=itsVxcFitBasis.get();
        itsScalarFitter->Grid().EmitGridReport();
    }
}
PWFittedVxc::~PWFittedVxc() = default;   // itsScalarFitter's abstract type is complete here

// rho(r) on the fit grid for cd -- one inverse FFT, recomputed only on a new density serial (newCD), so
// MakeMatrix and GetEnergy share it (whichever runs first this iteration pays; the other reuses).
void PWFittedVxc::RefreshRhoGrid(const cChargeDensity* cd) const
{
    if (!newCD(cd)) return;
    auto fd=dynamic_cast<const qchem::ChargeDensity::FourierDensity*>(cd);
    assert(fd && "PWFittedVxc requires a FourierDensity (periodic) charge density");
    // RAW-FIRST (doc/GPWPlan 0.5(f2)): a collocation-backed density answers with rho_DM(r) directly --
    // pointwise >= 0 for an aufbau D, so the rho>0 guard never bites and the grid calibration can relax
    // (the C=8 driver was the BALL path's Gibbs lobes).  A plane-wave density (or a matrix-free seed)
    // answers EMPTY and we take the ball round trip below -- bit-identical to the pre-f2 path.
    itsRhoGrid=fd->GetRhoOnGrid(*itsVxcFitBasis);
    itsRhoIsRaw=(itsRhoGrid.size()!=0);
    // ROUTE STABILITY (R2.16).  RAW and BALL are not two ways to compute one number: BALL's ball-projected
    // rho is non-variational (H_xc != dE_xc/dD), RAW's is exact, so switching route mid-SCF switches the
    // FUNCTIONAL BEING MINIMISED underneath the optimiser.  The capability is a property of the ORBITAL
    // BASIS lineage -- GPW's Overlap3CTensor always sets applyRaw, a plane-wave basis never does -- i.e. it
    // is fixed at construction, and so the route must be too.
    //
    // The ONE unavoidable exception is the SEED: a matrix-free density (iteration 0) has no D to collocate,
    // so rho_DM = phi^T D phi does not exist and it can only answer BALL.  That is inherent to seeding, not
    // a design choice, and iteration 0's energy is discarded anyway.  So: exempt the matrix-free seed,
    // LATCH the route on the first density-matrix-backed density, and THROW if it ever changes after that.
    // A silent mid-SCF functional change becomes a loud error.
    if (dynamic_cast<const ChargeDensity::cDM_CD*>(cd))    // matrix-backed: the route must be stable from here
    {
        if (!itsRouteLatched) {itsRouteLatched=true; itsLatchedRaw=itsRhoIsRaw;}
        else if (itsLatchedRaw!=itsRhoIsRaw)
            throw std::runtime_error(
                std::string("PWFittedVxc: the XC route changed mid-SCF (")
                + (itsLatchedRaw?"RAW -> BALL":"BALL -> RAW")
                + ").  These minimise DIFFERENT functionals -- BALL's ball-projected rho is non-variational "
                  "-- so the optimiser would be chasing a moving target.  The route is a property of the "
                  "orbital basis and must not change once the SCF is running.");
    }
    if (!itsRhoIsRaw)
        itsRhoGrid=itsScalarFitter->Grid().RhoOnGrid(fd->GetFourierDensity(*itsVxcFitBasis));   // rho-tilde via Overlap3C, onto the FIT grid
    // DIAGNOSTIC (env GPW_XCROUTE): which V_xc route fires this iteration -- RAW (applyRawAdjoint, FD-exact/
    // variational; gate GPW.RawXCConsistencyFD) vs BALL (DoFit/Overlap, non-variational under BallOnly;
    // gate GPW.XCPotentialConsistencyFD).  Answers whether GDM is fighting the ball non-variationality.
    if (std::getenv("GPW_XCROUTE"))
    {
        double rmin=1e300, rmax=-1e300;
        for (double r : itsRhoGrid) { rmin=std::min(rmin,r); rmax=std::max(rmax,r); }
        std::cout << "[xc route] " << (itsRhoIsRaw ? "RAW (applyRawAdjoint, variational)"
                                                   : "BALL (DoFit/Overlap, NON-variational under BallOnly)")
                  << "  rho grid=[" << std::scientific << std::setprecision(2) << rmin << ", " << rmax << "]"
                  << " npts=" << itsRhoGrid.size() << std::defaultfloat << std::endl;
    }
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
        // XC-collapse diagnostic (doc/GPWPlan §0e step 2): the collocated rho's min/max/negative content, and
        // the XC the rho>0 guard ZEROES.  If rho rings locally-negative near the sharp F peaks, GetEpsXc(rho<=0)
        // = 0 silently drops that eps_xc*rho -- the +7 Ha Exc collapse (fine -5.09 vs coarse -12.19).  Separates
        // lead (c) [Gibbs ringing + guard: rho_min very negative, big negCharge] from lead (b) [genuinely lower
        // peaks: rho_max small, rho_min ~ 0].
        double rmin=1e300, rmax=-1e300; size_t nneg=0;
        rvec_t negOnly(itsRhoGrid.size()), excLost(itsRhoGrid.size());
        for (size_t q=0;q<itsRhoGrid.size();q++)
        {
            const double r=itsRhoGrid[q];
            rmin=std::min(rmin,r); rmax=std::max(rmax,r);
            negOnly[q] = r<0.0 ? r : 0.0;                                  // negative-charge density
            // eps_xc*rho the guard drops at rho<0, estimated on |rho| (the magnitude the ringing would carry)
            excLost[q] = r<0.0 ? itsXc->GetEpsXc(-r)*(-r) : 0.0;
            if (r<0.0) ++nneg;
        }
        const double negCharge=itsScalarFitter->Grid().Integral(negOnly);
        const double excLostI =itsScalarFitter->Grid().Integral(excLost);
        std::cout << "[xc grid] rho_min=" << std::scientific << std::setprecision(3) << rmin
                  << " rho_max=" << rmax << " neg-frac=" << double(nneg)/double(itsRhoGrid.size())
                  << " negCharge=" << negCharge << " Exc_lost@rho<0~=" << excLostI
                  << std::defaultfloat << std::endl;
    }
}

// XC through the pre-built ortho scalar fitter, mirroring the molecular FittedVxc: the fitter batch-samples
// the v_xc(rho) field on the FIT basis's grid and forward-transforms it (the projection IS the fit on the
// orthonormal {G}); the ORBITAL basis then assembles <i|v_xc|j>.  No O(Npts*n^2) pointwise density sampling.
chmat_t PWFittedVxc::MakeMatrix(const cobs_t* bs, const Spin&, const cChargeDensity* cd) const
{
    RefreshRhoGrid(cd);
    if (itsRhoIsRaw)
    {
        // RAW route (0.5(f2)): H_xc through the raw adjoint of the SAME tensor whose applyRaw produced
        // itsRhoGrid -- box-truncation per level + the analytic gather -- so H_xc == dE_xc/dD of the raw
        // discrete functional to machine precision (gate: GPW.RawXCConsistencyFD).  No ball fit anywhere.
        const auto& orb=dynamic_cast<const BasisSet::Band_FT_IBS&>(*bs);   // genuine "is it?" cross-cast (throws)
        const G_ERI3& g=orb.Overlap3C(*itsVxcFitBasis);
        assert(g.applyRawAdjoint && "raw rho without a raw adjoint: Overlap3C must carry both");
        rvec_t v(itsRhoGrid.size());
        for (size_t q=0;q<itsRhoGrid.size();q++) v[q]=itsXc->GetVxc(itsRhoGrid[q]);
        return g.applyRawAdjoint(v);
    }
    itsScalarFitter->DoFit(PWVxcField(itsXc.get(), itsRhoGrid, &itsScalarFitter->Grid()));
    return itsScalarFitter->Overlap(bs);                                        // <i|v_xc|j> (no kernel)
}

void PWFittedVxc::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    RefreshRhoGrid(cd);   // reuses MakeMatrix's transform this iteration (same density serial)
    rvec_t exc(itsRhoGrid.size());
    for (size_t q=0;q<itsRhoGrid.size();q++) {double ro=itsRhoGrid[q]; exc[q]=itsXc->GetEpsXc(ro)*ro;}
    te.Exc += itsScalarFitter->Grid().Integral(exc);   // E_xc = integral eps_xc(rho) rho, on the fitter's grid
    // GPW health DIAGNOSTIC (item 2 increment 2): the signed grid-charge leak integral rho_grid - Tr(DS), the
    // electrons lost to collocation-grid truncation.  Surfaced as the solid SCF display's ρ_lost/N column.
    te.GridChargeLost = itsScalarFitter->Grid().Integral(itsRhoGrid) - cd->GetTotalCharge();
}

std::ostream& PWFittedVxc::Write(std::ostream& os) const
{
    return os << "    PW exchange-correlation potential v_xc(rho(r))." << std::endl;
}

// ---- XC_GridEngine: the pair-shared mesh + Phi tables + per-serial rho ----------------------------------

XC_GridEngine::XC_GridEngine(mesh_t mesh, Symmetry::Lattice_3D::Fold fold,
                             std::vector<Symmetry::SpinAction> sigmas, std::vector<char> flipFixed)
    : itsMesh(std::move(mesh))
    , itsFold(std::move(fold))
    , itsSigmas(std::move(sigmas))
    , itsFlipFixed(std::move(flipFixed))
{
    assert(itsMesh && itsMesh->size()>0);
    // A live fold must partition THIS mesh (empty = free run, no star-average).
    assert(itsFold.owner.empty() || itsFold.owner.size()==itsMesh->size());
    // Shubnikov (S3): σ tags require a live fold, and the zero flags cover the whole mesh.
    assert(itsSigmas.empty() || !itsFold.owner.empty());
    assert(itsFlipFixed.empty() || itsFlipFixed.size()==itsMesh->size());
}

// The (npts x n) basis table for one Bloch block: chi_i at every mesh point -- the ONE image-summed
// Bloch evaluation sweep this block ever pays (geometry-fixed, so it serves every SCF iteration).
// Keyed by the block's SPATIAL Irrep (never a pointer).  Irrep -- not BasisSetID -- because the table
// is a property of the irrep/k block within THIS run; BasisSetID's extra radial resolution exists for
// the cross-RUN disk cache, which this in-memory per-run table does not share.
const mat_t<dcmplx>& XC_GridEngine::Phi(const cobs_t* bs) const
{
    const Irrep id=bs->GetIrrep(Spin::None);
    auto it=itsPhi.find(id);
    if (it!=itsPhi.end()) return it->second;

    qchem::report::Timed timed("setup: XC-mesh Phi tables");
    const rvec3vec_t& R=itsMesh->Points();
    mat_t<dcmplx> P(R.size(), bs->GetVectorSize());
    // THE dominant SETUP bucket on an atom-centred XC mesh (MnO Γ, 48k points x 122 Bloch functions:
    // 56 s serial, and it DOUBLES on an imposed run's invariant mesh).  Each point is an independent
    // Bloch image sum writing its OWN row, so the loop parallelises with no reduction and no ordering
    // question -- the table is bit-identical at any thread count.  Opt-in (GPW_OMP_THREADS) like every
    // other parallel region here; see qchem.Parallel for why serial is the default.
    auto fill=[&](size_t g)
    {
        cvec_t phi=(*bs)(R[g]);
        for (size_t i=0; i<phi.size(); i++) P(g,i)=phi[i];
    };
#ifdef QCHEM_OPENMP
    if (const int nthreads=qchem::WorkerThreads(); nthreads>1)
    {
        std::exception_ptr firstEx;                 // throw containment (an escape = std::terminate)
        #pragma omp parallel for schedule(static) num_threads(nthreads)
        for (size_t g=0; g<R.size(); g++)
        {
            try { fill(g); }
            catch (...)
            {
                #pragma omp critical (xc_phi_throw)
                if (!firstEx) firstEx=std::current_exception();
            }
        }
        if (firstEx) std::rethrow_exception(firstEx);
        return itsPhi.emplace(id, std::move(P)).first->second;
    }
#endif
    for (size_t g=0; g<R.size(); g++) fill(g);
    return itsPhi.emplace(id, std::move(P)).first->second;
}

// rho at the mesh points, once per density serial for the WHOLE pair: the density GEMMs the cached
// tables against its private D (DM_RhoAtPoints; blocks not yet tabled self-evaluate pointwise -- first
// pass only).  A non-DM density (no DM face) falls back to the pointwise ScalarFunction sweep.
const rvec_t& XC_GridEngine::Rho(const cChargeDensity* cd, const cobs_t* ensureBlock) const
{
    assert(cd);
    // R2.9(i): the scalar and spin-resolved caches do not cross-invalidate (see the \warning on the
    // members).  An engine belongs to ONE xc/correlation pair, and a pair is either polarized or not, so
    // only one of the two routes is ever driven.  Pin it here rather than trusting the comment.
    assert(itsPolVersion==size_t(-1) && "XC_GridEngine: this engine already served RhoPol -- the scalar "
           "and spin-resolved rho caches have no cross-invalidation, so one of them would go stale");
    if (ensureBlock) Phi(ensureBlock);
    if (cd->Version()==itsRhoVersion) return itsRho;
    itsRhoVersion=cd->Version();
    qchem::report::Timed timed("scf: XC-mesh rho sampling (all iterations)");
    if (auto dm=dynamic_cast<const cDM_CD*>(cd))
        itsRho=dm->DM_RhoAtPoints(itsMesh->Points(), itsPhi);
    else
        itsRho=(*cd)(itsMesh->Points());   // non-DM (mixed rho-tilde / seed) densities batch here
    if (!itsFold.owner.empty())                    // §6a W1: the orbit-mean star-average (imposed runs only)
        Symmetry::Lattice_3D::SymmetrizeValues(itsFold, itsRho);
    return itsRho;
}

// The spin-resolved sibling of Rho: the {up,down} PAIR is cached under ONE density serial (a polarized
// density's Version() forwards to its Up child -- a single scalar cache would hand the Up raster to the
// Down channel).  A cPolarized_CD answers per channel (each channel composite GEMMs its own D against the
// SHARED Phi tables); a spin-agnostic density (the seed) collapses to rho/2 per channel, so the first
// iterations run the exact unpolarized collapse (v^sigma(rho/2,rho/2)=v^P(rho)).
const rvec_t& XC_GridEngine::RhoPol(const cChargeDensity* cd, const Spin& s, const cobs_t* ensureBlock) const
{
    assert(cd);
    assert(itsRhoVersion==size_t(-1) && "XC_GridEngine: this engine already served the scalar Rho -- the "
           "two rho caches have no cross-invalidation, so one of them would go stale");
    assert(s!=Spin::None && "XC_GridEngine::RhoPol: ask for a channel, not the total");
    if (ensureBlock) Phi(ensureBlock);
    if (cd->Version()!=itsPolVersion)
    {
        itsPolVersion=cd->Version();
        qchem::report::Timed timed("scf: XC-mesh rho sampling (all iterations)");
        if (auto pol=dynamic_cast<const ChargeDensity::cPolarized_CD*>(cd))
        {
            itsRhoUp=pol->GetChargeDensity(Spin::Up  )->DM_RhoAtPoints(itsMesh->Points(), itsPhi);
            itsRhoDn=pol->GetChargeDensity(Spin::Down)->DM_RhoAtPoints(itsMesh->Points(), itsPhi);
        }
        else if (auto sr=dynamic_cast<const ChargeDensity::cSpinResolved_CD*>(cd))
        {   // MATRIX-FREE spin-resolved density: the seed (PolarizedSeedCD, SCFSeedingPlan §10) at
            // iteration 0, and -- the expensive case -- the ρ̃-MIXED density (PolarizedMixCD over
            // FourierMixCD) on EVERY Kerker/Pulay iteration.  Neither carries a D, so both batch through
            // the plain tChargeDensity face instead of the Phi GEMM: an image-summed atomic sum for the
            // seed, a batched inverse FT over the whole {G} for the mixer.  Its OWN bucket because it is
            // a different algorithm at a different cadence -- lumping it into the GEMM hid the fact that
            // the mixed-density sampling, not the GEMM, was the iteration's largest XC cost.
            qchem::report::Timed seed("scf: XC-mesh rho sampling (matrix-free density)");
            itsRhoUp=(*sr->GetChannel(Spin::Up  ))(itsMesh->Points());
            itsRhoDn=(*sr->GetChannel(Spin::Down))(itsMesh->Points());
        }
        else
        {   // spin-agnostic seed: rho_up=rho_down=rho/2 (the molecular HalfDensity rule, cd85d13c)
            if (auto dm=dynamic_cast<const cDM_CD*>(cd))
                itsRhoUp=dm->DM_RhoAtPoints(itsMesh->Points(), itsPhi);
            else
                itsRhoUp=(*cd)(itsMesh->Points());
            itsRhoUp*=0.5;
            itsRhoDn=itsRhoUp;
        }
        if (!itsFold.owner.empty())
        {
            if (itsSigmas.empty())
            {   // §6a W1 grey imposition: the spatial star-average acts per channel independently
                Symmetry::Lattice_3D::SymmetrizeValues(itsFold, itsRhoUp);
                Symmetry::Lattice_3D::SymmetrizeValues(itsFold, itsRhoDn);
            }
            else
            {   // Shubnikov S3: the (ρ,m) pair diagonalizes σ -- the TOTAL is even (plain orbit mean),
                // the MAGNETIZATION is odd (χ-signed mean, zeroed first at the flip-fixed points where
                // the exact projector annihilates it).  Per-channel averaging would be WRONG here: a
                // σ=Flip op maps ρ_up onto ρ_dn, not onto itself.
                rvec_t rho = itsRhoUp + itsRhoDn;
                rvec_t m   = itsRhoUp - itsRhoDn;
                for (size_t g=0; g<itsFlipFixed.size(); ++g) if (itsFlipFixed[g]) m[g]=0.0;
                Symmetry::Lattice_3D::SymmetrizeValues      (itsFold, rho);
                Symmetry::Lattice_3D::SymmetrizeValuesSigned(itsFold, itsSigmas, m);
                itsRhoUp = 0.5*(rho + m);
                itsRhoDn = 0.5*(rho - m);
            }
        }
    }
    return s==Spin::Up ? itsRhoUp : itsRhoDn;
}

// <i|v|j> = Phi^dag diag(w v) Phi over the cached table: scale the rows, one zgemm, hermitize.  (The
// GEMM result is Hermitian up to roundoff; the explicit i<=j fill keeps chmat_t's invariant exactly.)
chmat_t XC_GridEngine::Matrix(const cobs_t* bs, const rvec_t& v) const
{
    const mat_t<dcmplx>& P=Phi(bs);
    const rvec_t&        w=itsMesh->Weights();
    assert(v.size()==P.rows());
    qchem::report::Timed timed("scf: XC-mesh quadrature H_xc (all iterations)");
    mat_t<dcmplx> WP(P.rows(), P.columns());
    auto scale=[&](size_t g)
    {
        const double wv=w[g]*v[g];
        for (size_t i=0; i<P.columns(); i++) WP(g,i)=wv*P(g,i);
    };
    const size_t n=P.columns(), npts=P.rows();
    mat_t<dcmplx> M(n, n, dcmplx(0.0));
    // The per-iteration quadrature GEMM (once per XC term per spin): O(npts n^2), and with the rho
    // sampling threaded it is the SCF loop's largest bucket.  Two things happen here:
    //  * TRIANGULAR: H is Hermitian and chmat_t stores i<=j, so only the UPPER triangle of M is
    //    computed -- column block [j0,j1) needs rows [0,j1) of the left operand, halving the flops.
    //    (The old full-M build then averaged M(i,j) with conj(M(j,i)); that average only symmetrised
    //    ROUNDOFF -- w and v are real, so M is Hermitian exactly in exact arithmetic.  Taking the
    //    upper triangle keeps chmat_t's invariant just as exactly, at ~1e-16 in the elements.)
    //  * THREADED by output column block -- every element M(i,j) is still ONE dot product over ALL
    //    mesh points, accumulated by ONE thread in the serial order, so the matrix is bit-identical at
    //    any thread count (a partition of the OUTPUT, not of the reduction).  Blocks are triangular =
    //    unequal work, hence dynamic scheduling over ~4 blocks per thread.  The row scaling above
    //    parallelises trivially (row g is private to g).
    auto triBlock=[&](size_t j0, size_t nj)
    {
        const size_t j1=j0+nj;
        blazem::submatrix(M,0,j0,j1,nj) = blazem::trans(blazem::conj(blazem::submatrix(P,0,size_t(0),npts,j1)))
                                        * blazem::submatrix(WP,0,j0,npts,nj);
    };
#ifdef QCHEM_OPENMP
    const int nthreads=qchem::WorkerThreads();
    if (nthreads>1)
    {
        #pragma omp parallel for schedule(static) num_threads(nthreads)
        for (size_t g=0; g<npts; g++) scale(g);
        const size_t blk=std::max<size_t>(1, (n+4*size_t(nthreads)-1)/(4*size_t(nthreads)));
        const size_t nb=(n+blk-1)/blk;
        #pragma omp parallel for schedule(dynamic,1) num_threads(nthreads)
        for (size_t b=0; b<nb; b++) triBlock(b*blk, std::min(blk, n-b*blk));
    }
    else
#endif
    {
        for (size_t g=0; g<npts; g++) scale(g);
        triBlock(0, n);
    }
    chmat_t H(n);
    for (size_t i=0; i<n; i++)
        for (size_t j=i; j<n; j++)
            H(i,j)=M(i,j);
    return H;
}

// ---- DeltaFittedVxc ------------------------------------------------------------------------------------------

// Built with the SHARED quadrature engine (the caller builds ONE engine per XC pair -- mesh + Phi tables
// + per-serial rho -- and hands it to both the exchange and the correlation term).
DeltaFittedVxc::DeltaFittedVxc(const xc_t& xc, engine_t engine)
    : itsXc(xc)
    , itsEngine(std::move(engine))
{
    assert(itsEngine);
}

// v_xc(rho_g) pointwise on the engine's shared rho, then the engine's Phi-table quadrature (one GEMM).
chmat_t DeltaFittedVxc::MakeMatrix(const cobs_t* bs, const Spin&, const cChargeDensity* cd) const
{
    const rvec_t& rho=itsEngine->Rho(cd, bs);
    rvec_t v(rho.size());
    for (size_t g=0; g<rho.size(); g++) v[g]=itsXc->GetVxc(rho[g]);
    return itsEngine->Matrix(bs, v);
}

void DeltaFittedVxc::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    const rvec_t& rho=itsEngine->Rho(cd);   // reuses the iteration's table (same density serial)
    const rvec_t& W=itsEngine->Mesh().Weights();
    double exc=0, q=0;
    for (size_t g=0; g<rho.size(); g++)
    {
        const double ro=rho[g];
        exc+=W[g]*itsXc->GetEpsXc(ro)*ro;   // E_xc = integral eps_xc(rho) rho
        q  +=W[g]*ro;
    }
    te.Exc += exc;
    // The mesh-charge leak (the Becke grid's health metric, mirroring PWFittedVxc's grid-charge-lost): the
    // quadrature integral of rho vs the analytic Tr(DS).
    te.GridChargeLost = q - cd->GetTotalCharge();
}

std::ostream& DeltaFittedVxc::Write(std::ostream& os) const
{
    return os << "    XC-mesh exchange-correlation potential v_xc(rho(r)) ("
              << itsEngine->Mesh().size() << " atom-centred points)." << std::endl;
}

// ---- DeltaFittedVxcPol (spin-native exchange, tier 4b) --------------------------------------------------------

DeltaFittedVxcPol::DeltaFittedVxcPol(const xc_t& xc, engine_t engine)
    : itsXc(xc)
    , itsEngine(std::move(engine))
{
    assert(itsXc);
    assert(itsEngine);
}

// v_x^sigma(rho_sigma) pointwise on this block's own channel raster, then the shared Phi quadrature.
chmat_t DeltaFittedVxcPol::MakeMatrix(const cobs_t* bs, const Spin& s, const cChargeDensity* cd) const
{
    assert(s!=Spin::None && "DeltaFittedVxcPol: a polarized term needs an Up/Down spin");
    const rvec_t& rho=itsEngine->RhoPol(cd, s, bs);
    rvec_t v(rho.size());
    for (size_t g=0; g<rho.size(); g++) v[g]=itsXc->GetVxc(rho[g]);
    return itsEngine->Matrix(bs, v);
}

void DeltaFittedVxcPol::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    const rvec_t& up=itsEngine->RhoPol(cd, Spin::Up  );
    const rvec_t& dn=itsEngine->RhoPol(cd, Spin::Down);
    const rvec_t& W =itsEngine->Mesh().Weights();
    double exc=0, q=0;
    for (size_t g=0; g<up.size(); g++)
    {
        exc+=W[g]*(itsXc->GetEpsXc(up[g])*up[g] + itsXc->GetEpsXc(dn[g])*dn[g]);   // E_x = Σ_σ ∫ ε_x(ρ_σ) ρ_σ
        q  +=W[g]*(up[g]+dn[g]);
    }
    te.Exc += exc;
    te.GridChargeLost = q - cd->GetTotalCharge();   // mesh-charge leak (same health metric as DeltaFittedVxc)
}

std::ostream& DeltaFittedVxcPol::Write(std::ostream& os) const
{
    return os << "    XC-mesh SPIN-NATIVE exchange v_x(rho_sigma(r)) ("
              << itsEngine->Mesh().size() << " atom-centred points)." << std::endl;
}

// ---- DeltaFittedVcorrPol (spin-native correlation, tier 4b) ---------------------------------------------------

DeltaFittedVcorrPol::DeltaFittedVcorrPol(const corr_t& corr, engine_t engine)
    : itsCorr(corr)
    , itsEngine(std::move(engine))
{
    assert(itsCorr);
    assert(itsEngine);
}

// v_c^sigma(rho_up,rho_down) couples BOTH channel rasters at every point (through r_s and zeta).
chmat_t DeltaFittedVcorrPol::MakeMatrix(const cobs_t* bs, const Spin& s, const cChargeDensity* cd) const
{
    assert(s!=Spin::None && "DeltaFittedVcorrPol: a polarized term needs an Up/Down spin");
    const rvec_t& up=itsEngine->RhoPol(cd, Spin::Up  , bs);
    const rvec_t& dn=itsEngine->RhoPol(cd, Spin::Down, bs);
    rvec_t v(up.size());
    for (size_t g=0; g<up.size(); g++) v[g]=itsCorr->GetVc(up[g], dn[g], s);
    return itsEngine->Matrix(bs, v);
}

void DeltaFittedVcorrPol::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    const rvec_t& up=itsEngine->RhoPol(cd, Spin::Up  );
    const rvec_t& dn=itsEngine->RhoPol(cd, Spin::Down);
    const rvec_t& W =itsEngine->Mesh().Weights();
    double ec=0;
    for (size_t g=0; g<up.size(); g++)
        ec+=W[g]*itsCorr->GetEpsC(up[g], dn[g])*(up[g]+dn[g]);   // E_c = ∫ ε_c(ρ↑,ρ↓) ρ_total
    te.Exc += ec;
}

std::ostream& DeltaFittedVcorrPol::Write(std::ostream& os) const
{
    return os << "    XC-mesh SPIN-NATIVE correlation v_c^sigma(rho_up,rho_down) ("
              << itsEngine->Mesh().size() << " atom-centred points)." << std::endl;
}

} //namespace
