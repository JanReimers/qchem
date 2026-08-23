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
#include <optional>    // the conditionally-charged sub-buckets of the H_xc quadrature
#include <stdexcept>
module qchem.Hamiltonian.Internal.PWTerms;
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.ChargeDensity.FourierDensity;   // cast cd UP to its reciprocal-space coefficients rho-tilde
import qchem.BasisSet.Orbital_DFT_IBS;         // cast bs UP to the reciprocal-space DFT capability (Hartree/XC)
import qchem.BasisSet.G_FieldEvaluator;    // G_RasterTransform: the fit basis's FFT pair (RhoOnGrid, the BALL route)
import qchem.Pseudopotential.Integrals_Pseudo;   // cast bs ACROSS to the external-PP operator-assembly mixin (Ven_PP_*)
import qchem.Fitting.FunctionFitter;        // Fitting::Factory (both PW fitters) + ProjectedDensity_G / ProjectedScalar_R
import qchem.Structure;                       // Structure::isFinite()/SumFormFactors() -- the G=0 alignment (term-side)
import qchem.Blaze;                            // blazem::zeroH<dcmplx> (the null-PP V_long block)
import qchem.Mesh.Quadrature;                 // qcMesh::Mesh (the Vxc_Quadrature engine's quadrature mesh)
import qchem.Reporting;                       // Timed (the setup/scf timing ledger)
import qchem.Parallel;                         // WorkerThreads (GPW_OMP_THREADS -- the XC-mesh table + quadrature loops)

namespace qchem::Hamiltonian
{

namespace
{
// EXACT narrow for the real-TRIM-block faces (doc/RealComplexPlan.md Step 3c): identity for U=dcmplx,
// bitwise-asserted real part for U=double -- Step 0 made the phases exactly ±1, so a TRIM block's
// complex-assembled term matrix has imag==0.0 EXACTLY (gate: GPW.TRIM_RealBlockMatchesComplexBitwise).
// A sibling of BasisSet::Lattice_3D::ToScalar, duplicated here because qcHamiltonian must not import
// qcLattice_BS (the DAG runs the other way); consolidate into qcMath if a third copy ever appears.
template <class U> hmat_t<U> NarrowExact(const chmat_t& m)
{
    if constexpr (std::is_same_v<U,dcmplx>) return m;
    else
    {
        hmat_t<U> r(m.rows());
        for (size_t i=0;i<m.rows();i++)
            for (size_t j=i;j<m.columns();j++)
            {
                assert(std::imag(m(i,j))==0.0 && "TRIM narrow: imaginary part must be EXACTLY zero (Step 0)");
                r(i,j)=std::real(m(i,j));
            }
        return r;
    }
}
} //anon


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
template <class U> hmat_t<U> Ven_PP_Short::MakeMatrixT(const tobs_t<U>* bs, const Spin&) const
{
    auto pw=dynamic_cast<const Pseudopotential::Integrals_Pseudo<U>*>(bs);
    assert(pw && "Ven_PP_Short requires an Integrals_Pseudo (e.g. plane-wave / GPW) basis");
    // SHORT-range local only.  The LONG (softened-Coulomb) half is Ven_PP_Long and the KB projectors are
    // Ven_PP_NonLocal (the CP2K local-PP split, doc/GPWPlan.md 0e-PP).
    return pw->MakeLocalPotentialShort(&*theStructure, *itsLocal);
}
chmat_t Ven_PP_Short::MakeMatrix (const cobs_t* bs, const Spin& s) const {return MakeMatrixT<dcmplx>(bs,s);}
rsmat_t Ven_PP_Short::MakeMatrixR(const robs_t* bs, const Spin& s) const {return MakeMatrixT<double>(bs,s);}

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

template <class U> hmat_t<U> Ven_PP_NonLocal::MakeMatrixT(const tobs_t<U>* bs, const Spin&) const
{
    auto pw=dynamic_cast<const Pseudopotential::Integrals_Pseudo<U>*>(bs);
    assert(pw && "Ven_PP_NonLocal requires an Integrals_Pseudo (e.g. plane-wave / GPW) basis");
    return pw->MakeSeparablePotential(&*theStructure, *itsSep);
}
chmat_t Ven_PP_NonLocal::MakeMatrix(const cobs_t* bs, const Spin& s) const
{
    if (std::getenv("GPW_NL_PER_L"))
    {   // I0 diagnostic (doc/SphericalLatticePlan.md): bank the per-l blocks once per irrep block.
        // Complex path only -- the itsByL bank is chmat_t; extend if the diagnostic ever needs real blocks.
        auto pw=dynamic_cast<const Pseudopotential::Integrals_Pseudo<dcmplx>*>(bs);
        assert(pw);
        const std::string id=bs->BasisSetID();
        if (itsByLSeen.insert(id).second)
            for (auto& lH : pw->MakeSeparablePotentialByL(&*theStructure, *itsSep))
                itsByL[lH.first].emplace(id, std::move(lH.second));
    }
    return MakeMatrixT<dcmplx>(bs,s);
}
rsmat_t Ven_PP_NonLocal::MakeMatrixR(const robs_t* bs, const Spin& s) const {return MakeMatrixT<double>(bs,s);}

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

template <class U> hmat_t<U> Ven_PP_Long::MakeMatrixT(const tobs_t<U>* bs, const Spin&) const
{
    auto pp=dynamic_cast<const Pseudopotential::Integrals_Pseudo<U>*>(bs);
    assert(pp && "Ven_PP_Long requires an Integrals_Pseudo (e.g. plane-wave / GPW) basis");
    return pp->MakeLocalPotentialLong(&*theStructure, *itsLocal);
}
chmat_t Ven_PP_Long::MakeMatrix (const cobs_t* bs, const Spin& s) const {return MakeMatrixT<dcmplx>(bs,s);}
rsmat_t Ven_PP_Long::MakeMatrixR(const robs_t* bs, const Spin& s) const {return MakeMatrixT<double>(bs,s);}

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

template <class U> hmat_t<U> Vee_Hartree::MakeMatrixT(const tobs_t<U>* bs, const Spin&, const cChargeDensity* cd) const
{
    newCD(cd);   // dirty the Irrep cache if cd is new (the cross-iteration freshness mechanism)
    auto fd=dynamic_cast<const qchem::ChargeDensity::FourierDensity*>(cd);
    assert(fd && "Vee_Hartree requires a FourierDensity (periodic) charge density");
    // The two-axis DFT face (V1.1): a real TRIM block is Orbital_DFT_IBS<double,dcmplx> -- same complex
    // fit side, so the SAME Repulsion3C tensor serves both; only the final block is narrowed.
    auto bft=dynamic_cast<const BasisSet::Orbital_DFT_IBS<U,dcmplx>*>(bs);
    assert(bft && "Vee_Hartree requires a Orbital_DFT_IBS (reciprocal-space DFT) orbital basis");
    // The density contracts D against the basis's D-free Coulomb tensor Repulsion3C (kernel baked) to give
    // V_H(dm) [FORWARD]; the KS matrix <i|V_H|j> = Σ_k V_H(G_k) <i|e^{iG_k}|j> is the BACKWARD contraction of the
    // SAME Repulsion3C tensor over the CD fit basis (its applyAdjoint -- the overlap integrate-back on the fit
    // grid; the Coulomb kernel is forward-only, already in V_H).  So forward AND backward run on the one fit grid
    // (doc/GPWPlan §0e step 2).
    ΔG_Map VH=fd->GetRepulsion3C(*itsFitBasis);
    return NarrowExact<U>(ContractAdjoint(bft->Repulsion3C(*itsFitBasis),
        [&VH](const ivec3_t& dm)->dcmplx { auto it=VH.find(dm); return it==VH.end()?dcmplx(0.0):it->second; }));
}
chmat_t Vee_Hartree::MakeMatrix (const cobs_t* bs, const Spin& s, const cChargeDensity* cd) const {return MakeMatrixT<dcmplx>(bs,s,cd);}
rsmat_t Vee_Hartree::MakeMatrixR(const robs_t* bs, const Spin& s, const cChargeDensity* cd) const {return MakeMatrixT<double>(bs,s,cd);}

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
// A field ALREADY SAMPLED at the quadrature's points, presented as the ProjectedScalar_R the ortho scalar
// fitter consumes (the BALL route's only client).  The values are v_xc(rho(r_g)), computed by the TERM --
// which is where the functional lives, so no qcBasisSet->qcHamiltonian library cycle appears here.
// It is GRID-BOUND: only the ortho fitter samples it, in bulk, on exactly the points it was built for.
class PWVxcField
    : public virtual ScalarFunction<double>
    , public         Fitting::ProjectedScalar_R
{
public:
    PWVxcField(const rvec_t& vals, const qcMesh::Mesh* grid) : itsVals(vals), itsGrid(grid) {}

    // Pointwise is NOT supported: this field carries only grid values, and nothing samples it pointwise (the
    // ortho fitter uses the bulk overload).  Make the grid-bound contract explicit rather than silently wrong.
    virtual double  operator()(const rvec3_t&) const override
        {throw std::logic_error("PWVxcField is grid-bound: sample it in bulk on the fit grid, not pointwise");}
    virtual rvec3_t Gradient  (const rvec3_t&) const override {return rvec3_t(0,0,0);}

    // Bulk: the precomputed values.  The points MUST be the ones they were computed at -- assert IDENTITY
    // (not merely size), so a future diagnostic that samples a different same-cardinality point set fails
    // loudly instead of pairing values to the wrong points (review #2).
    virtual rvec_t  operator()(const rvec3vec_t& rs) const override
    {
        assert(SampledOnGrid(rs) && "PWVxcField: must be sampled on the fit basis's own grid (identity)");
        return itsVals;
    }

    virtual const ScalarFunction<double>* GetScalarFunction() const override {return this;}
private:
    bool SampledOnGrid(const rvec3vec_t& rs) const   // debug-only identity check (the points are cached)
    {
        const rvec3vec_t& g=itsGrid->Points();
        if (rs.size()!=g.size()) return false;
        for (size_t q=0;q<rs.size();q++)
        {
            double dx=rs[q].x-g[q].x, dy=rs[q].y-g[q].y, dz=rs[q].z-g[q].z;
            if (dx*dx+dy*dy+dz*dz > 1e-24) return false;
        }
        return true;
    }
    const rvec_t&        itsVals;   // v_xc at the quadrature points (owned by the caller; field is transient)
    const qcMesh::Mesh*  itsGrid;   // those points, for the identity check
};
} // anonymous

// ---- XC_PairQuadrature: rho by collocation, H by the SAME tensor's raw adjoint -------------------------

// It carries the fit basis (quadrature + collocation key + Overlap3C key) and, for the BALL fallback only,
// an ortho scalar fitter over it.  No grid announcement here: the fit basis self-reports at ITS
// construction, role-labeled (user ruling 2026-08-16).
XC_PairQuadrature::XC_PairQuadrature(fbs_t fb)
    : itsFitBasis(std::move(fb))
    , itsScalarFitter(Fitting::Factory(itsFitBasis))   // the ortho (G-space) scalar fit -- BALL route only
{
    assert(itsFitBasis);
}
XC_PairQuadrature::~XC_PairQuadrature() = default;   // itsScalarFitter's abstract type is complete here

// My points+weights -- PRIVATE.  A raster answers with its corners i/N at weight Omega/Npts; what leaves
// this object is Integrate() and Matrix(), the same two questions the singles strategy answers, which is
// what lets ONE term run on either.
const qcMesh::Mesh& XC_PairQuadrature::Mesh() const
{
    auto* q=dynamic_cast<const BasisSet::Quadrature*>(itsFitBasis.get());
    assert(q && "XC_PairQuadrature: the fit basis must carry its quadrature (BasisSet::Quadrature)");
    return q->Mesh();
}

// rho(r) on the raster for cd -- recomputed only on a new density serial, so the XC pair's two terms and
// their two energies share ONE collocation per iteration (it used to be one per TERM: the same raw
// collocation ran twice an iteration because each term owned its own copy of this).
void XC_PairQuadrature::Refresh(const cChargeDensity* cd) const
{
    assert(cd);
    if (cd->Version()==itsRhoVersion) return;
    itsRhoVersion=cd->Version();
    auto fd=dynamic_cast<const qchem::ChargeDensity::FourierDensity*>(cd);
    assert(fd && "XC_PairQuadrature requires a FourierDensity (periodic) charge density");
    // RAW-FIRST (doc/GPWPlan 0.5(f2)): a collocation-backed density answers with rho_DM(r) directly --
    // pointwise >= 0 for an aufbau D, so the rho>0 guard never bites and the grid calibration can relax
    // (the C=8 driver was the BALL path's Gibbs lobes).  A plane-wave density (or a matrix-free seed)
    // answers EMPTY and we take the ball round trip below -- bit-identical to the pre-f2 path.
    itsRho=fd->GetRhoOnGrid(*itsFitBasis);
    itsRhoIsRaw=(itsRho.size()!=0);
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
                std::string("XC_PairQuadrature: the XC route changed mid-SCF (")
                + (itsLatchedRaw?"RAW -> BALL":"BALL -> RAW")
                + ").  These minimise DIFFERENT functionals -- BALL's ball-projected rho is non-variational "
                  "-- so the optimiser would be chasing a moving target.  The route is a property of the "
                  "orbital basis and must not change once the SCF is running.");
    }
    if (!itsRhoIsRaw)
    {
        // The BALL round trip is a RASTER operation by nature (an inverse FFT of rho-tilde onto the fit
        // basis's own grid), which is why this strategy -- and not the mesh-general singles one -- is
        // where it lives.
        auto* ge=dynamic_cast<const BasisSet::G_RasterTransform*>(itsFitBasis.get());
        assert(ge && "XC_PairQuadrature: the BALL route needs the fit basis's raster transforms");
        itsRho=ge->RhoOnGrid(fd->GetFourierDensity(*itsFitBasis));   // rho-tilde via Overlap3C, onto the FIT grid
    }
    // DIAGNOSTIC (env GPW_XCROUTE): which V_xc route fires this iteration -- RAW (applyRawAdjoint, FD-exact/
    // variational; gate GPW.RawXCConsistencyFD) vs BALL (DoFit/Overlap, non-variational under BallOnly;
    // gate GPW.XCPotentialConsistencyFD).  Answers whether GDM is fighting the ball non-variationality.
    if (std::getenv("GPW_XCROUTE"))
    {
        double rmin=1e300, rmax=-1e300;
        for (double r : itsRho) { rmin=std::min(rmin,r); rmax=std::max(rmax,r); }
        std::cout << "[xc route] " << (itsRhoIsRaw ? "RAW (applyRawAdjoint, variational)"
                                                   : "BALL (DoFit/Overlap, NON-variational under BallOnly)")
                  << "  rho grid=[" << std::scientific << std::setprecision(2) << rmin << ", " << rmax << "]"
                  << " npts=" << itsRho.size() << std::defaultfloat << std::endl;
    }
    if (ReportGridCharge())
    {
        // Grid charge vs analytic charge: the electrons LOST to grid truncation (high-G aliasing of rho).
        // == CP2K's "Electronic density on regular grids: <int rho> <error>" -- a controlled cutoff metric.
        const double qGrid=Integrate(itsRho);   // integral rho_grid d3r
        const double qDM  =cd->GetTotalCharge();                // Tr(D S) (analytic, ~ N)
        std::cout << "[grid charge] integral rho_grid=" << std::fixed << std::setprecision(6) << qGrid
                  << "  Tr(DS)=" << qDM
                  << "  lost=" << std::scientific << std::setprecision(3) << (qGrid-qDM)
                  << std::defaultfloat << std::endl;
        // XC-collapse diagnostic (doc/GPWPlan §0e step 2): the collocated rho's min/max/negative content.
        // If rho rings locally-negative near the sharp F peaks, GetEpsXc(rho<=0)=0 silently drops that
        // eps_xc*rho -- the +7 Ha Exc collapse.  Separates lead (c) [Gibbs ringing + guard: rho_min very
        // negative, big negCharge] from lead (b) [genuinely lower peaks: rho_max small, rho_min ~ 0].
        // NB the eps_xc estimate needs a functional and this object has none (the functional lives in the
        // TERM, which is the right place for it), so the lost-Exc column is reported by magnitude only.
        double rmin=1e300, rmax=-1e300; size_t nneg=0;
        rvec_t negOnly(itsRho.size());
        for (size_t q=0;q<itsRho.size();q++)
        {
            const double r=itsRho[q];
            rmin=std::min(rmin,r); rmax=std::max(rmax,r);
            negOnly[q] = r<0.0 ? r : 0.0;                                  // negative-charge density
            if (r<0.0) ++nneg;
        }
        std::cout << "[xc grid] rho_min=" << std::scientific << std::setprecision(3) << rmin
                  << " rho_max=" << rmax << " neg-frac=" << double(nneg)/double(itsRho.size())
                  << " negCharge=" << Integrate(negOnly)
                  << std::defaultfloat << std::endl;
    }
}

double XC_PairQuadrature::Integrate(const rvec_t& f) const {return qcMesh::Integrate(Mesh(), f);}
size_t XC_PairQuadrature::NumPoints() const {return Mesh().size();}

const rvec_t& XC_PairQuadrature::Rho(const cChargeDensity* cd) const {Refresh(cd); return itsRho;}

const rvec_t& XC_PairQuadrature::RhoPol(const cChargeDensity*, const Spin&) const
{
    throw std::logic_error("XC_PairQuadrature: the pair (collocation) route is not spin-native -- there is "
        "no per-channel rho_sigma collocation, and a per-spin pair route with its own rho caches is not "
        "designed.  Use the delta/singles quadrature (VxcFit::Delta), which is spin-native on any mesh; "
        "VxcFit::Auto already selects it for a polarized run.");
}

// <i|v|j>, THE EXACT ADJOINT of whichever route Refresh took (which is why both live on one object):
//   RAW  -- the raw adjoint of the SAME tensor whose applyRaw produced itsRho (box-truncation per level +
//           the analytic gather), so H_xc == dE_xc/dD of the ONE raw discrete functional to machine
//           precision (gate: GPW.RawXCConsistencyFD).  No ball fit anywhere.
//   BALL -- fit v_xc on the {G} ball and contract: the legacy, non-variational pairing.
// Two-axis face (V1.1): the tensor follows TFit==dcmplx for both block scalars; narrow at the end.
template <class U> hmat_t<U> XC_PairQuadrature::MatrixT(const tobs_t<U>* bs, const rvec_t& v) const
{
    assert(itsRhoVersion!=size_t(-1) && "XC_PairQuadrature::Matrix before Rho: the adjoint must follow the "
           "collocation that fixed the route");
    if (itsRhoIsRaw)
    {
        const auto& orb=dynamic_cast<const BasisSet::Orbital_DFT_IBS<U,dcmplx>&>(*bs);   // genuine "is it?" cross-cast (throws)
        const Projector3<dcmplx>& g=orb.Overlap3C(*itsFitBasis);
        assert(g.applyRawAdjoint && "raw rho without a raw adjoint: Overlap3C must carry both");
        return NarrowExact<U>(g.applyRawAdjoint(v));
    }
    if constexpr (std::is_same_v<U,dcmplx>)
    {
        itsScalarFitter->DoFit(PWVxcField(v, &Mesh()));   // our quadrature IS the fitter's: one owner, the fit basis
        // The COMPLEX contraction face (ISP): this fitter's raster basis serves Bloch blocks; the
        // real-block branch above never reaches here (it takes the raw adjoint or throws).
        return dynamic_cast<const Fitting::FitContraction<dcmplx>&>(*itsScalarFitter).Overlap(bs);
    }
    else
        // The legacy BALL-fit route's fitter Overlap is dcmplx-faced; nothing real-block reaches it (the
        // raw feed is the production default, and a real TRIM block implies a GPW lineage).  Loud, not silent.
        throw std::logic_error("XC_PairQuadrature: a real TRIM block on the legacy ball-fit XC route is not "
                               "wired -- use the raw collocated feed (the default) or the delta quadrature");
}
chmat_t XC_PairQuadrature::Matrix(const cobs_t* bs, const rvec_t& v) const {return MatrixT<dcmplx>(bs,v);}
rsmat_t XC_PairQuadrature::Matrix(const robs_t* bs, const rvec_t& v) const {return MatrixT<double>(bs,v);}

// CAPABILITY DECIDES (doc/OpenWork.md): a delta basis carries points and nothing else -> singles; a
// raster-backed one carries the FFT transforms and keys the 3-centre tensor -> pair.  One decision, taken
// once, and latched for the run by the simple fact that the Hamiltonian builds this object once.
std::shared_ptr<const XC_Quadrature>
MakeXCQuadrature(const std::shared_ptr<const BasisSet::cFIT_SF_ABS>& fb)
{
    assert(fb);
    // CAPABILITY decides.  A basis that carries the {r}<->{G} transforms is raster-backed, so its
    // collocation pair (Overlap3C's applyRaw/applyRawAdjoint) exists and the PAIR route is available --
    // and preferred, being the production GPW path.  Anything else can only be contracted through a Phi
    // table: SINGLES.  Note this asks what the basis CAN do, never what it IS.
    if (dynamic_cast<const BasisSet::G_RasterTransform*>(fb.get()))
        return std::make_shared<const XC_PairQuadrature>(fb);
    auto fit=std::dynamic_pointer_cast<const BasisSet::cFIT_SF_Delta>(fb);
    assert(fit && "MakeXCQuadrature: a non-raster Vxc fit basis must be able to answer its own quadrature "
                  "(integrate / sample / symmetrize) -- i.e. it must be a FIT_SF_Delta");
    return std::make_shared<const XC_SinglesQuadrature>(fit);
}

// ---- XC_SinglesQuadrature: the pair-shared mesh + Phi tables + per-serial rho ----------------------------------

XC_SinglesQuadrature::XC_SinglesQuadrature(fit_t fit)
    : itsFit(std::move(fit))
{
    assert(itsFit && "XC_SinglesQuadrature: the delta fit basis IS the quadrature -- it cannot be null");
    // The bundle's own invariants (fold partitions the mesh, sigmas need a fold, flags cover it) are
    // checked where the bundle becomes an object -- in the basis's ctor -- and the star-average it will
    // apply announces itself there too (providers self-report).  Nothing left to validate here.
}

// The quadrature questions go THROUGH the basis: it owns points, weights and the site partition, and
// answers integrals rather than handing them over.
double XC_SinglesQuadrature::Integrate(const rvec_t& f) const {return itsFit->Integrate(f);}
size_t XC_SinglesQuadrature::NumPoints() const {return itsFit->NumPoints();}

// THE XC DM-rho REPAIR (doc/OpenWork.md, the factored-rho section).  Under rho-tilde mixing the density
// driving the Fock build is a G-space FIELD, so XC has to inverse-transform a truncated series at every
// mesh point -- 5.19 s per sampling on MnO (51% of the run) AND 8.5-18% of the points come back with
// rho<0, where the functionals guard `if (rho>0)` and silently contribute nothing.  But the mixer BUILT
// that field from a density matrix and now retains it (cDM_Sourced_CD), so the exact density is one
// cross-cast away: 0.042 s and rho>=0 by construction.  Hartree keeps the preconditioned field -- Poisson
// is LINEAR and diagonal in G -- while XC, a NONLINEAR POINTWISE functional, gets the cusp.  At the fixed
// point they agree, so this changes the SCF TRAJECTORY, not the answer.
// GPW_XC_DM_SOURCE=1 to arm it.  OPT-IN: a trajectory change must earn its place against banked recipes.
bool UseDMSource()
{
    static const bool on=[]{ const char* e=std::getenv("GPW_XC_DM_SOURCE"); return e && std::atoi(e)!=0; }();
    return on;
}
//! The exact density behind a field-backed one + THE DAMPING THAT FIELD APPLIED -- always together, because
//! taking the first without the second is the half-damped map (see cDM_Sourced_CD::EffectiveAlpha).
struct ExactSource { const cDM_CD* cd=nullptr; double alpha=0.0; explicit operator bool() const {return cd;} };
ExactSource ExactSourceOf(const qchem::ChargeDensity::tChargeDensity<dcmplx>* cd)
{
    if (!UseDMSource() || !cd) return {};
    auto* src=dynamic_cast<const qchem::ChargeDensity::cDM_Sourced_CD*>(cd);
    if (!src) return {};
    return {src->DMSource().get(), src->EffectiveAlpha()};
}
// GPW_XC_DM_MIX overrides alpha_eff for CONTROLS only (=1 reproduces the undamped route); unset = use the
// mix's own.  GPW_XC_DM_BOOST scales it: alpha_eff came out ~0.20 on NaF and ~0.35 on MnO -- measured, but
// low against fractions those cells tolerate, and MnO converged 53 -> 39 iterations at boost 2.  NB f_K<=1
// forces alpha_eff<=alpha, so a boost >1/mean(f) puts XC ABOVE the mixer's own alpha -- defensible (XC's
// response kernel is finite at G->0, unlike Hartree's 4pi/G^2) but outside the preconditioner's bracket,
// hence a knob and not a default.
double DMSourceMixOverride()
{
    static const double a=[]{ const char* e=std::getenv("GPW_XC_DM_MIX"); return e ? std::atof(e) : -1.0; }();
    return a;
}
double DMSourceMixBoost()
{
    static const double b=[]{ const char* e=std::getenv("GPW_XC_DM_BOOST"); return e ? std::atof(e) : 1.0; }();
    return b;
}
//! Blend \a fresh into \a running at the mix's own alpha (outside (0,1) => passthrough, which is the right
//! bootstrap AND the right answer for an unmixed or overshooting step).  Non-negativity survives: a convex
//! combination of non-negative rasters is non-negative.
void DampXCChannel(rvec_t& running, const rvec_t& fresh, double alphaEff)
{
    const double ov=DMSourceMixOverride();
    const double a =(ov>=0.0) ? ov : DMSourceMixBoost()*alphaEff;
    static const bool trace=std::getenv("GPW_XC_ALPHA")!=nullptr;
    if (trace) std::cout<<"[XC alpha] alpha_eff="<<alphaEff<<(ov>=0.0?"  (OVERRIDDEN by GPW_XC_DM_MIX)":"")
                        <<"  boost="<<DMSourceMixBoost()
                        <<"  applied="<<((a>0.0&&a<1.0)?a:1.0)<<std::endl;
    if (a<=0.0 || a>=1.0 || running.size()!=fresh.size()) { running=fresh; return; }
    for (size_t g=0; g<fresh.size(); ++g) running[g]=(1.0-a)*running[g]+a*fresh[g];
}

// GPW_RHO_NEGATIVE=1: IS THE SAMPLED rho EVER NEGATIVE, AND WHAT DOES THAT COST?  A direct A/B between the
// two routes feeding this mesh, because they differ in KIND: the DM route gives rho = ||L^dag Phi||^2, a
// sum of squares, so rho>=0 BY CONSTRUCTION; the rho-tilde route inverse-transforms a TRUNCATED Fourier
// series, which RINGS, and rings NEGATIVE where the true density is sharpest -- the nuclear cusps this mesh
// exists to integrate.  And it is not loud: SlaterExchange::GetVxc guards `if (ro > 0.0)` and returns 0, so
// such a point contributes NOTHING to v_xc or E_xc with no diagnostic.  Hence the WEIGHTED MASS, not just
// the count -- that is the E_xc being silently dropped.
// Takes the QUADRATURE, not its mesh: the two masses are integrals, so it asks for them (R1.0) and never
// touches a weight -- which is what let the last Mesh() use out of this diagnostic.
void ReportNegativeRho(const BasisSet::cFIT_SF_Delta& q, const rvec_t& rho, const char* route)
{
    static const bool on=std::getenv("GPW_RHO_NEGATIVE")!=nullptr;
    if (!on || rho.size()==0) return;
    rvec_t negOnly(rho.size(), 0.0), absRho(rho.size());
    size_t cnt=0; double minRho=0.0;
    for (size_t g=0; g<rho.size(); g++)
    {
        absRho[g]=std::fabs(rho[g]);
        if (rho[g]>=0.0) continue;
        cnt++; negOnly[g]=rho[g]; minRho=std::min(minRho,rho[g]);
    }
    const double negMass=q.Integrate(negOnly), absMass=q.Integrate(absRho);
    std::cout<<"[rho<0] route="<<route<<"  points="<<cnt<<" of "<<rho.size()
             <<" ("<<100.0*double(cnt)/double(rho.size())<<"%)"
             <<"  negative mass="<<negMass<<" e ("<<100.0*std::fabs(negMass)/std::max(absMass,1e-30)
             <<"% of integral|rho|)  min rho="<<minRho
             <<"   [these points contribute ZERO to E_xc -- GetVxc guards rho>0]"<<std::endl;
}

// rho at the mesh points, once per density serial for the WHOLE pair: the density GEMMs the cached
// tables against its private D (DM_RhoAtPoints; blocks not yet tabled self-evaluate pointwise -- first
// pass only).  A non-DM density (no DM face) falls back to the pointwise ScalarFunction sweep.
const rvec_t& XC_SinglesQuadrature::Rho(const cChargeDensity* cd) const
{
    assert(cd);
    // R2.9(i): the scalar and spin-resolved caches do not cross-invalidate (see the \warning on the
    // members).  An engine belongs to ONE xc/correlation pair, and a pair is either polarized or not, so
    // only one of the two routes is ever driven.  Pin it here rather than trusting the comment.
    assert(itsPolVersion==size_t(-1) && "XC_SinglesQuadrature: this engine already served RhoPol -- the scalar "
           "and spin-resolved rho caches have no cross-invalidation, so one of them would go stale");
    if (cd->Version()==itsRhoVersion) return itsRho;
    itsRhoVersion=cd->Version();
    qchem::report::Timed timed("scf: XC-mesh rho sampling (all iterations)");
    if (auto dm=dynamic_cast<const cDM_CD*>(cd))
    {
        itsRho=dm->DM_RhoAtPoints(*itsFit);   // the density asks the basis for each block's table
        ReportNegativeRho(*itsFit, itsRho, "DM");
    }
    else if (auto ex=ExactSourceOf(cd))
    {   // THE REPAIR, unpolarized sibling of the RhoPol branch below: the field retains the D it was mixed
        // from, so XC samples that instead of inverse-transforming a truncated series.
        if (ex.cd->Version()==itsSrcVersion)
            std::cerr<<"[XC DM-source] ** STALE: the mixed field advanced to serial "<<itsRhoVersion
                     <<" but its retained density matrix did NOT (still "<<itsSrcVersion
                     <<") -- V_xc is being built from the PREVIOUS iteration's D."<<std::endl;
        itsSrcVersion=ex.cd->Version();
        DampXCChannel(itsXCMix, ex.cd->DM_RhoAtPoints(*itsFit), ex.alpha);
        itsRho=itsXCMix;   // the running mix lives in its OWN buffer -- see the \warning on itsXCMix
        ReportNegativeRho(*itsFit, itsRho, "DM-source");
    }
    else
    {
        itsRho=itsFit->Sample(*cd);   // non-DM (mixed rho-tilde / seed) densities batch through the basis
        ReportNegativeRho(*itsFit, itsRho, "matrix-free");
    }
    itsFit->Symmetrize(itsRho);   // §6a W1: the basis's own orbit-mean projector (no-op on a free run)
    return itsRho;
}

// The real-block ensure siblings (3c-3): build the real block's OWN typed table first (PhiR), then the
// shared sampling path -- the exact mirror of the complex ensureBlock argument.

// The spin-resolved sibling of Rho: the {up,down} PAIR is cached under ONE density serial (a polarized
// density's Version() forwards to its Up child -- a single scalar cache would hand the Up raster to the
// Down channel).  A cPolarized_CD answers per channel (each channel composite GEMMs its own D against the
// SHARED Phi tables); a spin-agnostic density (the seed) collapses to rho/2 per channel, so the first
// iterations run the exact unpolarized collapse (v^sigma(rho/2,rho/2)=v^P(rho)).
const rvec_t& XC_SinglesQuadrature::RhoPol(const cChargeDensity* cd, const Spin& s) const
{
    assert(cd);
    assert(itsRhoVersion==size_t(-1) && "XC_SinglesQuadrature: this engine already served the scalar Rho -- the "
           "two rho caches have no cross-invalidation, so one of them would go stale");
    assert(s!=Spin::None && "XC_SinglesQuadrature::RhoPol: ask for a channel, not the total");
    if (cd->Version()!=itsPolVersion)
    {
        itsPolVersion=cd->Version();
        qchem::report::Timed timed("scf: XC-mesh rho sampling (all iterations)");
        if (auto pol=dynamic_cast<const ChargeDensity::cPolarized_CD*>(cd))
        {
            itsRhoUp=pol->GetChargeDensity(Spin::Up  )->DM_RhoAtPoints(*itsFit);
            itsRhoDn=pol->GetChargeDensity(Spin::Down)->DM_RhoAtPoints(*itsFit);
            ReportNegativeRho(*itsFit, itsRhoUp, "DM(up)");
            ReportNegativeRho(*itsFit, itsRhoDn, "DM(dn)");
        }
        else if (auto sr=dynamic_cast<const ChargeDensity::cSpinResolved_CD*>(cd))
        {   // MATRIX-FREE spin-resolved density: the seed (PolarizedSeedCD, SCFSeedingPlan §10) at
            // iteration 0, and -- the expensive case -- the ρ̃-MIXED density (PolarizedMixCD over
            // FourierMixCD) on EVERY Kerker/Pulay iteration.  Neither carries a D, so both batch through
            // the plain tChargeDensity face instead of the Phi GEMM: an image-summed atomic sum for the
            // seed, a batched inverse FT over the whole {G} for the mixer.  Its OWN bucket because it is
            // a different algorithm at a different cadence -- lumping it into the GEMM hid the fact that
            // the mixed-density sampling, not the GEMM, was the iteration's largest XC cost.
            // THE REPAIR: each channel of a rho-tilde-mixed density retains the DM-backed channel it was
            // mixed from, so ask for the exact one FIRST and fall back per CHANNEL -- the seed has no source
            // on either while a mixed density has one on both, so a pair-level test would be right today and
            // wrong the first time they differ.
            const ExactSource exUp=ExactSourceOf(sr->GetChannel(Spin::Up  ));
            const ExactSource exDn=ExactSourceOf(sr->GetChannel(Spin::Down));
            if (exUp && exDn)
            {
                // STALENESS GUARD (see itsSrcVersion): the field's serial just advanced, so the SOURCE's
                // must have too, or XC gets last iteration's D under a cache that believes it is fresh.
                if (exUp.cd->Version()==itsSrcVersion)
                    std::cerr<<"[XC DM-source] ** STALE: the mixed field advanced to serial "<<itsPolVersion
                             <<" but its retained density matrix did NOT (still "<<itsSrcVersion
                             <<") -- V_xc is being built from the PREVIOUS iteration's D."<<std::endl;
                itsSrcVersion=exUp.cd->Version();
                rvec_t up=exUp.cd->DM_RhoAtPoints(*itsFit);
                rvec_t dn=exDn.cd->DM_RhoAtPoints(*itsFit);
                DampXCChannel(itsXCMixUp, up, exUp.alpha);   // match the damping Hartree gets, so the map
                DampXCChannel(itsXCMixDn, dn, exDn.alpha);   //   is not half-damped
                itsRhoUp=itsXCMixUp; itsRhoDn=itsXCMixDn;
                ReportNegativeRho(*itsFit, itsRhoUp, "DM-source(up)");
                ReportNegativeRho(*itsFit, itsRhoDn, "DM-source(dn)");
            }
            else
            {
            qchem::report::Timed seed("scf: XC-mesh rho sampling (matrix-free density)");
            itsRhoUp=itsFit->Sample(*sr->GetChannel(Spin::Up  ));
            itsRhoDn=itsFit->Sample(*sr->GetChannel(Spin::Down));
            ReportNegativeRho(*itsFit, itsRhoUp, "matrix-free(up)");
            ReportNegativeRho(*itsFit, itsRhoDn, "matrix-free(dn)");
            }
        }
        else
        {   // spin-agnostic seed: rho_up=rho_down=rho/2 (the molecular HalfDensity rule, cd85d13c)
            if (auto dm=dynamic_cast<const cDM_CD*>(cd))
                itsRhoUp=dm->DM_RhoAtPoints(*itsFit);
            else
                itsRhoUp=itsFit->Sample(*cd);
            itsRhoUp*=0.5;
            itsRhoDn=itsRhoUp;
        }
        {
            // The basis projects the (ρ,m) PAIR: it knows whether the run is magnetic (σ tags -> ρ even,
            // m odd with the flip-fixed zeros), grey (each argument averaged independently) or free (a
            // no-op).  This used to be a three-way branch here over the fold's raw contents.
            rvec_t rho = itsRhoUp + itsRhoDn;
            rvec_t m   = itsRhoUp - itsRhoDn;
            itsFit->SymmetrizeSpin(rho, m);
            itsRhoUp = 0.5*(rho + m);
            itsRhoDn = 0.5*(rho - m);
        }
        // THE OBSERVABLE, reported where it is free (doc/OpenWork.md Step 0a).  This is the ONE place that
        // knows a NEW density has just been sampled on an atom-partitioned mesh, so the per-site integrated
        // moments cost a block sum over data already in hand -- and reporting them here means EVERY
        // polarized atom-centred run gets them, not just the one test that used to fake it with a point
        // probe.  Units: electrons (x mu_B for the magnetic moment).  Named partition, because until
        // Bader's zero-flux basins land the number is partition-dependent.
        EmitSiteMoments();
    }
    return s==Spin::Up ? itsRhoUp : itsRhoDn;
}

// One line + one report entry per NEW density, from inside the serial-advance branch above.  Silent when
// the mesh carries no site partition (a uniform grid has no atomic basins to integrate over).
void XC_SinglesQuadrature::EmitSiteMoments() const
{
    const rvec_t mu=itsFit->SiteIntegrals(rvec_t(itsRhoUp-itsRhoDn));
    if (mu.size()==0) return;                      // no site partition on this quadrature
    double net=0.0, absSum=0.0;
    for (size_t a=0;a<mu.size();a++) { net+=mu[a]; absSum+=std::fabs(mu[a]); }
    if (absSum < 1e-8) return;                    // an unpolarized density has nothing to say
    qchem::report::json j;
    j["partition"]="Becke";                        // NOT canonical -- see the header's Bader note
    j["units"]="electrons";
    // Index loop, not the iterator-pair ctor: std cannot see Blaze's exported iterator op==/op!=
    // across the module boundary (CLAUDE.md "Includes & types").
    std::vector<double> v(mu.size());
    for (size_t a=0;a<mu.size();a++) v[a]=mu[a];
    j["mu"]=v;  j["net"]=net;
    qchem::report::EmitAt("scf", "siteMoments", j);
    if (std::getenv("QCHEM_SITE_MOMENTS"))
    {
        std::cout<<"[site moments] Becke-partitioned Integral w_A (rho_up-rho_dn) d3r [e]:";
        for (size_t a=0;a<mu.size();a++) std::cout<<"  "<<a<<":"<<mu[a];
        std::cout<<"   net="<<net<<std::endl;
    }
}

// The per-site INTEGRATED moment (see the header): mu_A = Integral w_A(r) [rho_up - rho_dn] d3r, in
// electrons.  Free -- RhoPol has already sampled both channels for this density serial (and cached them),
// and the mesh's weights already carry each site's Becke partition w_A, so this is a block sum over data
// in hand.  An unpolarized density gives exactly zero (rho_up == rho_dn by the HalfDensity collapse),
// which is the honest answer, not a special case.
rvec_t XC_SinglesQuadrature::SiteMoments(const cChargeDensity* cd) const
{
    assert(cd);
    const rvec_t& up=RhoPol(cd, Spin::Up);
    const rvec_t& dn=RhoPol(cd, Spin::Down);
    return itsFit->SiteIntegrals(rvec_t(up-dn));   // empty when the quadrature has no site partition
}

// <i|v|j>: ASK THE BASIS.  It owns the points, the weights and the Phi table, so the whole quadrature
// -- Phi^dag diag(w v) Phi -- is its operation; this strategy only decides WHICH v to hand it.  That is
// what closed the last weight/coordinate escape (doc/CleanupCandidates.md R1.0 increment 2).
chmat_t XC_SinglesQuadrature::Matrix(const cobs_t* bs, const rvec_t& v) const {return itsFit->Quadrature(*bs, v);}
rsmat_t XC_SinglesQuadrature::Matrix(const robs_t* bs, const rvec_t& v) const {return itsFit->Quadrature(*bs, v);}

// ---- Vxc_Quadrature ------------------------------------------------------------------------------------------

// Built with the SHARED quadrature engine (the caller builds ONE engine per XC pair -- mesh + Phi tables
// + per-serial rho -- and hands it to both the exchange and the correlation term).
Vxc_Quadrature::Vxc_Quadrature(const xc_t& xc, quad_t quad)
    : itsXc(xc)
    , itsQuad(std::move(quad))
{
    assert(itsQuad);
}

// v_xc(rho_g) pointwise on the engine's shared rho, then the engine's Phi-table quadrature (one GEMM).
// ONE body per term for both block scalars (Step 3c): the rho raster is block-independent; only the
// ensure hint (complex map) and the final quadrature (typed Phi table) differ, both handled below.
template <class U> hmat_t<U> Vxc_Quadrature::MakeMatrixT(const tobs_t<U>* bs, const Spin&, const cChargeDensity* cd) const
{
    const rvec_t& rho=itsQuad->Rho(cd);
    rvec_t v(rho.size());
    for (size_t g=0; g<rho.size(); g++) v[g]=itsXc->GetVxc(rho[g]);
    return itsQuad->Matrix(bs, v);
}
chmat_t Vxc_Quadrature::MakeMatrix (const cobs_t* bs, const Spin& s, const cChargeDensity* cd) const {return MakeMatrixT<dcmplx>(bs,s,cd);}
rsmat_t Vxc_Quadrature::MakeMatrixR(const robs_t* bs, const Spin& s, const cChargeDensity* cd) const {return MakeMatrixT<double>(bs,s,cd);}

void Vxc_Quadrature::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    const rvec_t& rho=itsQuad->Rho(cd);   // reuses the iteration's table (same density serial)
    rvec_t exc(rho.size());
    for (size_t g=0; g<rho.size(); g++) exc[g]=itsXc->GetEpsXc(rho[g])*rho[g];
    const double q=itsQuad->Integrate(rho);
    te.Exc += itsQuad->Integrate(exc);     // E_xc = integral eps_xc(rho) rho, on the quadrature's weights
    // The mesh-charge leak (the quadrature's health metric -- CP2K's grid-charge-lost readout): the
    // quadrature integral of rho vs the analytic Tr(DS).
    te.GridChargeLost = q - cd->GetTotalCharge();
}

std::ostream& Vxc_Quadrature::Write(std::ostream& os) const
{
    return os << "    XC-mesh exchange-correlation potential v_xc(rho(r)) ("
              << itsQuad->NumPoints() << " atom-centred points)." << std::endl;
}

// ---- Vxc_QuadraturePol (spin-native exchange, tier 4b) --------------------------------------------------------

Vxc_QuadraturePol::Vxc_QuadraturePol(const xc_t& xc, quad_t quad)
    : itsXc(xc)
    , itsQuad(std::move(quad))
{
    assert(itsXc);
    assert(itsQuad);
}

// v_x^sigma(rho_sigma) pointwise on this block's own channel raster, then the shared Phi quadrature.
template <class U> hmat_t<U> Vxc_QuadraturePol::MakeMatrixT(const tobs_t<U>* bs, const Spin& s, const cChargeDensity* cd) const
{
    assert(s!=Spin::None && "Vxc_QuadraturePol: a polarized term needs an Up/Down spin");
    const rvec_t& rho=itsQuad->RhoPol(cd, s);
    rvec_t v(rho.size());
    for (size_t g=0; g<rho.size(); g++) v[g]=itsXc->GetVxc(rho[g]);
    return itsQuad->Matrix(bs, v);
}
chmat_t Vxc_QuadraturePol::MakeMatrix (const cobs_t* bs, const Spin& s, const cChargeDensity* cd) const {return MakeMatrixT<dcmplx>(bs,s,cd);}
rsmat_t Vxc_QuadraturePol::MakeMatrixR(const robs_t* bs, const Spin& s, const cChargeDensity* cd) const {return MakeMatrixT<double>(bs,s,cd);}

void Vxc_QuadraturePol::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    const rvec_t& up=itsQuad->RhoPol(cd, Spin::Up  );
    const rvec_t& dn=itsQuad->RhoPol(cd, Spin::Down);
    rvec_t exc(up.size());
    for (size_t g=0; g<up.size(); g++)
        exc[g]=itsXc->GetEpsXc(up[g])*up[g] + itsXc->GetEpsXc(dn[g])*dn[g];   // E_x = Σ_σ ∫ ε_x(ρ_σ) ρ_σ
    const double q=itsQuad->Integrate(rvec_t(up+dn));
    te.Exc += itsQuad->Integrate(exc);
    te.GridChargeLost = q - cd->GetTotalCharge();   // mesh-charge leak (same health metric as Vxc_Quadrature)
}

std::ostream& Vxc_QuadraturePol::Write(std::ostream& os) const
{
    return os << "    XC-mesh SPIN-NATIVE exchange v_x(rho_sigma(r)) ("
              << itsQuad->NumPoints() << " atom-centred points)." << std::endl;
}

// ---- Vcorr_QuadraturePol (spin-native correlation, tier 4b) ---------------------------------------------------

Vcorr_QuadraturePol::Vcorr_QuadraturePol(const corr_t& corr, quad_t quad)
    : itsCorr(corr)
    , itsQuad(std::move(quad))
{
    assert(itsCorr);
    assert(itsQuad);
}

// v_c^sigma(rho_up,rho_down) couples BOTH channel rasters at every point (through r_s and zeta).
template <class U> hmat_t<U> Vcorr_QuadraturePol::MakeMatrixT(const tobs_t<U>* bs, const Spin& s, const cChargeDensity* cd) const
{
    assert(s!=Spin::None && "Vcorr_QuadraturePol: a polarized term needs an Up/Down spin");
    const rvec_t& up=itsQuad->RhoPol(cd, Spin::Up  );
    const rvec_t& dn=itsQuad->RhoPol(cd, Spin::Down);
    rvec_t v(up.size());
    for (size_t g=0; g<up.size(); g++) v[g]=itsCorr->GetVc(up[g], dn[g], s);
    return itsQuad->Matrix(bs, v);
}
chmat_t Vcorr_QuadraturePol::MakeMatrix (const cobs_t* bs, const Spin& s, const cChargeDensity* cd) const {return MakeMatrixT<dcmplx>(bs,s,cd);}
rsmat_t Vcorr_QuadraturePol::MakeMatrixR(const robs_t* bs, const Spin& s, const cChargeDensity* cd) const {return MakeMatrixT<double>(bs,s,cd);}

void Vcorr_QuadraturePol::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    const rvec_t& up=itsQuad->RhoPol(cd, Spin::Up  );
    const rvec_t& dn=itsQuad->RhoPol(cd, Spin::Down);
    rvec_t ec(up.size());
    for (size_t g=0; g<up.size(); g++) ec[g]=itsCorr->GetEpsC(up[g], dn[g])*(up[g]+dn[g]);
    te.Exc += itsQuad->Integrate(ec);   // E_c = ∫ ε_c(ρ↑,ρ↓) ρ_total
}

std::ostream& Vcorr_QuadraturePol::Write(std::ostream& os) const
{
    return os << "    XC-mesh SPIN-NATIVE correlation v_c^sigma(rho_up,rho_down) ("
              << itsQuad->NumPoints() << " atom-centred points)." << std::endl;
}

} //namespace
