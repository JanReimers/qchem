// File: Hamiltonian/Internal/Imp/PWTerms_Hartree.C  the Hartree term V_H[rho_elec] -- the one density-dependent Coulomb term.
//
// One implementation unit of module qchem.Hamiltonian.Internal.PWTerms.  Split 2026-09-08 out of a
// single 1213-line Imp/PWTerms.C (user: "PWTerms.C is huge, again doing too many things") into the
// interface-plus-many-Imp-units shape Internal/Terms.C has always had.  Helpers shared by more than
// one unit (NarrowExact, SampledField) live in the module INTERFACE's non-exported section, which is
// exactly what module-internal linkage is for -- they are visible to every unit of this module and to
// nothing outside it.
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
import qchem.Diagnostics;                     // the ONE diagnostics registry (D-ENV)
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.ChargeDensity.FourierDensity;   // cast cd UP to its reciprocal-space coefficients rho-tilde
import qchem.BasisSet.Orbital_DFT_IBS;         // cast bs UP to the reciprocal-space DFT capability (Hartree/XC)
import qchem.BasisSet.Transition_DFT_IBS;      // the (k+q, k) pair's adjoint gather (the response, R3)
import qchem.ChargeDensity.TransitionDensity;  // TransitionFourierDensity (the response's Hartree face)
import qchem.BasisSet.Orbital_PP_IBS;             // cast bs ACROSS to the species-field integral service (Ven_PP_*)
import qchem.Fitting.FunctionFitter;        // Fitting::Factory (both PW fitters) + ProjectedDensity_G / ProjectedScalar_R
import qchem.Blaze;                            // blazem::zeroH<dcmplx> (the null-PP V_long block)
import qchem.Mesh.Quadrature;                 // qcMesh::Mesh (the Vxc_Quadrature engine's quadrature mesh)
import qchem.Reporting;                       // Timed (the setup/scf timing ledger)
import qchem.BasisSet.G_FieldEvaluator;


namespace qchem::Hamiltonian
{

// CD fit basis and calling IrrepCD::GetRepulsion3C(fbs).  Pure V_H[rho_elec]: no structure, no PP model.
Vee_Hartree::Vee_Hartree(fbs_t fb)
    : itsFitBasis(fb)
{}

template <class U> hmat_t<U> Vee_Hartree::MakeMatrixT(const tobs_t<U>* bs, const Spin&, DensityID id) const
{
    // The two-axis DFT face (V1.1): a real TRIM block is Orbital_DFT_IBS<double,dcmplx> -- same complex
    // fit side, so the SAME Repulsion3C tensor serves both; only the final block is narrowed.
    auto bft=dynamic_cast<const BasisSet::Orbital_DFT_IBS<U,dcmplx>*>(bs);
    assert(bft && "Vee_Hartree requires a Orbital_DFT_IBS (reciprocal-space DFT) orbital basis");
    // The density contracts D against the basis's D-free Coulomb tensor Repulsion3C (kernel baked) to give
    // V_H(dm) [FORWARD]; the KS matrix <i|V_H|j> = Σ_k V_H(G_k) <i|e^{iG_k}|j> is the BACKWARD contraction of the
    // SAME Repulsion3C tensor over the CD fit basis (its applyAdjoint -- the overlap integrate-back on the fit
    // grid; the Coulomb kernel is forward-only, already in V_H).  So forward AND backward run on the one fit grid
    // (doc/GPWPlan §0e step 2).  The FIELD is the one the refresh built (V1.41 stage 2): a block READS it, and
    // \a cd only names which density it must be.
    //
    // The adjoint is the G-space integrator's, a view of the same cached Repulsion3C tensor.  (3b moves V_H out of the
    // term and into the plane-wave fitter, doc/Records/VOneOfN_D2Design.md §2a row 7.)
    return blazem::NarrowExact<U>(GSpacePairProjectionIntegrator<dcmplx>(bft->Repulsion3C(*itsFitBasis))
                                  .Adjoint(RefreshedField(id)));
}
chmat_t Vee_Hartree::MakeMatrix(const cobs_t* bs, const Spin& s, DensityID id) const {return MakeMatrixT<dcmplx>(bs,s,id);}
rsmat_t Vee_Hartree::MakeMatrix(const robs_t* bs, const Spin& s, DensityID id) const {return MakeMatrixT<double>(bs,s,id);}

// V_H(dm) for this density, kept until the density's serial moves (declaration doc in PWTerms.C).
const ΔG_Map& Vee_Hartree::CoulombField(const cChargeDensity* cd) const
{
    assert(cd);
    if (!(DensityID(cd)==itsFieldID))
    {
        qchem::report::Timed miss("scf: E_H V_H field build (memo miss)");
        auto* fd=dynamic_cast<const qchem::ChargeDensity::FourierDensity*>(cd);
        assert(fd && "Vee_Hartree requires a FourierDensity (periodic) charge density");
        itsField       =fd->GetRepulsion3C(*itsFitBasis);
        itsFieldID=DensityID(cd);
    }
    return itsField;
}

// The block path's read (V1.41 stage 2, ruling C): one current field, filled by RefreshForDensity (or the energy
// pass), checked against the density's identity.  A miss is a caller that skipped the refresh, a broken contract
// no caller can act on, so it throws instead of building the field inside the block loop.
const ΔG_Map& Vee_Hartree::RefreshedField(DensityID id) const
{
    if (id.IsNone() || !(id==itsFieldID))
        throw std::logic_error("Vee_Hartree: block asked for a density the refresh did not see (RefreshForDensity first)");
    return itsField;
}

// Omega = Integral(1) through the raster's OWN quadrature rule -- so the constant comes from the same face
// (and the same summation order) every other raster integral in this file uses, and the term still holds no
// Structure.  Geometry-fixed, so it is asked once and memoized.
double Vee_Hartree::Volume() const
{
    if (itsVolume==0.0)
    {
        auto* rt=dynamic_cast<const BasisSet::G_RasterTransform*>(itsFitBasis.get());
        assert(rt && "Vee_Hartree: the CD fit basis is not a raster (G_RasterTransform) -- no volume to ask for");
        itsVolume=rt->Integral(rvec_t(rt->RasterSize(),1.0));
    }
    return itsVolume;
}

// E_H = 1/2 integral rho V_H[rho] -- the 1/2 is the electron-electron double-counting factor.
//
// ★ IN G SPACE, AND THAT IS THE WHOLE POINT (2026-09-05).  The obvious form is 0.5*Tr(D <i|V_H|j>), and it
// is what this did -- but <i|V_H|j> is a full real-space GATHER, and the ENERGY is evaluated at rho_new
// while the FOCK was built at rho_mix, so the two are different fields and no memo can join them.  Measured
// on the MnO parity row: the Hartree matrix was built TWICE per SCF iteration, 2 of the 4 gathers, against
// CP2K's 2 integrate_v_rspace calls per step for the WHOLE KS matrix (doc/Benchmark.md §5f).
//
// The pairing needs no matrix.  The gather is the EXACT ADJOINT of the collocation on this same fit grid,
// so by Parseval, with rho-tilde and V_H both given over the fit ball {G},
//     Tr(D <i|V_H|j>) = integral rho V_H = Omega * Sum_{dm} conj(rho-tilde(dm)) V_H(dm),
// and V_H is the density's OWN answer -- one collocation (a CollocMemo replay, since the Fock already
// collocated this D) plus an FFT: ~0.01 s where the gather was ~1.9 s.
// ★ AND ONE MAP IS ENOUGH.  V_H = k rho-tilde with the Poisson kernel k REAL, so
// conj(rho-tilde) V_H = |V_H|^2/k -- the field pairs with its own source and rho-tilde is never fetched.
// That matters: a second fetch costs a second IBZ STAR-AVERAGE of the whole map (measured 16 ms/call on
// Si Gamma with 48 point ops -- five times the gather it was replacing there).  dm=0 drops out on its own
// (k=0, the neutralising background) and the sum is manifestly REAL, which the trace form was only by
// symmetry.
//
// ⚠ NOT BIT-IDENTICAL with the trace form: it is the same integral summed in a different order, so pinned
// energies move at roundoff scale.  GPW_EH_TRACE=1 computes BOTH and prints the difference -- the A/B that
// says so, kept as an instrument rather than deleted with the evidence.
// THE EAGER REFRESH PHASE (doc/OpenWork.md item KP).  V_H[rho] is a function of the density alone, so one
// evaluation serves every Bloch block of the iteration -- but it used to be computed inside whichever
// block's MakeMatrix ran first, which is a WRITE during what is otherwise a read-only per-block loop.
// CoulombField() keeps its own serial guard, so this is a PRE-WARM and not a replacement: a caller that
// reaches the term without a prologue still gets the right field, just built later.
void Vee_Hartree::RefreshForDensity(const cChargeDensity* cd) const
{
    assert(cd);
    Volume();              // geometry-fixed; memoised on first ask
    CoulombField(cd);      // the density-dependent half -- the one that matters here
}

void Vee_Hartree::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    newCD(cd);
    // The kernel is asked of the FIT BASIS (its G-space Poisson seam), so this term still holds no Structure.
    auto* pk=dynamic_cast<const BasisSet::G_PoissonKernel*>(itsFitBasis.get());
    assert(pk && "Vee_Hartree: the CD fit basis carries no Poisson kernel (G_PoissonKernel) face");
    qchem::report::Timed timer("scf: E_H G-space pairing (no matrix)");
    const ΔG_Map& VH=CoulombField(cd);                     // 4pi rho-tilde/|G|^2 (Poisson, kernel baked)
    double e=0.0;
    for (const auto& [dm,v] : VH)
        if (const double k=pk->CoulombKernel(dm); k>0.0) e+=std::norm(v)/k;   // conj(rho-tilde) V_H = |V_H|^2/k
    e*=0.5*Volume();
    if (qchem::Diagnostics::Enabled("eh_trace"))
    {
        const double eTrace=0.5*cd->DM_Contract(this,cd);
        std::cout << "[E_H A/B] G-space=" << std::setprecision(12) << e << " trace=" << eTrace
                  << " diff=" << e-eTrace << " rel=" << (eTrace!=0.0 ? (e-eTrace)/eTrace : 0.0)
                  << std::setprecision(6) << std::endl;
    }
    te.Add("Eee", e, EnergyRole::Potential);   // Tr(D V_H[rho_mix]) under mixing is NOT 2e: not claimed (see FittedVee)
}

// THE RESPONSE: the transition density's own G+q Poisson -- linear, so no D0 -- then each pair's exact adjoint.
void Vee_Hartree::RefreshForDensity(const cbs_t*, const cChargeDensity*, const TransitionDensity<dcmplx>& delta) const
{
    if (delta.Version()==itsDeltaVersion) return;
    auto* fd=dynamic_cast<const qchem::ChargeDensity::TransitionFourierDensity*>(&delta);   // abstract -> abstract
    if (!fd) throw std::logic_error("Vee_Hartree: this transition density has no G+q (TransitionFourierDensity) face");
    qchem::report::Timed miss("response: delta V_H field");
    itsDeltaField  =fd->GetTransitionRepulsion(*itsFitBasis);
    itsDeltaVersion=delta.Version();
}

mat_t<dcmplx> Vee_Hartree::GetMatrix(const cobs_t* bra, const cobs_t* ket, const Spin&, const TransitionDensity<dcmplx>& delta) const
{
    if (delta.Version()!=itsDeltaVersion) RefreshForDensity(nullptr, nullptr, delta);
    auto* k=dynamic_cast<const BasisSet::Transition_DFT_IBS*>(ket);   // abstract -> abstract
    auto* b=dynamic_cast<const BasisSet::Transition_DFT_IBS*>(bra);
    if (!k || !b) throw std::logic_error("Vee_Hartree: the response needs Transition_DFT_IBS (collocating) blocks");
    return k->TransitionPotential(*itsFitBasis, *b, itsDeltaField);
}

std::ostream& Vee_Hartree::Write(std::ostream& os) const
{
    return os << "    PW electron-electron: Hartree V_H[rho] (G-space Poisson)." << std::endl;
}


} //namespace
