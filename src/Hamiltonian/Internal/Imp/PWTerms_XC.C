// File: Hamiltonian/Internal/Imp/PWTerms_XC.C  the XC TERMS -- the thin layer that owns the FUNCTIONAL and nothing else.
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
import qchem.RunPolicy;   // theRunPolicy().XCFromDM() -- the declared XC-feed deviation (N5)
import qchem.Energy;
import qchem.ChargeDensity;
import qchem.ChargeDensity.FourierDensity;   // cast cd UP to its reciprocal-space coefficients rho-tilde
import qchem.BasisSet.Orbital_DFT_IBS;         // cast bs UP to the reciprocal-space DFT capability (Hartree/XC)
import qchem.BasisSet.G_FieldEvaluator;    // G_RasterTransform: the fit basis's FFT pair (RhoOnGrid, the BALL route)
import qchem.BasisSet.Orbital_PP_IBS;             // cast bs ACROSS to the species-field integral service (Ven_PP_*)
import qchem.Fitting.FunctionFitter;        // Fitting::Factory (both PW fitters) + ProjectedDensity_G / ProjectedScalar_R
import qchem.Structure;                       // Structure::isFinite()/SumFormFactors() -- the G=0 alignment (term-side)
import qchem.Blaze;                            // blazem::zeroH<dcmplx> (the null-PP V_long block)
import qchem.Mesh.Quadrature;                 // qcMesh::Mesh (the Vxc_Quadrature engine's quadrature mesh)
import qchem.Reporting;                       // Timed (the setup/scf timing ledger)
import qchem.Parallel;                         // WorkerThreads (GPW_OMP_THREADS -- the XC-mesh table + quadrature loops)


namespace qchem::Hamiltonian
{

// ---- Vxc_Quadrature ------------------------------------------------------------------------------------------

// Built with the SHARED quadrature engine (mesh + Phi tables + per-serial rho) and the imposed spin subgroup.
Vxc_Quadrature::Vxc_Quadrature(const xc_t& xc, sampler_t quad, SpinGroup g)
    : itsXc(xc)
    , itsSampler(std::move(quad))
    , itsGroup(g)
{
    assert(itsXc);
    assert(itsSampler);
}

// THE CHANNEL RASTERS -- the one place the imposed subgroup is read on the Fock/energy path.  Polarized:
// the sampler's {up,dn} pair (a spin-agnostic seed collapses to rho/2 per channel INSIDE RhoPol).  SU(2):
// the folded doublet's one raster, halved into `scratch` and handed over as BOTH channels -- the exact
// zeta=0 collapse, so the functional's spin-native face reproduces the old scalar path bit for bit
// (2*(rho/2)==rho; VWN routes exact zeta=0 through its scalar formula).
Vxc_Quadrature::Rasters Vxc_Quadrature::ChannelRasters(const cChargeDensity* cd, rvec_t& scratch) const
{
    if (itsGroup==SpinGroup::Polarized)
        return {itsSampler->RhoPol(cd,Spin::Up), itsSampler->RhoPol(cd,Spin::Down)};
    scratch=0.5*itsSampler->Rho(cd);
    return {scratch, scratch};
}

// v^sigma(rho_up, rho_dn) pointwise on the channel rasters, then the engine's Phi-table quadrature (one
// GEMM).  ONE body for both block scalars (Step 3c): the rasters are block-independent; only the final
// quadrature (typed Phi table) differs.  A folded-doublet block (Spin::None) is the same field either
// channel would see, so it asks for Up by convention.
template <class U> hmat_t<U> Vxc_Quadrature::MakeMatrixT(const tobs_t<U>* bs, const Spin& s, const cChargeDensity* cd) const
{
    rvec_t half;
    const Rasters r=ChannelRasters(cd, half);
    const Spin fs = s==Spin::None ? Spin::Up : s;
    rvec_t v(r.up.size());
    for (size_t g=0; g<v.size(); g++) v[g]=itsXc->GetVxc(r.up[g], r.dn[g], fs);
    return itsSampler->Matrix(bs, v);
}
chmat_t Vxc_Quadrature::MakeMatrix (const cobs_t* bs, const Spin& s, const cChargeDensity* cd) const {return MakeMatrixT<dcmplx>(bs,s,cd);}
rsmat_t Vxc_Quadrature::MakeMatrixR(const robs_t* bs, const Spin& s, const cChargeDensity* cd) const {return MakeMatrixT<double>(bs,s,cd);}

void Vxc_Quadrature::GetEnergy(EnergyBreakdown& te, const cDM_CD* cd) const
{
    rvec_t half;
    const Rasters r=ChannelRasters(cd, half);   // reuses the iteration's rasters (same density serial)
    rvec_t exc(r.up.size()), rho(r.up.size());
    // E_xc = Integral e(rho_up,rho_dn), the PER-VOLUME energy density, because that is the form that
    // COMPOSES: exchange contributes Sum_sigma eps_x(rho_sigma) rho_sigma and correlation eps_c rho_tot, and
    // those share no denominator.
    for (size_t g=0; g<exc.size(); g++) { exc[g]=itsXc->GetExcDensity(r.up[g], r.dn[g]); rho[g]=r.up[g]+r.dn[g]; }
    te.Add("Exc", itsSampler->Integrate(exc), EnergyRole::Potential);
    // The charge accounting of the same density: the mesh-charge leak (the quadrature's health metric --
    // CP2K's grid-charge-lost readout) is the quadrature integral of rho vs the analytic Tr(DS); the
    // integrated site moments are the observable, where both rasters are in hand (R1.0h).
    te.charge.lost        = itsSampler->Integrate(rho) - cd->GetTotalCharge();
    te.charge.siteMoments = SiteMoments(cd);
}

// THE OBSERVABLE (doc/OpenWork.md N1/T2; R1.0h): mu_A = Integral w_A (rho_up - rho_dn), in electrons, over the
// sampler's atom-centred partition.  FREE -- both channel rasters are cached for this density serial -- and
// the term's, not the sampler's, because the term is what knows the field is a spin difference; the
// sampler only knows how to partition-integrate a field at its points.  Empty when the run resolves no
// spin (m == 0 identically) or the quadrature has no site blocks -- the caller treats empty as "this run
// cannot answer", never as "zero".
rvec_t Vxc_Quadrature::SiteMoments(const cChargeDensity* cd) const
{
    if (itsGroup!=SpinGroup::Polarized) return rvec_t();
    return itsSampler->SiteIntegrals(rvec_t(itsSampler->RhoPol(cd,Spin::Up)-itsSampler->RhoPol(cd,Spin::Down)));
}

std::ostream& Vxc_Quadrature::Write(std::ostream& os) const
{
    return os << "    XC-mesh " << (itsGroup==SpinGroup::Polarized ? "SPIN-NATIVE v_xc^sigma(rho_up,rho_down)"
                                                                   : "v_xc(rho(r))")
              << ", exchange+correlation in ONE gather (" << itsSampler->NumPoints() << " points)." << std::endl;
}

// ---- THE EAGER REFRESH PHASE (doc/OpenWork.md item KP) -------------------------------------------------
// rho on the quadrature's points is k-INDEPENDENT -- one array (or one pair), correct for every Bloch block
// -- but it was filled by whichever block asked first.  Warming it here makes the block loop read-only in
// the ordinary path, which is the precondition for running the blocks concurrently.
void Vxc_Quadrature::RefreshForDensity(const cChargeDensity* cd) const
{
    itsSampler->WarmForDensity(cd, itsGroup==SpinGroup::Polarized);
}

} //namespace
