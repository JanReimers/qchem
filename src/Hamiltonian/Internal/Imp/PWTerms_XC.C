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
import qchem.Pseudopotential.Integrals_Pseudo;   // cast bs ACROSS to the external-PP operator-assembly mixin (Ven_PP_*)
import qchem.Fitting.FunctionFitter;        // Fitting::Factory (both PW fitters) + ProjectedDensity_G / ProjectedScalar_R
import qchem.Structure;                       // Structure::isFinite()/SumFormFactors() -- the G=0 alignment (term-side)
import qchem.Blaze;                            // blazem::zeroH<dcmplx> (the null-PP V_long block)
import qchem.Mesh.Quadrature;                 // qcMesh::Mesh (the Vxc_Quadrature engine's quadrature mesh)
import qchem.Reporting;                       // Timed (the setup/scf timing ledger)
import qchem.Parallel;                         // WorkerThreads (GPW_OMP_THREADS -- the XC-mesh table + quadrature loops)


namespace qchem::Hamiltonian
{

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
// T2 (doc/OpenWork.md N1): forward to the quadrature, which owns the atom-centred partition.  Empty when
// it has none -- the caller (SolidCalculation) treats empty as "this run cannot answer", never as "zero".
rvec_t Vxc_QuadraturePol::SiteMoments(const cChargeDensity* cd) const {return itsQuad->SiteMoments(cd);}
rvec_t Vcorr_QuadraturePol::SiteMoments(const cChargeDensity* cd) const {return itsQuad->SiteMoments(cd);}

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
    // The PER-VOLUME energy density, because that is the form that COMPOSES: exchange contributes
    // Σ_σ ε_x(ρ_σ)ρ_σ and correlation ε_c·ρ_tot, and those share no denominator (see GetExcDensity).
    // For a plain correlation functional the default IS ε_c·(ρ↑+ρ↓), so this line is bit-identical to the
    // one it replaces; for the composite it is the only correct sum.
    for (size_t g=0; g<up.size(); g++) ec[g]=itsCorr->GetExcDensity(up[g], dn[g]);
    te.Exc += itsQuad->Integrate(ec);   // E_xc = ∫ e_xc(ρ↑,ρ↓)
}

std::ostream& Vcorr_QuadraturePol::Write(std::ostream& os) const
{
    // It carries the WHOLE spin-native functional since 2026-09-04 (MakeVxcTerms hands it a
    // CompositeExFunctional summing exchange AND correlation into one gather), so the line must not still
    // say "correlation" -- the console is how a run states what it built.
    return os << "    XC-mesh SPIN-NATIVE v_xc^sigma(rho_up,rho_down), exchange+correlation in ONE gather ("
              << itsQuad->NumPoints() << " atom-centred points)." << std::endl;
}


// ---- THE EAGER REFRESH PHASE (doc/OpenWork.md item KP) -------------------------------------------------
// rho on the quadrature's points is k-INDEPENDENT -- one array, correct for every Bloch block -- but it was
// filled by whichever block asked first.  Warming it here makes the block loop read-only in the ordinary
// path, which is the precondition for running the blocks concurrently.  All three delegate to the shared
// engine, so an XC PAIR (exchange + correlation over one quadrature) warms exactly once between them: the
// engine's serial guard turns the second call into a lookup.
void Vxc_Quadrature      ::RefreshForDensity(const cChargeDensity* cd) const {itsQuad->WarmForDensity(cd,false);}
void Vxc_QuadraturePol   ::RefreshForDensity(const cChargeDensity* cd) const {itsQuad->WarmForDensity(cd,true );}
void Vcorr_QuadraturePol ::RefreshForDensity(const cChargeDensity* cd) const {itsQuad->WarmForDensity(cd,true );}

} //namespace
