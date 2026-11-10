// File: Calculation/XCRhoPolicy.C  The XC ρ-axis policy (user ruling 2026-11-09, V-XC-ONE-TERM): which ρ XC sees, ρ[D] or the Vxc fit's
// projection Pρ, and what the policy layer says about the combination stated.
//
// The two values are one axis (Hamiltonian::XCρ).  They coincide where the fit's P is exact (the δ basis on the XC points), and split
// for the Gaussian auxiliary (A1) fit:
//   - A1 x ρ[D]  the Gaussian x δ cell: v(ρ[D]) is FITTED for H, so H is not dE/dD (doc/Pins.md pin 32) and the SCF stops off the
//                stationary point of any energy.  Ruling (b): it RUNS WITH A WARNING -- an experienced user's cheap better-than-SAD
//                start for a more accurate run.
//   - A1 x Pρ    variational (E = ∫ e(Pρ), H = dE/dD), but XC sees the fitted density: an explicit departure from pin 28, reported.
//   - δ  x Pρ    accepted: P is exact on the XC points, so Pρ = ρ[D] -- the row says so.
// A policy of the policy layer: each decision files one `policy` row where it fires (the caller's report cursor).
module;
#include <string>
export module qchem.XCRhoPolicy;
import qchem.Hamiltonian.Types;   // VxcFit, XCρ
import qchem.Policy;

export namespace qchem
{

//! File the ρ-axis decision for a DFT run stating \a fit and \a rho.  Called by the setups before the Hamiltonian is built (a solid
//! passes the fit its Hamiltonian RESOLVES from Auto, so the row names the cell that runs); reports, never refuses (the cells it sees
//! all run).
inline void ReportXCRho(Hamiltonian::VxcFit fit, Hamiltonian::XCρ rho)
{
    using Hamiltonian::VxcFit; using Hamiltonian::XCρ;
    if (fit==VxcFit::Auxiliary && rho==XCρ::Exact)
        policy::Report("XCRhoFittedPotential", "the Gaussian auxiliary Vxc fit at xc.rho exact: v(rho[D]) is fitted for H", policy::Verdict::Warn,
            "H is then not dE/dD (doc/Pins.md pin 32): the SCF stops off the stationary point of any energy, so its energy is not "
            "variational and nothing that needs dE/dD (forces, a variational restart claim) may be read from it.  It is meant as a cheap "
            "better-than-SAD start; state xc.rho fit for the variational A1 route, or the default delta fit for rho[D] exactly");
    else if (fit==VxcFit::Auxiliary && rho==XCρ::Fit)
        policy::Report("XCRhoProjected", "XC sees P rho, the Gaussian auxiliary fit's projection of rho[D]", policy::Verdict::Choose,
            "stated (xc.rho fit): variational (E = integral e(P rho), H = dE/dD), but the A1 basis is incomplete, so XC is not evaluated on "
            "rho[D] (an explicit departure from doc/Pins.md pin 28)");
    else if (rho==XCρ::Fit && (fit==VxcFit::Delta || fit==VxcFit::Auto))
        policy::Report("XCRhoProjected", "xc.rho fit on the delta Vxc fit: fit is exact here, P rho = rho[D]", policy::Verdict::Choose,
            "the delta basis IS the XC quadrature's points, so its projection is the identity there");
    else if (rho==XCρ::Fit && fit==VxcFit::PlaneWave)
        policy::Report("XCRhoProjected", "xc.rho fit on the plane-wave Vxc fit: XC sees the G-ball's projection of rho[D]", policy::Verdict::Choose,
            "stated: the plane-wave fit is the density's Fourier series on its ball, so P rho = rho[D] wherever the ball holds every "
            "orbital-pair wave vector");
}

//! \brief REFUSE what needs H = dE/dD on the Gaussian x δ cell (ruling (b): "everything that needs H = dE/dD refuses on it by name").
//! \a accelerator is the stated accelerator type ("DIIS" | "GDM" | "Ladder"): GDM, and a Ladder that hands off to its GDM rung, descend
//! E along a geodesic whose direction IS H, so on this cell their line search would minimise a different function from the one whose
//! gradient they follow.  (Forces refuse on their own: the fitted-v term has no geometry pullback, and every such term names itself.)
inline void RefuseGradientConsumers(Hamiltonian::VxcFit fit, Hamiltonian::XCρ rho, const std::string& accelerator)
{
    using Hamiltonian::VxcFit; using Hamiltonian::XCρ;
    if (fit==VxcFit::Auxiliary && rho==XCρ::Exact && (accelerator=="GDM" || accelerator=="Ladder"))
        policy::Refuse("XCRhoFittedPotentialGradient", "the "+accelerator+" accelerator with the Gaussian auxiliary Vxc fit at xc.rho exact",
            "GDM descends E along the direction H, and on this cell H is the fitted v(rho[D]), not dE/dD (doc/Pins.md pin 32), so its line "
            "search would minimise a different function from the one it follows -- state xc.rho fit (the variational A1 route) or the DIIS "
            "accelerator");
}

} // namespace qchem
