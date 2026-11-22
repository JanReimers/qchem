// File: Calculation/XCRhoPolicy.C  The XC ρ-axis policy (user ruling 2026-11-09, V-XC-ONE-TERM): which ρ XC sees, ρ[D] or the Vxc fit's
// projection Pρ, and what the policy layer says about the combination stated.
//
// THE ρ AXIS IS RETIRED (user 2026-11-15, "Your plan sounds good"): XC sees ρ[D] (doc/Pins.md pin 28), so xc.rho fit is accepted only where
// the fit's P is EXACT (Pρ = ρ[D]) and REFUSED, citing pin 28, wherever it is not.  The SCF enforces the same rule at run time, fail closed
// (SCFIterator RequirePin28 on the XC term's certificate).  The cells:
//   - A1 x ρ[D]  the Gaussian x δ cell: v(ρ[D]) is FITTED for H, so H is not dE/dD (doc/Pins.md pin 32).  Ruling (b): it RUNS WITH A
//                WARNING -- an experienced user's cheap better-than-SAD start for a more accurate run.  (XC SEES ρ[D]: pin 28 holds.)
//   - A1 / multipole x Pρ   REFUSED (pin 28): the fit is incomplete, so Pρ is not ρ[D].
//   - δ  x Pρ    accepted: P is exact on the XC points, so Pρ = ρ[D] -- the row says so.
//   - PW x Pρ    exact over plane-wave orbitals (the 4 Ecut ball); over GPW orbitals REFUSED (pin 28), on the raster and at Becke points.
// A policy of the policy layer: each decision files one `policy` row where it fires (the caller's report cursor).
module;
#include <string>
export module qchem.XCRhoPolicy;
import qchem.Hamiltonian.Types;   // VxcFit, XCρ
import qchem.Hamiltonian.Factory; // XCFunctional, ConsumesDensityGradient, Describe (the GGA route refusal)
import qchem.Policy;

export namespace qchem
{

//! The pin-28 refusal of a stated xc.rho fit on an inexact fit, \a which naming the fit.
[[noreturn]] inline void RefuseProjectedRho(const std::string& which, const std::string& why)
{
    policy::Refuse("XCRhoProjectedRefused", "xc.rho fit on "+which,
        why+", so its projection P rho is not rho[D], and XC must see rho[D] (doc/Pins.md pin 28; the rho axis was retired 2026-11-15) -- "
        "state xc.rho exact (the default), or a fit whose projection is exact: xc.fit Delta");
}

//! File the ρ-axis decision for a DFT run stating \a fit and \a rho.  Called by the setups before the Hamiltonian is built (a solid
//! passes the fit its Hamiltonian RESOLVES from Auto, so the row names the cell that runs).  REFUSES xc.rho fit on an inexact fit (pin 28).
inline void ReportXCRho(Hamiltonian::VxcFit fit, Hamiltonian::XCρ rho)
{
    using Hamiltonian::VxcFit; using Hamiltonian::XCρ;
    if (fit==VxcFit::Auxiliary && rho==XCρ::Exact)
        policy::Report("XCRhoFittedPotential", "the Gaussian auxiliary Vxc fit at xc.rho exact: v(rho[D]) is fitted for H", policy::Verdict::Warn,
            "H is then not dE/dD (doc/Pins.md pin 32): the SCF stops off the stationary point of any energy, so its energy is not "
            "variational and nothing that needs dE/dD (forces, a variational restart claim) may be read from it.  It is meant as a cheap "
            "better-than-SAD start; state the default delta fit for H = dE/dD");
    else if (fit==VxcFit::Auxiliary && rho==XCρ::Fit)
        RefuseProjectedRho("the Gaussian auxiliary (A1) Vxc fit", "the A1 basis is incomplete");
    else if (fit==VxcFit::Multipole && rho==XCρ::Fit)
        RefuseProjectedRho("the multipole Vxc fit", "the site multipole channels are an incomplete fit");
    else if (rho==XCρ::Fit && (fit==VxcFit::Delta || fit==VxcFit::Auto))
        policy::Report("XCRhoProjected", "xc.rho fit on the delta Vxc fit: fit is exact here, P rho = rho[D]", policy::Verdict::Choose,
            "the delta basis IS the XC quadrature's points, so its projection is the identity there");
}

//! The SOLID ρ-axis decision, on the fit the Hamiltonian RESOLVES from Auto and the XC points \a becke.  The plane-wave fit's cells split by
//! the ORBITAL basis:
//!   - plane-wave orbitals: the Vxc ball is >= 4 Ecut by construction (\c PlaneWave_IBS::CreateVxcFitBasisSet, pin 28), so it holds
//!     every orbital-pair difference and P rho = rho[D] EXACTLY, on its raster or synthesized at Becke points (I3).  Both values of
//!     xc.rho are the same run; a stated fit is reported as such.
//!   - Gaussian (GPW) orbitals: a Gaussian pair product is not band-limited, so no ball holds rho[D].
//!       xc.rho fit is REFUSED (pin 28, 2026-11-15) on the raster and at Becke points alike: P rho is the raster's rho[D] truncated to the
//!       ball.  On the raster the realized route COLLOCATES rho[D] (\c PairDensitySampler), so xc.rho exact is honoured; on a Becke mesh
//!       xc.rho exact is REFUSED by name (user ruling 2026-11-14, option (ii)): XC of rho[D] at the Becke points IS the delta cell.
//! Every other cell is the fit-agnostic decision above.
inline void ReportSolidXCRho(Hamiltonian::VxcFit fit, Hamiltonian::XCρ rho, bool planeWaveOrbitals, bool becke)
{
    using Hamiltonian::VxcFit; using Hamiltonian::XCρ;
    if (fit!=VxcFit::PlaneWave) return ReportXCRho(fit, rho);
    if (planeWaveOrbitals)
    {
        if (rho==XCρ::Fit)
            policy::Report("XCRhoProjected", "xc.rho fit on the plane-wave Vxc fit over plane-wave orbitals: fit is exact here, P rho = rho[D]",
                policy::Verdict::Choose, "the Vxc ball is at least 4 Ecut, so it holds every orbital-pair difference G_i - G_j (doc/Pins.md pin 28)");
        return;
    }
    if (rho==XCρ::Fit)
        RefuseProjectedRho("the plane-wave Vxc fit over Gaussian (GPW) orbitals",
            std::string("a Gaussian pair product is not band-limited, so no G-ball holds rho[D]")+(becke ? " (its synthesis at the Becke points "
            "is the raster's rho[D] truncated to the ball)" : ""));
    if (becke && rho==XCρ::Exact)
        policy::Refuse("XCRhoBallOnBeckeExact", "the plane-wave Vxc fit on a Becke XC mesh at xc.rho exact over Gaussian (GPW) orbitals",
            "no ball holds a Gaussian pair product, so XC of rho[D] at the Becke points is the delta cell -- state xc.fit Delta (that cell, "
            "exactly), or a uniform XC mesh (user ruling 2026-11-14)");
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

//! THE GGA ROUTE REFUSAL (GGA stage 1; GGA design note §5; user ruling 2026-11-20).  A functional that reads \f$\nabla\rho\f$ is built and
//! named, but no XC sampling route samples the gradient yet (stages 2-4: the molecular point route, the periodic point route, the raster),
//! so every route REFUSES it here, before anything is built -- never a silent LDA in its place.  \a where names the run kind.
inline void RefuseUnroutedGradient(const Hamiltonian::XCFunctional& xc, const std::string& where)
{
    if (!Hamiltonian::ConsumesDensityGradient(xc)) return;
    policy::Refuse("GGARoute", Hamiltonian::Describe(xc)+" on "+where,
        "a GGA reads the density gradient, and no XC sampling route samples it yet (GGA stages 2-4: the molecular point route, the periodic "
        "point route, the GPW raster) -- state an LDA (model LDA or PADE) until the route lands");
}
//! \brief THE MOLECULAR GGA ROUTE (GGA stage 2): the all-electron δ Vxc fit samples \f$\nabla\rho\f$ (the density-gradient strategy), so a GGA
//! there RUNS and the policy says which route it took; a Gaussian (A1) or multipole Vxc fit has no gradient route and REFUSES, named.
//! (A pseudopotential molecule refuses any non-LDA model on its own, \c PseudopotentialIsLDA.)  \a model is the run's; a non-GGA passes silently.
inline void RouteMolecularGradient(Hamiltonian::Model model, Hamiltonian::VxcFit fit)
{
    if (!Hamiltonian::IsDFT(model) || !Hamiltonian::ConsumesDensityGradient(Hamiltonian::ModelToXC(model))) return;
    const std::string what=Hamiltonian::Describe(Hamiltonian::ModelToXC(model))+" on a molecule";
    if (fit!=Hamiltonian::VxcFit::Auto && fit!=Hamiltonian::VxcFit::Delta)
        policy::Refuse("GGARoute", what+" with a non-delta Vxc fit",
            "only the delta Vxc fit samples the density gradient on a molecule (GGA stage 2); the Gaussian auxiliary and multipole fits "
            "have no gradient route -- state the delta fit (the default)");
    policy::Report("GGARoute", what+": the delta Vxc fit samples grad rho at the XC points (GGA stage 2)", policy::Verdict::Choose,
        "the molecular point route: grad rho through the delta basis's gradient 3-centre overlap, H = dE/dD by the {v, w} gather");
}
//! The same for a stated MODEL (a non-DFT model reads no functional, so it passes).
inline void RefuseUnroutedGradient(Hamiltonian::Model m, const std::string& where)
{
    if (Hamiltonian::IsDFT(m)) RefuseUnroutedGradient(Hamiltonian::ModelToXC(m), where);
}

} // namespace qchem
