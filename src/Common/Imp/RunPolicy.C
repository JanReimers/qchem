// File: Common/Imp/RunPolicy.C  Resolving the deviation set once (see the interface for WHY).
//
// THE SITES THAT CONSULT THIS POLICY -- kept here so the list is checkable against a grep:
//   DMLowRank   -> ChargeDensity/Imp/Factory.C          (which rho route the density factory builds)
//   StreamFold  -> BasisSet/Lattice/Imp/BasisSet.C   (whether the collocation streams are orbit-folded)
//   MixRhoM     -> ChargeDensity/Imp/DensityMixer.C     (which channel basis the G-space factories compose in)
//   XCFromDM    -> Hamiltonian/Internal/Imp/PWTerms.C   (which rho the XC term is fed)
//   SymmetryImposition -> Calculation/Imp/SolidCalculation.C (ANDed with SolidCalcOptions::imposeSymmetry)
//   BeckeXC     -> Calculation/Imp/SolidCalculation.C          (passed to qcMesh::ResolveXCMesh as allowBecke)
//   DAwareScreen-> BasisSet/Lattice/Evaluators/GPW/Imp/Evaluator.C (WHICH LatticeScreener the evaluator builds)
// Two further CP2K deviations are TYPED OPTIONS rather than env flags and are therefore NOT here:
// SolidCalcOptions::raster (BallOnly -- which IS CP2K's bet, vindicated by doc/OpenWork.md N2) and
// SolidCalcOptions::cutoffFactor (C=2).  They are chosen by the caller and reported by the run banner
// beside these, so the printed table is still complete even though the mechanism differs.
// (xcMesh WAS a third such typed option; 2026-08-28 promoted it to the table above -- it is 43% of the
// MnO row, which is far too large a difference to leave resting on caller discipline.)
module;
#include <optional>
#include <string>
#include <sstream>
#include <vector>
module qchem.RunPolicy;

namespace qchem
{

RunPolicy::RunPolicy(const RunPolicySpec& spec)
    : itsSpec(spec), itsCP2KCompat(spec.cp2kCompat)
{
    itsDMLowRank  = Resolve("policy.dmLowRank", "factored/low-rank rho route",
                            /*cp2k*/false, /*qchem default*/true, spec.dmLowRank);
    itsStreamFold = Resolve("policy.streamFold",  "orbit fold on the GPW collocation pair streams",
                            /*cp2k*/false, /*qchem default*/true, spec.streamFold);
    // N3 PROMOTED 2026-09-20 (OpenWork section 4): measured on MnO AFM-II with the deck's loop shape, (rho,m)
    // takes 18 / 21 iterations where (up,dn) takes 22 / 24 (U=0 / U=4 eV), energies identical to 1e-10 Ha --
    // Kerker's 4pi/G^2 has no business damping the magnetisation channel.  Still off under CP2K_COMPAT.
    itsMixRhoM    = Resolve("policy.mixRhoM",  "(rho,m) mixing channels instead of (up,dn)",
                            /*cp2k*/false, /*qchem default*/true, spec.mixRhoM);
    itsXCFromDM   = Resolve("policy.xcFromDM", "Vxc fed rho[D] wholesale instead of rho_mix",
                            /*cp2k*/false, /*qchem default*/false, spec.xcFromDM);
    // NB the qchem default here is TRUE meaning "obey the caller", not "impose": the option itself
    // defaults off in SolidCalcOptions.  What CP2K parity forbids is the CAPABILITY, so that is what is
    // tabled -- and the facade ANDs this with the caller's own flag.
    itsImpose     = Resolve("policy.imposeSymmetry", "space-group imposition available to the caller",
                            /*cp2k*/false, /*qchem default*/true, spec.imposeSymmetry);
    itsBeckeXC    = Resolve("policy.beckeXC", "atom-centred (Becke) XC quadrature instead of the uniform grid",
                            /*cp2k*/false, /*qchem default*/true, spec.beckeXC);
    // (It selects a screener OBJECT, not a branch.)
    itsDAware     = Resolve("policy.dAwareScreen", "D-aware collocation box tolerance eps/|c_ij| instead of flat eps",
                            /*cp2k*/false, /*qchem default*/true, spec.dAwareScreen);
    itsUEigen     = Resolve("policy.hubbardEigen", "DFT+U on the Lowdin block's eigenvalues (Dudarev) instead of its diagonal populations",
                            /*cp2k*/false, /*qchem default*/true, spec.hubbardEigen);
}

// EXPLICIT BEATS THE UMBRELLA (see the interface): if the knob was named at all, that is the answer,
// and `stated` records it so the banner can say the umbrella did not get its way.
Deviation RunPolicy::Resolve(const char* knob, const char* what, bool cp2kValue, bool qchemDefault, const std::optional<bool>& statedValue)
{
    const bool stated=statedValue.has_value();
    const bool value = stated       ? *statedValue
                     : itsCP2KCompat ? cp2kValue
                     :                 qchemDefault;
    return Deviation{knob, what, cp2kValue, value, stated};
}

bool RunPolicy::AtParity() const
{
    for (const Deviation& d : Deviations()) if (d.Deviates()) return false;
    return true;
}

std::string RunPolicy::Banner() const
{
    std::ostringstream os;
    os<<"policy.cp2kCompat="<<(itsCP2KCompat?"1":"0")<<" -> "<<(AtParity()?"AT PARITY":"DEVIATING")<<";";
    for (const Deviation& d : Deviations())
        os<<"  "<<d.knob<<"="<<(d.value?"on":"off")<<(d.Deviates()?"*":"")<<(d.stated?"(stated)":"");
    os<<"   [* = differs from CP2K]";
    return os.str();
}

// ONE object for the process's lifetime, ASSIGNED (never replaced) by ReresolveRunPolicy, so a
// reference handed out earlier stays valid across an A/B flip.
static RunPolicy& thePolicy()
{
    static RunPolicy p;   // the DEFAULT policy until a facade installs a stated one
    return p;
}
const RunPolicy& theRunPolicy()  { return thePolicy(); }
void SetRunPolicy(const RunPolicySpec& spec) { thePolicy() = RunPolicy(spec); }

} //namespace
