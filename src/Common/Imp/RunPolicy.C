// File: Common/Imp/RunPolicy.C  Resolving the deviation set once (see the interface for WHY).
//
// THE SITES THAT CONSULT THIS POLICY -- kept here so the list is checkable against a grep:
//   DMLowRank   -> ChargeDensity/Imp/Factory.C          (which rho route the density factory builds)
//   StreamFold  -> BasisSet/Lattice_3D/Imp/BasisSet.C   (whether the collocation streams are orbit-folded)
//   MixRhoM     -> ChargeDensity/DensityMixer.C         (which channel basis MakePeriodicMixer mixes in)
//   XCFromDM    -> Hamiltonian/Internal/Imp/PWTerms.C   (which rho the XC term is fed)
//   SymmetryImposition -> Calculation/Imp/SolidCalculation.C (ANDed with SolidCalcOptions::imposeSymmetry)
//   BeckeXC     -> Calculation/Imp/SolidCalculation.C          (passed to qcMesh::ResolveXCMesh as allowBecke)
//   DAwareScreen-> BasisSet/Lattice_3D/Evaluators/GPW/Imp/Evaluator.C (WHICH LatticeScreener the evaluator builds)
// Two further CP2K deviations are TYPED OPTIONS rather than env flags and are therefore NOT here:
// SolidCalcOptions::raster (BallOnly -- which IS CP2K's bet, vindicated by doc/OpenWork.md N2) and
// SolidCalcOptions::cutoffFactor (C=2).  They are chosen by the caller and reported by the run banner
// beside these, so the printed table is still complete even though the mechanism differs.
// (xcMesh WAS a third such typed option; 2026-08-28 promoted it to the table above -- it is 43% of the
// MnO row, which is far too large a difference to leave resting on caller discipline.)
module;
#include <cstdlib>
#include <string>
#include <sstream>
#include <vector>
module qchem.RunPolicy;

namespace qchem
{

// A knob is SET when the variable exists at all; its VALUE is the usual 0/non-0.  Distinguishing
// "set to 0" from "not set" is the whole point -- it is what lets an explicit knob outrank CP2K_COMPAT.
static bool IsSet (const char* n) {return std::getenv(n)!=nullptr;}
static bool AsBool(const char* n) {const char* s=std::getenv(n); return s && std::atoi(s)!=0;}

RunPolicy::RunPolicy()
    : itsCP2KCompat(AsBool("CP2K_COMPAT"))
{
    itsDMLowRank  = Resolve("QCHEM_DM_LOWRANK", "factored/low-rank rho route",
                            /*cp2k*/false, /*qchem default*/true);
    itsStreamFold = Resolve("GPW_STREAM_FOLD",  "orbit fold on the GPW collocation pair streams",
                            /*cp2k*/false, /*qchem default*/true);
    itsMixRhoM    = Resolve("QCHEM_MIX_RHO_M",  "(rho,m) mixing channels instead of (up,dn)",
                            /*cp2k*/false, /*qchem default*/false);
    itsXCFromDM   = Resolve("GPW_XC_DM_SOURCE", "Vxc fed rho[D] wholesale instead of rho_mix",
                            /*cp2k*/false, /*qchem default*/false);
    // NB the qchem default here is TRUE meaning "obey the caller", not "impose": the option itself
    // defaults off in SolidCalcOptions.  What CP2K parity forbids is the CAPABILITY, so that is what is
    // tabled -- and the facade ANDs this with the caller's own flag.
    itsImpose     = Resolve("QCHEM_IMPOSE_SYMMETRY", "space-group imposition available to the caller",
                            /*cp2k*/false, /*qchem default*/true);
    itsBeckeXC    = Resolve("QCHEM_BECKE_XC", "atom-centred (Becke) XC quadrature instead of the uniform grid",
                            /*cp2k*/false, /*qchem default*/true);
    // The knob NAME is unchanged from the experiment it grew out of, so every measurement banked against
    // GPW_DAWARE_SCREEN=0 still reproduces -- but it now selects a screener OBJECT, not a branch.
    itsDAware     = Resolve("GPW_DAWARE_SCREEN", "D-aware collocation box tolerance eps/|c_ij| instead of flat eps",
                            /*cp2k*/false, /*qchem default*/true);
}

// EXPLICIT BEATS THE UMBRELLA (see the interface): if the knob was named at all, that is the answer,
// and `stated` records it so the banner can say the umbrella did not get its way.
Deviation RunPolicy::Resolve(const char* knob, const char* what, bool cp2kValue, bool qchemDefault)
{
    const bool stated=IsSet(knob);
    const bool value = stated       ? AsBool(knob)
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
    os<<"CP2K_COMPAT="<<(itsCP2KCompat?"1":"0")<<" -> "<<(AtParity()?"AT PARITY":"DEVIATING")<<";";
    for (const Deviation& d : Deviations())
        os<<"  "<<d.knob<<"="<<(d.value?"on":"off")<<(d.Deviates()?"*":"")<<(d.stated?"(stated)":"");
    os<<"   [* = differs from CP2K]";
    return os.str();
}

// ONE object for the process's lifetime, ASSIGNED (never replaced) by ReresolveRunPolicy, so a
// reference handed out earlier stays valid across an A/B flip.
static RunPolicy& thePolicy()
{
    static RunPolicy p;   // first use, after main() has had its chance to set the environment
    return p;
}
const RunPolicy& theRunPolicy()  { return thePolicy(); }
void ReresolveRunPolicy()        { thePolicy() = RunPolicy(); }

} //namespace
