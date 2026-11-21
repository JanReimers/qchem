// File: Calculation/Imp/SetupDecisions.C  The bodies of the structure-neutral setup decisions (interface and docs: ../SetupDecisions.C).
module;
#include <optional>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
module qchem.SetupDecisions;
import qchem.Policy;                         // policy::Refuse / RefuseAs: each refusal files its row

namespace qchem
{

SpinChannels ResolveSpinChannels(int multiplicity, int nel, std::optional<bool> spinsShareFermi, const std::string& who)
{
    const bool singlet = multiplicity==1 && !spinsShareFermi.value_or(false);
    if (singlet && nel%2!=0)
        policy::RefuseAs("SpinMultiplicity", who+": multiplicity 1 is the SINGLET (nUp = nDown), impossible with an odd Nelec "
                                            +std::to_string(nel)+" -- state multiplicity 2 (a doublet) or 0 (minimal spin)"
                                            +(spinsShareFermi ? ", or spinsShareFermi true so one Fermi level sets the channels (a metal)" : ""));
    const int twoS = multiplicity>1 ? multiplicity-1 : (singlet ? 0 : nel%2);
    if (twoS>nel || (nel-twoS)%2!=0)
        policy::RefuseAs("SpinMultiplicity", who+": multiplicity "+std::to_string(multiplicity)+" (2S="+std::to_string(twoS)
                                            +") parity disagrees with Nelec "+std::to_string(nel)+" -- need 0 <= 2S <= Nelec and Nelec-2S even");
    return {(nel+twoS)/2, (nel-twoS)/2};
}

SpinGroup ResolveSpinGroup(SpinGroup stated, int multiplicity)
{
    return multiplicity>1 ? SpinGroup::Polarized : stated;
}

SCFAccelerators::AcceleratorParams ResolveAcceleratorParams(const AcceleratorOptions& acc, SCFAccelerators::AcceleratorParams p)
{
    auto gdm=[](const AcceleratorOptions::GDM& g, SCFAccelerators::GDMParams& o)
    {
        if (g.engageBelowFD) o.EngageBelowFD=*g.engageBelowFD;
        if (g.trust)         o.Trust=*g.trust;
        if (g.trustBackoff)  o.TrustBackoff=*g.trustBackoff;
        if (g.trustMin)      o.TrustMin=*g.trustMin;
        if (g.pcFloor)       o.PCFloor=*g.pcFloor;
        if (g.aufbauTol)     o.AufbauTol=*g.aufbauTol;
        if (g.aufbauMaxRepicks) o.AufbauMaxRepicks=*g.aufbauMaxRepicks;
    };
    if (acc.diis.nProj)         p.diis.Nproj=*acc.diis.nProj;
    if (acc.diis.engageBelowFD) p.diis.EngageBelowFD=*acc.diis.engageBelowFD;
    if (acc.diis.svTol)         p.diis.SVTol=*acc.diis.svTol;
    if (acc.diis.svTolRel)      p.diis.SVTolRel=*acc.diis.svTolRel;
    gdm(acc.gdm, p.gdm);
    gdm(acc.ladder.gdm, p.ladderGdm);
    if (acc.ladder.eThresh)  p.ladderEThresh=*acc.ladder.eThresh;
    if (acc.ladder.stall)    p.ladderStall=*acc.ladder.stall;
    if (acc.ladder.floor)    p.ladderFloor=*acc.ladder.floor;
    if (acc.ladder.switchAt) p.ladderSwitchAt=*acc.ladder.switchAt;
    if (acc.ladder.signal)   p.ladderSignal=*acc.ladder.signal;
    p.ladderStartOnGDM = p.ladderStartOnGDM || acc.ladder.startOnGDM;
    return p;
}

SCFAccelerators::AcceleratorParams MolecularAcceleratorDefaults(double eMax)
{
    SCFAccelerators::AcceleratorParams p;
    p.diis      = {.Nproj=4, .EngageBelowFD=eMax, .SVTol=0.0, .SVTolRel=1e-4};
    p.gdm       = {.EngageBelowFD=eMax, .Trust=0.1};
    p.ladderGdm = {.EngageBelowFD=1e10, .Trust=0.1};    // the rung always steps once it is the active one
    p.ladderEThresh=1e-8; p.ladderStall=5; p.ladderFloor=1e-8; p.ladderSwitchAt=1e-6;
    p.ladderSignal=SCFAccelerators::ScheduleSignal::EnergyChange;
    return p;
}

nlohmann::json ToAcceleratorJson(const SCFAccelerators::AcceleratorParams& p)
{
    return {{"NProj",p.diis.Nproj},{"EngageBelowFD",p.diis.EngageBelowFD},{"SVTol",p.diis.SVTol},{"SVTolRel",p.diis.SVTolRel},{"Trust",p.ladderGdm.Trust},
            {"ethresh",p.ladderEThresh},{"stall",p.ladderStall},{"floor",p.ladderFloor},{"switchat",p.ladderSwitchAt},
            {"signal",p.ladderSignal==SCFAccelerators::ScheduleSignal::EnergyChange ? "EnergyChange" : "Error"}};
}

RunSeed ResolveRunSeed(ChargeDensity::SeedStrategy stated, const Structure& st, const Pseudopotential::PseudopotentialSet* pp,
                              bool densityFock, const std::vector<int>& siteCharges, const std::vector<int>& siteSpins)
{
    RunSeed r;
    if (pp) r.key=SeedKeyOf(*pp, st);
    r.key.siteCharges=siteCharges;
    r.key.siteSpins=siteSpins;
    r.seed=ResolveSeed(stated, st, {.densityFock = densityFock, .periodic = st.MakeSymmetry(Symmetry::kAdmissionTolerance)->HasLattice(),
                                    .chargesStated = !siteCharges.empty(),
                                    .covered = [key=r.key](int Z){ return SeedLibraryCovers(key, Z); }});
    return r;
}

SCFAccelerators::Type AcceleratorTypeByName(const std::string& name, const std::string& who)
{
    using T=SCFAccelerators::Type;
    if (name=="DIIS")   return T::DIIS;
    if (name=="GDM")    return T::GDM;
    if (name=="Ladder") return T::Ladder;
    policy::RefuseAs("AcceleratorType", who+": unknown accelerator type \""+name+"\" (expected DIIS | GDM | Ladder)");
}

SCFAccelerators::Type MolecularAcceleratorType(SCFAccelerators::Type asked, const std::string& who)
{
    using T=SCFAccelerators::Type;
    if (asked==T::Auto) return T::DIIS;
    if (asked==T::Null)
        policy::RefuseAs("AcceleratorType", who+": the Null accelerator -- a molecular SCF has no density mixer to pair it with (expected "
                                            "DIIS | GDM | Ladder, or Auto = DIIS)");
    return asked;
}

double DIISEngageGate(double stated, int Z, bool densityFock)
{
    return stated>0.0 ? stated : (densityFock ? 100.0 : Z*Z*0.08/32);
}

void RequireVxcAdmitted(const Structure& st, Hamiltonian::VxcFit fit, const qcMesh::MeshParams& mesh, const std::string& where)
{
    // THE PIN-30 PRE-FLIGHT (Group face step 3b; controller ruling 2026-11-21: a policy, so it reports under `policy` either way).  It asks
    // the structure's group alone, nothing built; the Vxc build then checks the built objects' Forces, and a disagreement there THROWS.
    if (const auto a=Hamiltonian::AdmitVxc(st, fit, mesh); !a)
        policy::Refuse("AdmitVxc", "the stated Vxc fit and XC mesh on this "+where, a.Error());
    policy::Report("AdmitVxc", "the stated Vxc fit and XC mesh on this "+where+": admitted", policy::Verdict::Choose,
                   "pin 30 pre-flight, asked of the structure's group G alone: what the fit and the XC points force lies within G; the Vxc "
                   "build checks the same of the BUILT objects (their Forces) and throws if the two ever disagree");
}
} //namespace qchem
