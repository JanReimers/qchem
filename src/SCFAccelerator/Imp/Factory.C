// File: SCFAccelerator/Factory.H  Create SCFAccelerators
module;
#include <cassert>
#include <memory>
#include <vector>
#include <iostream>
#include <stdexcept>
#include <string>
#include <nlohmann/json.hpp>
module qchem.SCFAccelerator.Factory;
import qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS;
import qchem.SCFAccelerator.Internal.SCFAcceleratorGDM;
import qchem.SCFAccelerator.Internal.SCFAcceleratorLadder;
import qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull;
import qchem.Reporting;   // each Ladder rung files itself into a `rungs` row opened here, its builder
using json = nlohmann::json;

namespace qchem::SCFAccelerators
{

SCFAccelerator* Factory(Type type,const nlohmann::json& js)
{
    SCFAccelerator* acc=0;
    switch (type)
    {
        case Type::DIIS:
        {
            size_t Nproj=js["NProj"].template get<size_t>();
            // "EMin" (the old DIISParams::FDMin) is still accepted in the json but was never read (D-DIIS-FDMIN).  DESIGNATED, so a
            // field deleted from DIISParams cannot shift the rest: the positional {Nproj,EMax,EMin,SVTol} silently became
            // SVTol=EMin, SVTolRel=SVTol when FDMin went, and moved every molecular DIIS run (Si2_PP_U, A_SL_HF_ion, 2026-10-30).
            double EMax=js["EMax"].template get<double>(),SVTol=js["SVTol"].template get<double>();
            acc=new SCFAcceleratorDIIS({.Nproj=Nproj, .EngageBelowFD=EMax, .SVTol=SVTol});
            break;
        }
        case Type::GDM:
        {
            double EMax  = js.contains("EMax")  ? js["EMax"].template get<double>()  : 1e10;
            double Trust = js.contains("Trust") ? js["Trust"].template get<double>() : 0.1;
            acc=new SCFAcceleratorGDM({.EngageBelowFD=EMax, .Trust=Trust});
            break;
        }
        case Type::Ladder:
        {
            // DIIS does the heavy lifting; GDM finishes once DIIS runs out of steam.
            size_t Nproj = js.contains("NProj") ? js["NProj"].template get<size_t>() : 4;
            double EMax  = js.contains("EMax")  ? js["EMax"].template get<double>()  : 1.0;
            double SVTol = js.contains("SVTol") ? js["SVTol"].template get<double>() : 5e-9;
            double Trust   = js.contains("Trust")   ? js["Trust"].template get<double>()   : 0.1;
            double ethresh = js.contains("ethresh") ? js["ethresh"].template get<double>() : 1e-8;
            int    stall   = js.contains("stall")   ? js["stall"].template get<int>()      : 5;
            double floor   = js.contains("floor")   ? js["floor"].template get<double>()   : 1e-8;
            double swat    = js.contains("switchat")? js["switchat"].template get<double>(): 0.0;
            // The tail trigger's signal: "Error" ([F,D], the door's default) or "EnergyChange" (|ΔE/E|) -- anything else is refused.
            ScheduleSignal signal=ScheduleSignal::Error;
            if (js.contains("signal"))
            {
                const std::string sg=js["signal"].template get<std::string>();
                if      (sg=="EnergyChange") signal=ScheduleSignal::EnergyChange;
                else if (sg!="Error") throw std::invalid_argument("SCFAccelerators::Factory: Ladder signal '"+sg+"' is not Error or EnergyChange");
            }
            std::vector<std::unique_ptr<SCFAccelerator>> rungs;
            { report::Row r("rungs"); rungs.emplace_back(new SCFAcceleratorDIIS({.Nproj=Nproj, .EngageBelowFD=EMax, .SVTol=SVTol})); }
            { report::Row r("rungs"); rungs.emplace_back(new SCFAcceleratorGDM ({.EngageBelowFD=1e10, .Trust=Trust})); } //always steps once it is the active rung
            // "startOnGDM": the seed is a RESTORED state (a molecule's or an atom's exact resume): start on the GDM rung, as the solid
            // door's ladderStartOnGDM does, so DIIS does not step the restored density away (NaF, 2026-10-02).
            const size_t start = js.contains("startOnGDM") && js["startOnGDM"].template get<bool>() ? rungs.size()-1 : 0;
            if (start>0) std::cout << "  *** SCF accelerator ladder: starting on rung " << start
                                   << " (GDM) -- the seed is a restored state ***" << std::endl;
            acc=new SCFAcceleratorLadder(std::move(rungs),ethresh,stall,floor,swat,signal,start);
            break;
        }
        case Type::Null: acc=new SCFAcceleratorNull(); break;
        case Type::Auto: throw std::invalid_argument("SCFAccelerators::Factory: Type::Auto must be resolved by the caller's policy first");
    }
    
    assert(acc);
    return acc;
}

// The SOLID front door: the typed twin (see the interface doc).  Same scalar-agnostic managers as the
// json door -- only the knob transport and the solid-tuned defaults differ.
SCFAccelerator* Factory(Type type, const SolidAcceleratorOptions& o)
{
    switch (type)
    {
        case Type::Auto: throw std::invalid_argument("SCFAccelerators::Factory: Type::Auto must be resolved by the caller's policy first");
        case Type::Null: return new SCFAcceleratorNull();
        case Type::DIIS: return new SCFAcceleratorDIIS(o.diis);
        case Type::GDM : return new SCFAcceleratorGDM (o.gdm);
        case Type::Ladder:
        {
            // DIIS does the heavy lifting; GDM finishes once DIIS runs out of steam.  The hand-off watches
            // the ENERGY CHANGE, not the residual: a non-variational collocation SCF limit-cycles [F,D]
            // above the fit floor while E settles, so the residual never signals the switch.
            std::vector<std::unique_ptr<SCFAccelerator>> rungs;
            { report::Row r("rungs"); rungs.emplace_back(new SCFAcceleratorDIIS(o.diis)); }
            { report::Row r("rungs"); rungs.emplace_back(new SCFAcceleratorGDM (o.ladderGdm)); }
            const size_t start=o.ladderStartOnGDM ? rungs.size()-1 : 0;
            if (start>0) std::cout << "  *** SCF accelerator ladder: starting on rung " << start
                                   << " (GDM) -- the seed is a restored state ***" << std::endl;
            return new SCFAcceleratorLadder(std::move(rungs), o.ladderEThresh, o.ladderStall,
                                            o.ladderFloor, o.ladderSwitchAt, o.ladderSignal, start);
        }
    }
    assert(false && "SCFAccelerators::Factory: unhandled Type");
    return nullptr;
}

} //namespace 
