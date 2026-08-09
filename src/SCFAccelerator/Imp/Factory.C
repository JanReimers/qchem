// File: SCFAccelerator/Factory.H  Create SCFAccelerators
module;
#include <cassert>
#include <memory>
#include <vector>
#include <nlohmann/json.hpp>
module qchem.SCFAccelerator.Factory;
import qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS;
import qchem.SCFAccelerator.Internal.SCFAcceleratorGDM;
import qchem.SCFAccelerator.Internal.SCFAcceleratorLadder;
import qchem.SCFAccelerator.Internal.SCFIrrepAcceleratorNull;
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
            double EMax=js["EMax"].template get<double>(),EMin=js["EMin"].template get<double>(),SVTol=js["SVTol"].template get<double>();
            acc=new SCFAcceleratorDIIS({Nproj,EMax,EMin,SVTol});
            break;
        }
        case Type::GDM:
        {
            double EMax  = js.contains("EMax")  ? js["EMax"].template get<double>()  : 1e10;
            double Trust = js.contains("Trust") ? js["Trust"].template get<double>() : 0.1;
            acc=new SCFAcceleratorGDM({EMax,Trust});
            break;
        }
        case Type::Ladder:
        {
            // DIIS does the heavy lifting; GDM finishes once DIIS runs out of steam.
            size_t Nproj = js.contains("NProj") ? js["NProj"].template get<size_t>() : 4;
            double EMax  = js.contains("EMax")  ? js["EMax"].template get<double>()  : 1.0;
            double EMin  = js.contains("EMin")  ? js["EMin"].template get<double>()  : 1e-7;
            double SVTol = js.contains("SVTol") ? js["SVTol"].template get<double>() : 5e-9;
            double Trust   = js.contains("Trust")   ? js["Trust"].template get<double>()   : 0.1;
            double ethresh = js.contains("ethresh") ? js["ethresh"].template get<double>() : 1e-8;
            int    stall   = js.contains("stall")   ? js["stall"].template get<int>()      : 5;
            double floor   = js.contains("floor")   ? js["floor"].template get<double>()   : 1e-8;
            double swat    = js.contains("switchat")? js["switchat"].template get<double>(): 0.0;
            std::vector<std::unique_ptr<SCFAccelerator>> rungs;
            rungs.emplace_back(new SCFAcceleratorDIIS({Nproj,EMax,EMin,SVTol}));
            rungs.emplace_back(new SCFAcceleratorGDM ({1e10,Trust})); //always steps once it is the active rung
            acc=new SCFAcceleratorLadder(std::move(rungs),ethresh,stall,floor,swat);
            break;
        }
        case Type::Null: acc=new SCFAcceleratorNull(); break;
    }
    
    assert(acc);
    return acc;
}

// The SOLID front door: the cSCFAccelerator twin.  Typed rather than json (see the interface doc).
cSCFAccelerator* Factory(Type type, const SolidAcceleratorOptions& o)
{
    switch (type)
    {
        case Type::Null: return new tSCFAcceleratorNull<dcmplx>();
        case Type::DIIS: return new cSCFAcceleratorDIIS({o.nProj,o.fdMax,o.fdMin,o.svTol});
        case Type::GDM : return new cSCFAcceleratorGDM ({o.gdmFDMax,o.gdmTrust});
        case Type::Ladder:
        {
            // DIIS does the heavy lifting; GDM finishes once DIIS runs out of steam.  The hand-off watches
            // the ENERGY CHANGE, not the residual: a non-variational collocation SCF limit-cycles [F,D]
            // above the fit floor while E settles, so the residual never signals the switch.
            std::vector<std::unique_ptr<cSCFAccelerator>> rungs;
            rungs.emplace_back(new cSCFAcceleratorDIIS({o.nProj,o.fdMax,o.fdMin,o.svTol}));
            rungs.emplace_back(new cSCFAcceleratorGDM ({o.gdmFDMax,o.ladderGdmTrust}));
            return new cSCFAcceleratorLadder(std::move(rungs), o.ladderEThresh, o.ladderStall,
                                             o.ladderFloor, o.ladderSwitchAt, ScheduleSignal::EnergyChange);
        }
    }
    assert(false && "SCFAccelerators::Factory: unhandled Type");
    return nullptr;
}

} //namespace 
