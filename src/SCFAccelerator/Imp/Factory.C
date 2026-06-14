// File: SCFAccelerator/Factory.H  Create SCFAccelerators
module;
#include <cassert>
#include <vector>
#include <nlohmann/json.hpp>
module qchem.SCFAccelerator.Factory;
import qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS;
import qchem.SCFAccelerator.Internal.SCFAcceleratorGDM;
import qchem.SCFAccelerator.Internal.SCFAcceleratorLadder;
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
            double Trust = js.contains("Trust") ? js["Trust"].template get<double>() : 0.1;
            double floor = js.contains("floor") ? js["floor"].template get<double>() : 1e-4;
            int    stall = js.contains("stall") ? js["stall"].template get<int>()    : 5;
            std::vector<SCFAccelerator*> rungs;
            rungs.push_back(new SCFAcceleratorDIIS({Nproj,EMax,EMin,SVTol}));
            rungs.push_back(new SCFAcceleratorGDM ({1e10,Trust})); //always steps once it is the active rung
            acc=new SCFAcceleratorLadder(rungs,floor,stall);
            break;
        }
        
    }
    
    assert(acc);
    return acc;
}

} //namespace 
