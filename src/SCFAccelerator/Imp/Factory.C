// File: SCFAccelerator/Factory.H  Create SCFAccelerators
module;
#include <cassert>
#include <nlohmann/json.hpp>
module qchem.SCFAccelerator.Factory;
import qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS;
import qchem.SCFAccelerator.Internal.SCFAcceleratorGDM;
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
            double EMax = js.contains("EMax") ? js["EMax"].template get<double>() : 1e10;
            acc=new SCFAcceleratorGDM({EMax});
            break;
        }
        
    }
    
    assert(acc);
    return acc;
}

} //namespace 
