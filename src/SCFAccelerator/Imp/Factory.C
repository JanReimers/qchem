// File: SCFAccelerator/Factory.H  Create SCFAccelerators
module;
#include <cassert>
#include <nlohmann/json.hpp>
using json = nlohmann::json;
module qchem.SCFAccelerator.Factory;
import qchem.SCFAccelerator.Internal.SCFAcceleratorDIIS;

namespace SCFAcceleratorF
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
            assert(false);
            break;
        }
        
    }
    
    assert(acc);
    return acc;
}

} //namespace 
