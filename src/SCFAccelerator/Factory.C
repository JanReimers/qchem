// File: SCFAccelerator/Factory.H  Create SCFAccelerators

#include <SCFAccelerator/Factory.H>
#include "SCFAccelerator_DIIS.H"
#include "SCFAccelerator_Null.H"
#include <cassert>
#include <nlohmann/json.hpp>
using json = nlohmann::json;

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
            acc=new SCFAccelerator_DIIS({Nproj,EMax,EMin,SVTol});
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
