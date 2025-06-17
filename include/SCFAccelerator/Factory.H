// File: SCFAccelerator/Factory.H  Create SCFAccelerators
#ifndef _SCFAccelerator_Factory_H_
#define _SCFAccelerator_Factory_H_

#include <types.H>
#include <nlohmann/json.hpp>

namespace SCFAcceleratorF
{
    enum class Type {DIIS,GDM};
    SCFAccelerator* Factory(Type,const nlohmann::json& js);
}

#endif //_SCFAccelerator_Factory_H_