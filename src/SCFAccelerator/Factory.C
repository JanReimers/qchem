// File: SCFAccelerator/Factory.H  Create SCFAccelerators
module;
#include <nlohmann/json_fwd.hpp>
export module qchem.SCFAccelerator.Factory;
export import qchem.SCFAccelerator;

export namespace SCFAcceleratorF
{
    enum class Type {DIIS,GDM};
    SCFAccelerator* Factory(Type,const nlohmann::json& js);
}

