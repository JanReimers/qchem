// File: Mesh/Factory.C  Create various mesh types.
module;
#include <nlohmann/json_fwd.hpp>
export module qchem.Mesh.Factory;
export import Mesh;
export import RadialMesh;

export namespace MeshF
{
    RadialMesh* Factory(qchem::RadialType,const nlohmann::json& js);
          Mesh* Factory(qchem::AngleType,const nlohmann::json& js); 
}
