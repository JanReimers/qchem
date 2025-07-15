// File: Mesh/Factory.H  Create various mesh types.
#ifndef _Mesh_Factory_H_
#define _Mesh_Factory_H_

#include <Mesh/fwd.H>
#include <nlohmann/json.hpp>
#include <Mesh/MeshParams.H>
import Mesh;
import RadialMesh;

namespace MeshF
{
    RadialMesh* Factory(qchem::RadialType,const nlohmann::json& js);
          Mesh* Factory(qchem::AngleType,const nlohmann::json& js); 
}

#endif //_Mesh_Factory_H_