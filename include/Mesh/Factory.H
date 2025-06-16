// File: Mesh/Factory.H  Create various mesh types.
#ifndef _Mesh_Factory_H_
#define _Mesh_Factory_H_

#include <types.H>
#include <nlohmann/json.hpp>

namespace MeshF
{
    enum class  RadialType {Log,MHL};
    enum class AngularType {EulerMaclaren,Gauss,GaussLegendre};
    RadialMesh* Factory( RadialType,const nlohmann::json& js);
          Mesh* Factory(AngularType,const nlohmann::json& js); 
}

#endif //_Mesh_Factory_H_