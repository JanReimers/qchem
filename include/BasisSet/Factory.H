// File:: Factory.H  Interfaces for various fctories.
#ifndef _Factory_H_
#define _Factory_H_

#include <types.H>
#include <nlohmann/json.hpp>

namespace BasisSetAtom
{
    enum class Type {Slater,Gaussian,BSpline,Slater_RKB,Gaussian_RKB,BSpline_RKB};
    BasisSet* Factory(Type,const nlohmann::json& js,const ElectronConfiguration& ec);
    BasisSet* Factory(Type,const nlohmann::json& js,size_t Z); //Use atomic number Z to deduce the ElectronConfiguration
    BasisSet* Factory(const nlohmann::json& js,size_t Z); //Type is in the json object
}

namespace BasisSetMolecule
{
    BasisSet* Factory(const nlohmann::json& js,const Cluster* cl);
}

#endif //_Factory_H_