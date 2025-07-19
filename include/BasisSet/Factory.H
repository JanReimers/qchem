// File:: BasisSet/Factory.H  Interfaces for various fctories.
#ifndef _BasisSet_Factory_H_
#define _BasisSet_Factory_H_

#include <nlohmann/json.hpp>

import qchem.BasisSet;
import qchem.Cluster;
import qchem.Symmetry.ElectronConfiguration;

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

#endif //_BasisSet_Factory_H_