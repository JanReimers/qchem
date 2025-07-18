// File: Hamiltonian/Factory.H  Create some charge densitytypes.
#ifndef _Hamiltonian_Factory_H_
#define _Hamiltonian_Factory_H_

#include <Mesh/fwd.H>
#include <BasisSet/fwd.H>
#include <Hamiltonian/fwd.H>
#include <memory>

import Mesh;

typedef std::shared_ptr<const Cluster> cl_t;
namespace HamiltonianF
{
    enum class Model {E1,HF,DE1,DHF}; //E1 is 1 electron. DE1 is Dirac 1 electron.
    enum class Pol   {UnPolarized,Polarized};
    Hamiltonian* Factory(Model,Pol,const cl_t& cl);
    Hamiltonian* Factory(Pol,const cl_t& cl,ExFunctional* , const MeshParams&, const BasisSet*); //DFT version
    Hamiltonian* Factory(Pol,const cl_t& cl,double alpha  , const MeshParams&, const BasisSet*); //DFT version

}

#endif //_Hamiltonian_Factory_H_
