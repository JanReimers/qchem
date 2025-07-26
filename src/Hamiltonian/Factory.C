// File: Hamiltonian/Factory.C Construct and return various Hamiltonian types.
module;
#include <memory>
export module qchem.Hamiltonian.Factory;
export import qchem.Hamiltonian;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Mesh;
import qchem.Cluster;
import qchem.BasisSet;

typedef std::shared_ptr<const Cluster> cl_t;
export namespace HamiltonianF
{
    enum class Model {E1,HF,DE1,DHF}; //E1 is 1 electron. DE1 is Dirac 1 electron.
    enum class Pol   {UnPolarized,Polarized};
    Hamiltonian* Factory(Model,Pol,const cl_t& cl);
    Hamiltonian* Factory(Pol,const cl_t& cl,ExFunctional* , const MeshParams&, const BasisSet*); //DFT version
    Hamiltonian* Factory(Pol,const cl_t& cl,double alpha  , const MeshParams&, const BasisSet*); //DFT version

}

