// File: Hamiltonian/Factory.C Construct and return various Hamiltonian types.
module;
#include <memory>
export module qchem.Hamiltonian.Factory;
export import qchem.Hamiltonian;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Hamiltonian.Types;
import qchem.Mesh;
import qchem.Structure;


export namespace qchem::Hamiltonian
{
    typedef std::shared_ptr<const Structure> st_t;
    enum class Model {E1,HF,DE1,DHF}; //E1 is 1 electron. DE1 is Dirac 1 electron.
    enum class Pol   {UnPolarized,Polarized};
    Hamiltonian* Factory(Model,Pol,const st_t& st);
    Hamiltonian* Factory(Pol,const st_t& st,ExFunctional* , const qcMesh::MeshParams&, const bs_t*); //DFT version
    Hamiltonian* Factory(Pol,const st_t& st,double alpha  , const qcMesh::MeshParams&, const bs_t*); //DFT version

} // namespace
