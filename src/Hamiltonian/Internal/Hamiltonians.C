// File:: Hamiltonians.H  Create fully implemented Hamiltonians
export module qchem.Hamiltonian.Internal.Hamiltonians;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Hamiltonian.Internal.Hamiltonian;
import qchem.BasisSet;
import Mesh;

export {
//
//  1 Electron
//
class Ham_1E : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_1E(const cl_t& cl);
};

//
//  Un-polarized Hartree-Fock
//
class Ham_HF_U : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_HF_U(const cl_t& cl);
};

//
// Un-polarized DFT, fitted Vee, and fitted DFT like Vxc
//
class Ham_DFT_U : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_DFT_U(const cl_t& cl,double alpha_ex, const MeshParams&, const BasisSet* fit_bs);
    Ham_DFT_U(const cl_t& cl,ExFunctional*  , const MeshParams&, const BasisSet* fit_bs);
};

//
//  Polarized Hartree-Fock.
//
class Ham_HF_P : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_HF_P(const cl_t& cl);
};


//
// Un-polarized DFT, fitted Vee, and fitted DFT like Vxc
//
class Ham_DFT_P : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_DFT_P(const cl_t& cl,double alpha_ex, const MeshParams&, const BasisSet* or_bs);
    Ham_DFT_P(const cl_t& cl,ExFunctional*  , const MeshParams&, const BasisSet* or_bs);
};


class Ham_DHF_1E : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_DHF_1E(const cl_t& cl);
};

//
//  Dirac-Hartree-Fock (which is polarized by definition).
//
class Ham_DHF : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_DHF(const cl_t& cl);
};

} //export block