// File:: Hamiltonians.H  Create fully implemented Hamiltonians

#include "Imp/Hamiltonian/Hamiltonian.H"

//
//  Un-polarized Hartree-Fock
//
class Ham_HF_U : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_HF_U(cl_t& cl);
};

//
// Un-polarized Semi Hartree-Fock, full Vee, and fitted DFT like Vxc
//
class Ham_SHF_U : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_SHF_U(cl_t& cl,double alpha_ex, const MeshParams&, const BasisSet* or_bs);
    Ham_SHF_U(cl_t& cl,ExFunctional*  , const MeshParams&, const BasisSet* or_bs);
};

//
// Un-polarized DFT, fitted Vee, and fitted DFT like Vxc
//
class Ham_DFT_U : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_DFT_U(cl_t& cl,double alpha_ex, const MeshParams&, const BasisSet* fit_bs);
    Ham_DFT_U(cl_t& cl,ExFunctional*  , const MeshParams&, const BasisSet* fit_bs);
};

//
//  Polarized Hartree-Fock.
//
class Ham_HF_P : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_HF_P(cl_t& cl);
};

//
// Un-polarized Semi Hartree-Fock, full Vee, and fitted DFT like Vxc
//
class Ham_SHF_P : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_SHF_P(cl_t& cl,double alpha_ex, const MeshParams&, const BasisSet* or_bs);
    Ham_SHF_P(cl_t& cl,ExFunctional*  , const MeshParams&, const BasisSet* or_bs);
};

//
// Un-polarized DFT, fitted Vee, and fitted DFT like Vxc
//
class Ham_DFT_P : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_DFT_P(cl_t& cl,double alpha_ex, const MeshParams&, const BasisSet* or_bs);
    Ham_DFT_P(cl_t& cl,ExFunctional*  , const MeshParams&, const BasisSet* or_bs);
};

//
//  Dirac-Hartree-Fock (which is polarized by definition).
//
class Ham_DHF : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_DHF(cl_t& cl);
};
