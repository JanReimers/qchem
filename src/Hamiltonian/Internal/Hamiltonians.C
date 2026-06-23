// File:: Hamiltonians.H  Create fully implemented Hamiltonians
export module qchem.Hamiltonian.Internal.Hamiltonians;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Hamiltonian.Internal.Hamiltonian;
import qchem.Hamiltonian.Types;
import qchem.Mesh;

export namespace qchem::Hamiltonian
{
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
    Ham_DFT_U(const cl_t& cl,double alpha_ex, const MeshParams&, const bs_t* bs);
    Ham_DFT_U(const cl_t& cl,ExFunctional*  , const MeshParams&, const bs_t* bs);
};

// Un-polarized LDA: Dirac exchange + VWN5 correlation as SEPARATE terms, so the correlation energy is the
// correct E_c = integral eps_c rho (not the exchange virial 3/4<rho|Vc>).  Exchange and correlation share
// one Vxc fit basis.  This is the "real" LSDA Hamiltonian for the NIST atomic oracle.
class Ham_DFTcorr_U : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_DFTcorr_U(const cl_t& cl, const MeshParams&, const bs_t* bs);
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
    Ham_DFT_P(const cl_t& cl,double alpha_ex, const MeshParams&, const bs_t* or_bs);
    Ham_DFT_P(const cl_t& cl,ExFunctional*  , const MeshParams&, const bs_t* or_bs);
};


class Ham_DHF_1E : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_DHF_1E(const cl_t& cl);
};

//
//  Dirac-Hartree-Fock.
//
class Ham_DHF_U : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_DHF_U(const cl_t& cl);
};

class Ham_DHF_P : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_DHF_P(const cl_t& cl);
};

} //namespace