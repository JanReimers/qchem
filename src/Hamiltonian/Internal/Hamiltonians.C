// File:: Hamiltonian/Internal/Hamiltonians.C  Create fully implemented Hamiltonians
export module qchem.Hamiltonian.Internal.Hamiltonians;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Hamiltonian.Internal.Hamiltonian;
import qchem.Hamiltonian.Types;
import qchem.Mesh1;

export namespace qchem::Hamiltonian
{
//
//  1 Electron
//
class Ham_1E : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_1E(const st_t& st);
};

//
//  Un-polarized Hartree-Fock
//
class Ham_HF_U : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_HF_U(const st_t& st);
};

//
// Un-polarized DFT, fitted Vee, and fitted DFT like Vxc
//
class Ham_DFT_U : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_DFT_U(const st_t& st,double alpha_ex, const qcMesh1::MeshParams&, const bs_t* bs);
    Ham_DFT_U(const st_t& st,ExFunctional*  , const qcMesh1::MeshParams&, const bs_t* bs);
};

// Un-polarized LDA: Dirac exchange + VWN5 correlation as SEPARATE terms, so the correlation energy is the
// correct E_c = integral eps_c rho (not the exchange virial 3/4<rho|Vc>).  Exchange and correlation share
// one Vxc fit basis.  This is the "real" LSDA Hamiltonian for the NIST atomic oracle.
class Ham_DFTcorr_U : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_DFTcorr_U(const st_t& st, const qcMesh1::MeshParams&, const bs_t* bs);
};

//
// Plane-wave LDA Kohn-Sham (dcmplx): kinetic + external(pseudo) + Hartree + Dirac exchange + VWN5
// correlation, assembled from the qcHamiltonian plane-wave terms (PWTerms).  Unlike the molecular DFT
// Hamiltonians it needs NO fit basis / mesh -- the plane-wave basis assembles every matrix in G-space
// and carries the pseudopotential (configured by the BasisSet factory), so this only needs the
// structure for the external term's structure factor.  Pair with a plane-wave BasisSet + cSCFIterator.
//
class Ham_PW_DFT : public virtual cHamiltonian, private cHamiltonianImp
{
public:
    Ham_PW_DFT(const st_t& st);
};

//
//  Polarized Hartree-Fock.
//
class Ham_HF_P : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_HF_P(const st_t& st);
};


//
// Un-polarized DFT, fitted Vee, and fitted DFT like Vxc
//
class Ham_DFT_P : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_DFT_P(const st_t& st,double alpha_ex, const qcMesh1::MeshParams&, const bs_t* or_bs);
    Ham_DFT_P(const st_t& st,ExFunctional*  , const qcMesh1::MeshParams&, const bs_t* or_bs);
};


class Ham_DHF_1E : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_DHF_1E(const st_t& st);
};

//
//  Dirac-Hartree-Fock.
//
class Ham_DHF_U : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_DHF_U(const st_t& st);
};

class Ham_DHF_P : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_DHF_P(const st_t& st);
};

} //namespace