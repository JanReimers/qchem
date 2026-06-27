// File:: Hamiltonian/Internal/Hamiltonians.C  Create fully implemented Hamiltonians
module;
#include <memory>
#include <string>
export module qchem.Hamiltonian.Internal.Hamiltonians;
import qchem.Hamiltonian.Internal.ExFunctional;
import qchem.Hamiltonian.Internal.Hamiltonian;
import qchem.Hamiltonian.Types;
import qchem.Mesh;
import qchem.Pseudopotential.LocalPotential;      // the PW pseudopotential model the term owns (Ham_PW_DFT ctor)
import qchem.Pseudopotential.SeparablePotential;
import qchem.Pseudopotential.GTH_Potentials;      // GTH_PP / GetGTH (the convenience ctor looks up + OWNS the model)

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
    Ham_DFT_U(const st_t& st,double alpha_ex, const qcMesh::MeshParams&, const bs_t* bs);
    Ham_DFT_U(const st_t& st,ExFunctional*  , const qcMesh::MeshParams&, const bs_t* bs);
};

// Un-polarized LDA: Dirac exchange + VWN5 correlation as SEPARATE terms, so the correlation energy is the
// correct E_c = integral eps_c rho (not the exchange virial 3/4<rho|Vc>).  Exchange and correlation share
// one Vxc fit basis.  This is the "real" LSDA Hamiltonian for the NIST atomic oracle.
class Ham_DFTcorr_U : public virtual Hamiltonian, private HamiltonianImp
{
public:
    Ham_DFTcorr_U(const st_t& st, const qcMesh::MeshParams&, const bs_t* bs);
};

//
// Plane-wave LDA Kohn-Sham (dcmplx): kinetic + external(pseudo) + Hartree + Dirac exchange + VWN5
// correlation, assembled from the qcHamiltonian plane-wave terms (PWTerms).  Unlike the molecular DFT
// Hamiltonians it needs NO fit basis / mesh -- the plane-wave basis assembles every matrix in G-space.
// The pseudopotential MODEL (local form factor + optional KB nonlocal) is handed to the external term,
// which asks the basis to assemble the matrix from it (the pseudo-wall).  Pair with a plane-wave
// BasisSet + cSCFIterator.  Two ways to supply the model:
//   - explicit ctor: the caller owns the local + optional KB nonlocal models (non-owning here; they must
//     outlive the SCF run);
//   - convenience ctor: name the element/functional/valence and the database (GetGTH) is looked up and
//     OWNED here -- the one-call plane-wave LDA Hamiltonian.
//
class Ham_PW_DFT : public virtual cHamiltonian, private cHamiltonianImp
{
public:
    Ham_PW_DFT(const st_t& st, const Pseudopotential::LocalPotential* loc,
               const Pseudopotential::SeparablePotential* nl=nullptr);
    Ham_PW_DFT(const st_t& st, const std::string& element,
               const std::string& functional="LDA", int valence=0);
private:
    void BuildTerms(const st_t& st, const Pseudopotential::LocalPotential* loc,
                    const Pseudopotential::SeparablePotential* nl);
    std::unique_ptr<Pseudopotential::GTH_PP> itsOwnedPP;  //!< owned model (convenience ctor); null for the explicit ctor
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
    Ham_DFT_P(const st_t& st,double alpha_ex, const qcMesh::MeshParams&, const bs_t* or_bs);
    Ham_DFT_P(const st_t& st,ExFunctional*  , const qcMesh::MeshParams&, const bs_t* or_bs);
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