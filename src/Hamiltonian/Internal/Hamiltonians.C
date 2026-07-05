// File:: Hamiltonian/Internal/Hamiltonians.C  Create fully implemented Hamiltonians
module;
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <initializer_list>
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
class Ham_1E : public virtual rHamiltonian, private rHamiltonianImp
{
public:
    Ham_1E(const st_t& st);
};

//
//  Un-polarized Hartree-Fock
//
class Ham_HF_U : public virtual rHamiltonian, private rHamiltonianImp
{
public:
    Ham_HF_U(const st_t& st);
};

//
// Un-polarized DFT, fitted Vee, and fitted DFT like Vxc
//
class Ham_DFT_U : public virtual rHamiltonian, private rHamiltonianImp
{
public:
    Ham_DFT_U(const st_t& st,double alpha_ex, const qcMesh::MeshParams&, const rbs_t* bs);
    Ham_DFT_U(const st_t& st,ExFunctional*  , const qcMesh::MeshParams&, const rbs_t* bs);
};

// Un-polarized LDA: Dirac exchange + VWN5 correlation as SEPARATE terms, so the correlation energy is the
// correct E_c = integral eps_c rho (not the exchange virial 3/4<rho|Vc>).  Exchange and correlation share
// one Vxc fit basis.  This is the "real" LSDA Hamiltonian for the NIST atomic oracle.
class Ham_DFTcorr_U : public virtual rHamiltonian, private rHamiltonianImp
{
public:
    Ham_DFTcorr_U(const st_t& st, const qcMesh::MeshParams&, const rbs_t* bs);  //!< Dirac exchange + VWN5
    //! Generic separate-terms LSDA: FittedVxc(\a exchange) [3/4 virial energy, exact for exchange] +
    //! FittedVcorr(\a correlation) [E_c = integral eps_c rho, the correct correlation energy].  The
    //! correlation functional MUST provide GetEpsXc (the energy density).  Takes ownership of both.  The
    //! libxc LSDA path uses this with Libxc_LDA exchange/correlation; the default ctor with Slater+VWN.
    Ham_DFTcorr_U(const st_t& st, ExFunctional* exchange, ExFunctional* correlation,
                  const qcMesh::MeshParams&, const rbs_t* bs);
};

// Spin-native (polarized) LSDA: Dirac exchange (FittedVxcPol) + spin-native VWN5 correlation
// (FittedVcorrPol) as separate terms sharing one Vxc fit basis.  The U = ζ=0 collapse; this is the
// open-shell / magnetism path (OpenWork B).  Correlation couples both spin channels -- see FittedVcorrPol.
class Ham_DFTcorr_P : public virtual rHamiltonian, private rHamiltonianImp
{
public:
    Ham_DFTcorr_P(const st_t& st, const qcMesh::MeshParams&, const rbs_t* bs);
};

//! Spin-native LSDA PSEUDO-atom/molecule: kinetic + V_loc(r) (the pseudized replacement for the bare
//! nuclear attraction, mesh-quadratured) + the KB-separable non-local projectors + Hartree + Dirac exchange
//! + VWN5 correlation + ion-ion (Zion cores; a direct pair sum for a molecule, ZERO for a lone atom).  NO
//! Ven.  The valence electron count comes from the structure (build the Atom(s) with charge = Z - valence)
//! and the angular channels from the electron configuration.  The non-local model may be null (local-only,
//! the over-bound stepping stone).  The SAME term list serves a single pseudo-atom and a pseudo-molecule --
//! the PP terms loop over all atoms and the ion-ion pair sum vanishes when there is only one.
//!
//! \a polarized selects the spin-native (open-shell) exchange-correlation -- FittedVxcPol + FittedVcorrPol,
//! the open-shell/magnetism path -- vs the \f$\zeta=0\f$ unpolarized collapse (FittedVxc + FittedVcorr, the
//! closed-shell efficiency case).  Everything else (kinetic/PP/Hartree/ion-ion) is spin-agnostic.
class Ham_PP : public virtual rHamiltonian, private rHamiltonianImp
{
public:
    //! Explicit models (shared with the caller).  \a vloc is the COMBINED local model: its real-space view
    //! (Vloc) feeds PP_Local and its \f$Z_{ion}\f$ feeds the ion-ion term.  \a sep may be null (local-only).
    Ham_PP(const st_t& st, std::shared_ptr<const Pseudopotential::LocalPotential> vloc,
           std::shared_ptr<const Pseudopotential::SeparablePotential_R> sep,
           const qcMesh::MeshParams&, const rbs_t* bs, bool polarized);
    //! Convenience: look up + OWN the GTH local + KB non-local models for \a element at valence \a q (LDA).
    Ham_PP(const st_t& st, const std::string& element, int q, const qcMesh::MeshParams&, const rbs_t* bs,
           bool polarized);
    //! Multi-species convenience: name each \a (element, valence); the GTH database is looked up per species
    //! and a per-Z router (MultiSpecies_Local/Separable) is built + OWNED, so each atom gets its OWN
    //! pseudopotential (the assembly + ion-ion already route on the atoms' itsZ).  The single-species case is
    //! a 1-element list.  E.g. Ham_PP(st, {{"Si",4},{"O",6}}, mesh, bs, polarized).
    Ham_PP(const st_t& st, const std::vector<std::pair<std::string,int>>& species,
           const qcMesh::MeshParams&, const rbs_t* bs, bool polarized);
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
//   - convenience ctors: name the element(s)/functional/valence and the database (GetGTH) is looked up
//     and OWNED here -- the one-call plane-wave LDA Hamiltonian, single- or multi-species.
//
class Ham_PW_DFT : public virtual cHamiltonian, private cHamiltonianImp
{
public:
    //! \a bs (the composite plane-wave basis) is the density-fit-basis source: BuildTerms builds the Hartree
    //! fit basis from it ONCE, mirroring the molecular Ham DFT ctors that take \c bs for FittedVee.
    //! Explicit-models ctor: the caller owns the local + optional KB nonlocal models (non-owning here).
    Ham_PW_DFT(const st_t& st, const cbs_t* bs, const Pseudopotential::LocalPotential* loc,
               const Pseudopotential::SeparablePotential* nl=nullptr);
    //! Single-species convenience ctor: look up + own the GTH PP for \a element.
    Ham_PW_DFT(const st_t& st, const cbs_t* bs, const std::string& element,
               const std::string& functional="LDA", int valence=0);
    //! Multi-species convenience ctor: name each (element, valence); the database is looked up per species
    //! and a per-Z router model (MultiSpecies_*) is built + OWNED -- e.g. Ham_PW_DFT(st, bs, {{"Na",1},{"F",7}}).
    Ham_PW_DFT(const st_t& st, const cbs_t* bs, std::initializer_list<std::pair<std::string,int>> species,
               const std::string& functional="LDA");
private:
    void BuildTerms(const st_t& st, const cbs_t* bs, const Pseudopotential::LocalPotential* loc,
                    const Pseudopotential::SeparablePotential* nl);
    //! Look up each (element, valence) from the GTH database, build + OWN the (per-Z router) local +
    //! separable models, and assemble the terms against them.  The single-species ctor is the 1-species case.
    void BuildFromGTH(const st_t& st, const cbs_t* bs, const std::vector<std::pair<std::string,int>>& species,
                      const std::string& functional);
    std::shared_ptr<const Pseudopotential::LocalPotential>     itsOwnedLocal;  //!< owned model (convenience ctors); null for explicit
    std::shared_ptr<const Pseudopotential::SeparablePotential> itsOwnedSep;
};

//
//  Polarized Hartree-Fock.
//
class Ham_HF_P : public virtual rHamiltonian, private rHamiltonianImp
{
public:
    Ham_HF_P(const st_t& st);
};


//
// Un-polarized DFT, fitted Vee, and fitted DFT like Vxc
//
class Ham_DFT_P : public virtual rHamiltonian, private rHamiltonianImp
{
public:
    Ham_DFT_P(const st_t& st,double alpha_ex, const qcMesh::MeshParams&, const rbs_t* or_bs);
    Ham_DFT_P(const st_t& st,ExFunctional*  , const qcMesh::MeshParams&, const rbs_t* or_bs);
};


class Ham_DHF_1E : public virtual rHamiltonian, private rHamiltonianImp
{
public:
    Ham_DHF_1E(const st_t& st);
};

//
//  Dirac-Hartree-Fock.
//
class Ham_DHF_U : public virtual rHamiltonian, private rHamiltonianImp
{
public:
    Ham_DHF_U(const st_t& st);
};

class Ham_DHF_P : public virtual rHamiltonian, private rHamiltonianImp
{
public:
    Ham_DHF_P(const st_t& st);
};

} //namespace