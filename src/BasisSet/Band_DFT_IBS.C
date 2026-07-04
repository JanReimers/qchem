// File: BasisSet/Band_DFT_IBS.C  Abstract DFT potential/energy assembly capability (mesh-integrated).
//
// PROVISIONAL interface (may be unified with the fit-based Orbital_DFT_IBS later).  It INVERTS the
// dependency between the Hamiltonian DFT terms and the concrete basis: the DFT terms live in
// qcHamiltonian and program against THIS abstract interface (in qcBasisSet); the concrete basis
// (PlaneWave_IBS in qcLattice_BS today; atom/molecule could implement it later) implements it.  Neither
// library depends on the other -- both depend only on qcBasisSet, exactly as the molecular terms depend
// on the abstract Orbital_DFT_IBS<T> (a term holds the abstract obs_t and dynamic_casts UP to the
// richer abstract capability; the concrete basis implements it).
//
// DESIGN PRINCIPLE -- "tell, don't ask", and do NOT conflate evaluation with integration:  a DFT term
// must not interrogate the basis for its internals (G-vectors, reciprocal lattice, mesh) and then do the
// integration itself.  It asks the basis high-level questions -- "the matrix of THIS potential", "the
// Hartree matrix for THIS density", "the integral of THIS function" -- and the basis owns HOW (radial
// quadrature for an atom, Becke grid for a molecule, FFT for plane waves; an atom has no G-vectors at
// all).  Hence the only inputs are real-space ScalarFunctions (the density IS one; the XC potential is a
// composed one) and there are deliberately NO getters here.
module;
export module qchem.BasisSet.Band_DFT_IBS;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Mesh_Integrated_IBS; // the field-operator base: Overlap(f) / Integral(f)
export import qchem.ScalarFunction;   // ScalarFunction<double> -- the real-space fields the term hands in
import qchem.Types;            // hmat_t<T>

export namespace qchem::BasisSet
{

//! \brief A basis that can assemble DFT potential matrices and energies by integrating real-space
//! scalar fields on its OWN mesh.  The capability the Kohn-Sham potential terms require, with the
//! integration scheme entirely the basis's business (see the design note above).
template <class T> class Band_DFT_IBS
    : public virtual Orbital_1E_IBS<T>
    , public virtual Mesh_Integrated_IBS<T>   // Overlap(f) = <i|f|j> and Integral(f) = integral f (inherited)
{
public:
    using Orbital_1E_IBS<T>::Overlap;         // keep the cached no-arg Overlap() visible beside the weighted form
    using Mesh_Integrated_IBS<T>::Overlap;    // ...and the field-operator Overlap(f), so both resolve unambiguously

    //! Coulomb repulsion matrix \f$\langle i|V_{Coul}[\rho]|j\rangle\f$ for the density \a rho (the
    //! \f$1/r_{12}\f$ Poisson solve is the basis's business), with the repulsion energy
    //! \f$E=\tfrac12\int\rho\,V_{Coul}\f$ returned by reference.  Not cached (the density changes each cycle).
    //! (The field-operator pair \c Overlap(f) / \c Integral(f) is inherited from \c Mesh_Integrated_IBS<T>;
    //! only the Poisson-kernel \c Repulsion is specific to this DFT face.)
    virtual hmat_t<T> Repulsion(const ScalarFunction<double>& rho, double& Eh) const=0;
};

} //namespace
