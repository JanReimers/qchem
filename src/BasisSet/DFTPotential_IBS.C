// File: BasisSet/DFTPotential_IBS.C  Abstract DFT potential/energy assembly capability (mesh-integrated).
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
export module qchem.BasisSet.DFTPotential_IBS;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.ScalarFunction;   // ScalarFunction<double> -- the real-space fields the term hands in
import qchem.Structure;        // Structure -- the external-potential structure factor
import qchem.Types;            // hmat_t<T>

export namespace BasisSet
{

//! \brief A basis that can assemble DFT potential matrices and energies by integrating real-space
//! scalar fields on its OWN mesh.  The capability the Kohn-Sham potential terms require, with the
//! integration scheme entirely the basis's business (see the design note above).
template <class T> class DFTPotential_IBS
    : public virtual Orbital_1E_IBS<T>
{
public:
    //! Matrix \f$\langle i|V|j\rangle\f$ of a real-space scalar potential \a V (e.g. the XC term passes
    //! \f$V(r)=v_{xc}(\rho(r))\f$); the basis integrates it on its own mesh.
    virtual hmat_t<T> IntegralPotential(const ScalarFunction<double>& V) const=0;

    //! Hartree matrix \f$\langle i|V_H[\rho]|j\rangle\f$ for the density \a rho (the Poisson solve is the
    //! basis's business), with the Hartree energy \f$E_H=\tfrac12\int\rho V_H\f$ returned by reference.
    virtual hmat_t<T> IntegralHartree(const ScalarFunction<double>& rho, double& Eh) const=0;

    //! Scalar integral \f$\int f\,d^3r\f$ over the cell, on the basis's own mesh (XC energy, double-counts).
    virtual double Integral(const ScalarFunction<double>& f) const=0;

    //! The configured external (pseudo)potential matrix \f$\langle i|V_{ext}|j\rangle\f$ for this
    //! structure -- the basis owns its model (bare \f$-Z/r\f$, norm-conserving pseudopotential, ...).
    virtual hmat_t<T> MakeExternalPotential(const Structure*) const=0;
};

} //namespace
