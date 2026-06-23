// File: BasisSet/Lattice_3D/LocalPotential.C  One-body LOCAL external potentials for plane-wave work.
//
// A local external potential, in a plane-wave basis, is fully specified by its reciprocal-space radial
// form factor v(|G|^2) per nuclear species.  The plane-wave basis set folds in the cell volume 1/Omega,
// the structure factor Sum_a e^{-iG.tau_a}, and the G=0 handling (see PlaneWave_IBS::MakeLocalPotential);
// a LocalPotential supplies ONLY the species form factor.
//
// This is the open/closed extension point for "lineage A" (plane waves + pseudopotentials, see
// doc/PlaneWavePlan.md): the bare Coulomb nucleus, the Gaussian-smeared nucleus (rung-1 local
// pseudopotential), and -- later -- tabulated norm-conserving pseudopotential form factors are all
// just implementations of FormFactor(); the assembler never changes.
module;
#include <cmath>

export module qchem.BasisSet.Lattice_3D.LocalPotential;
import qchem.Math; // FourPi

export namespace BasisSet::Lattice_3D
{

//! \brief A one-body local external potential, defined by its reciprocal-space form factor.
class LocalPotential
{
public:
    virtual ~LocalPotential() {}
    //! \brief \f$v(Z,|G|^2)\f$: the Fourier transform of a single atom's potential (\f$|G|>0\f$),
    //! excluding the \f$1/\Omega\f$ and the structure-factor phase (those are geometry, applied by
    //! the basis set).  [energy x volume]
    virtual double FormFactor(int Z, double G2) const=0;
};

//! \brief Bare nuclear Coulomb \f$v(G) = -4\pi Z/G^2\f$.  Physically exact but the 1s cusp makes the
//! plane-wave energy converge very slowly in \f$E_{cut}\f$.
class BareCoulomb : public LocalPotential
{
public:
    virtual double FormFactor(int Z, double G2) const {return -FourPi*Z/G2;}
};

//! \brief Gaussian-smeared nucleus (rung-1 "local pseudopotential"): the point charge is spread into a
//! Gaussian of width \f$\sigma\f$, giving \f$v(G) = -4\pi Z\,e^{-\sigma^2 G^2/2}/G^2\f$.  The smooth
//! charge removes the cusp, so the plane-wave energy converges rapidly with \f$E_{cut}\f$;
//! \f$\sigma\to 0\f$ recovers BareCoulomb.
class GaussianSmearedNucleus : public LocalPotential
{
public:
    explicit GaussianSmearedNucleus(double sigma) : itsSigma(sigma) {}
    virtual double FormFactor(int Z, double G2) const
    {
        return -FourPi*Z*std::exp(-0.5*itsSigma*itsSigma*G2)/G2;
    }
private:
    double itsSigma; //!< Smearing width (Bohr).
};

} //namespace
