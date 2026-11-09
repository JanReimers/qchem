// File: BasisSet/GSpaceDensityProjector.C  A density's ρ̃(G) on a plane-wave fit basis: the OVERLAP projector over {G}.
//
// V-ONE-OF-N 3c (doc/Records/VOneOfN_D2Design.md §2a row 10, §3): what every G-space client of a density asks for -- the Hartree
// fitter (GSpaceDensityFitter), the ρ̃ mixers, the XC and local-PP terms -- is the density's overlap projections on an orthonormal
// {G} basis, ρ̃(G).  Before 3c each asked the density through its own face (FourierDensity::GetFourierDensity), a second entry
// beside tChargeDensity::Project.  Now the density projects itself onto THIS object and nobody names a Fourier face:
//   - a density with a matrix contracts its D into the per-block OVERLAP forward (GSpacePairProjectionIntegrator);
//   - a density with no D is asked for its field: its own Fourier components when it knows them (AnalyticFourier), else sampled on
//     the raster and transformed;
//   - a k-composite star-averages the sum once (StarAverage), with the ops it owns.
// No kernel anywhere: the Coulomb kernel is the Hartree physics round's (GSpaceDensityFitter::Repulsion), never a projection's.
module;
#include <memory>
#include <string>
export module qchem.BasisSet.GSpaceDensityProjector;
export import qchem.BasisSet.PairProjector;     // what the factory returns: the one projection face a density answers
export import qchem.BasisSet.Orbital_DFT_IBS;   // cFIT_SF_ABS (the fit basis it projects onto)

export namespace qchem::BasisSet
{

//! \brief The ρ̃ projector on the plane-wave fit basis \a sf: the overlap projection on its orthonormal {G}, with the k-composite's
//! star-average timed under \a starAverageBucket (each client keeps its own ledger row).  Clients hold the FACE: they only ever ask
//! \c density.Project(*p).  THROWS when \a sf is not orthonormal or has no raster.
std::unique_ptr<PairProjector<dcmplx>> MakeGSpaceDensityProjector(std::shared_ptr<const cFIT_SF_ABS> sf,
                                                                  std::string starAverageBucket="scf: rho-tilde IBZ star-average (SymmetrizeGMap)");

} //namespace
