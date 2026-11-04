// File: ChargeDensity/DensityMatrixDistance.C  How far apart two density matrices are, in the basis's own metric.
//
// A capability face, not part of tMixableDensity: only an SCF whose Hamiltonian owns no real-space quadrature asks, and the
// densities that can answer are the matrix-backed ones (an irrep leaf, a composite of them).  The SCF cross-casts to it.
module;
export module qchem.ChargeDensity.DensityMatrixDistance;
export import qchem.ChargeDensity;

export namespace qchem::ChargeDensity
{

//! \brief \f$\mathrm{Tr}(\Delta D\,S\,\Delta D\,S)\f$: the squared S-metric size of a density-matrix change, summed over blocks and spins.
//!
//! It equals the squared Frobenius norm of \f$\Delta D\f$ in an orthonormal basis, so it does not depend on the basis the
//! density is written in: this is the invariant form of the RMS(\f$\Delta D\f$) other codes print.
template <class T> class tDensityMatrixDistance
{
public:
    virtual ~tDensityMatrixDistance() = default;
    //! \a cd must be the same representation: the same blocks over the same basis.
    virtual double SquaredDistanceTo(const tMixableDensity<T>& cd) const=0;
};

} //namespace
