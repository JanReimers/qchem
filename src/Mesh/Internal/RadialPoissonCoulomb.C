// File: Mesh/Internal/RadialPoissonCoulomb.C  The Coulomb potential of a mesh's charges, by a radial Poisson solve per (site, l, m).
//
// Becke & Dickson 1988.  qchem.Mesh.SiteMultipoles turns the charges into each site's shell moments Q_lm(i) and carries channels
// back to the points; what is here is the solve between the two.  Each (A,l,m) channel has the 1-D potential
//     v_lm(r) = 4 pi/(2l+1) [ r^-(l+1) int_0^r h r'^l dr' + r^l int_r^inf h r'^-(l+1) dr' ],       h dr = the channel's charge,
// and v(r_g) sums the channels of every site at every point.
//
// THE RADIAL INTEGRALS MUST BE DONE THIS WAY.  Measured on the TeF₆ prototype: differencing cumulative integrals blows up by
// l = 8, splining the r^(1-l) integrand does too, and splining rho alone gets the charge wrong.  What works is ONE spline, of
// the charge per unit node index (which is what Q_lm(i) is, for any radial rule uniform in its index: MHL, Log, Linear), with the
// l-dependence applied inside each interval as ratios <= 1, by 4-point Gauss-Legendre, in the two recursions
//     U_(i+1) = (r_i/r_(i+1))^(l+1) U_i + sum_q (r_q/r_(i+1))^l / r_(i+1) b_q        (outward)
//     T_i     = (r_i/r_(i+1))^l T_(i+1) + sum_q (r_i/r_q)^l / r_q b_q                (inward)
//
// ACCURACY (same prototype): truncation in l does not depend on the rung, 4e-5 Ha at l = 8-10, 6e-6 at 12, 2e-6 at 14-16 on a PP
// density; an all-electron density wants l about 20 and about twice the radial points the contraction needs.  Past that the
// floor is the final mesh contraction, not this operator.
//
// The bodies are in Internal/Imp/RadialPoissonCoulomb.C.
module;
#include <cstddef>
#include <memory>
#include <utility>
#include <vector>
export module qchem.Mesh.Internal.RadialPoissonCoulomb;
export import qchem.Mesh.CoulombOperator;
import qchem.Mesh.SiteMultipoles;
import qchem.Types;

export namespace qchem::qcMesh
{

//! \brief The radial Poisson solve of one site's channels, in place of a density's charge channels their potential's.
class RadialPoissonSolver
{
public:
    typedef qcMesh::SiteMultipoles::channels_t channels_t;
    explicit RadialPoissonSolver(std::shared_ptr<const qcMesh::SiteMultipoles> sites);

    //! Each (site, l) channel's radial solve: charge channels \a Q in, potential channels out.
    channels_t Solve(const channels_t& Q) const;
    //! The TRANSPOSE of \c Solve: \f$\partial/\partial Q\f$ of \f$\sum W\cdot V(Q)\f$.  Each (site, l) block is a linear map, read backwards
    //! (\c SolveAdjoint of one block).
    channels_t SolveAdjoint(const channels_t& W) const;

private:
    static std::vector<double> MakeLog(const std::vector<double>& r);
    //! The charge spline's nodes: index 0 is r = 0, nR+1 is r = infinity, where the charge per unit index is zero.
    static std::vector<double> MakeIndexNodes(size_t nR);
    //! The rule's own smooth part: log r less its power at the origin, whose log t no polynomial in t follows.
    std::vector<double> MakeSmoothLog() const;
    //! r at index \a t in (0, nR]: a spline of that smooth part, and the power law inside the first node, where the rule starts.
    //! THIS MAP CARRIES THE ENERGY.  A quadratic in log r here put 95 μHa into TeF₆'s PP Hartree (2.7e-7 of E_H), at any rung.
    void RadialAt(double t, double& r) const;

    //! Per interval i (index i to i+1, r_(i-1) to r_i with r_(-1) = 0) and Gauss point k: its r_q, and its charge weights on the index
    //! spline, the Gauss weight folded in (A, B on h_i, h_(i+1); C, D on their second derivatives).  l-independent: built once.
    void MakeGaussPoints();

    //! One (site, l) block of channels: \a Q in, this l's columns of \a V out.  See the file header.
    void Solve(const std::vector<double>& Q, int l, std::vector<double>& V) const;

    //! \c Solve of one block read backwards, step by step: the potential's weights \a W in, this l's columns of the charges' \a G out.
    void SolveAdjoint(const std::vector<double>& W, int l, std::vector<double>& G) const;

    std::shared_ptr<const qcMesh::SiteMultipoles> itsSites;
    const std::vector<double>& itsR;   //!< the sites' own nodes
    int    itsLMax;
    size_t itsNLM;
    std::vector<double> itsLogR, itsIndex;
    double itsInnerPower;
    qcMesh::ColumnSplines itsNodes;              //!< the radial nodes, by index
    std::vector<double> itsSmooth, itsSmoothM;   //!< the rule's smooth part, and its spline's second derivatives
    std::vector<double> itsRq, itsGW;            //!< per (interval, Gauss point): r_q, and its four charge weights (MakeGaussPoints)
};

//! \brief The potential at the mesh's points of charges at those points: the \f$\delta\f$ route's operator.
//!
//! IT ANSWERS THE SYMMETRIZED \f$\bar K=\tfrac12(K+K^T)\f$ (user ruling 2026-11-23, "A is the only correct option").  \f$K\f$ = Synthesize o
//! Solve o Moments expands each source about its own site and reads it at another's points, so \f$K_{gh}\ne K_{hg}\f$ across sites, by
//! the l truncation and then a mesh floor.  The energy \f$\tfrac12q^TKq\f$ sees only the symmetric part, so the Fock \f$\partial E/\partial D\f$
//! needs \f$\bar Kq\f$, and so does the fixed-D force.  MEASURED with \f$Kq\f$ in their place (lane S 2026-11-23, water LDA dzvp all-electron,
//! lMax 16): the force 7e-5 Ha/bohr off its energy's finite difference, 2e-5 of the Hartree force; the synthetic charge estimate had been
//! 100x smaller.  \f$K^T\f$ = Moments\f$^T\f$ o Solve\f$^T\f$ o Synthesize\f$^T\f$, each read backwards, about the forward pass's cost again.
class RadialPoissonCoulomb : public virtual CoulombOperator
{
public:
    explicit RadialPoissonCoulomb(std::shared_ptr<const qcMesh::SiteMultipoles> sites) : itsSites(sites), itsSolver(sites) {}
    //! \copydoc CoulombOperator::Potential
    //! \f$\bar Kq=\tfrac12(Kq+K^Tq)\f$, see the class.
    rvec_t Potential(const rvec_t& q) const override;
    //! \copydoc CoulombOperator::MakeSelfRepulsionPullback
    //! Each site's moments and radial solve are rigid with it, so only the synthesis's cross-site reads move (\c SiteMultipoles::CrossSitePullback).
    //! \f$q^TdKq=q^Td\bar Kq\f$, so the forward channels serve.
    bool       CanPullbackSelfRepulsion() const override {return true;}
    rvec3vec_t MakeSelfRepulsionPullback(const rvec_t& q) const override;
private:
    std::shared_ptr<const qcMesh::SiteMultipoles> itsSites;
    RadialPoissonSolver itsSolver;
};

//! \brief The same solve on a multipole fit's own coefficient axis: charge channels in, potential channels out.
//!
//! THE SECOND LEG IS THE ONE THING THIS ROUTE ADDS OVER \f$\delta\f$ x Poisson.  \f$v\f$ and \f$J\f$ are the same numbers.  But here the
//! fit is a function, so \f$\langle\tilde\rho|\tilde\rho\rangle\f$ is available, taken on each site's own shells (see
//! \c SiteMultipoles::OwnShellOverlap), and the two-leg \f$2E_{fit}-E_{fitfit}\f$ cancels the \f$l\f$ truncation to second order.  The
//! mesh's own charge error stays, first order, so at a given rung this matches \f$\delta\f$ and does not beat it.  What it buys is a lower
//! \f$l\f$: at fixed density, all-electron TeF₆ 150/47 is -6.6e-5 at l = 8, where \f$\delta\f$ needs l = 20 for -7.0e-5.
class MultipoleCoulomb : public virtual CoulombOperator
{
public:
    //! \a metric is the fit basis's diagonal.
    MultipoleCoulomb(std::shared_ptr<const qcMesh::SiteMultipoles> sites, rvec_t metric)
        : itsSites(sites), itsSolver(sites), itsMetric(std::move(metric)) {}
    //! \copydoc CoulombOperator::Potential
    rvec_t Potential(const rvec_t& q) const override;
    //! \copydoc CoulombOperator::SelfRepulsion
    double SelfRepulsion(const rvec_t& q, const rvec_t& v) const override;
    //! \copydoc CoulombOperator::HartreeEnergy
    //! Here the legs differ, and the second is what cancels the \f$l\f$ truncation: first order in the mesh, second in \f$l\f$
    //! (accepted 2026-11-06).  \a trDV is NOT \f$q\cdot v\f$: each site's channels hold only that site's own potential, so \f$q\cdot v\f$
    //! misses every cross-site term (23.5 Ha on water).
    double HartreeEnergy(double trDV, const rvec_t& q, const rvec_t& v) const override;
    bool UsesTrDV() const override {return true;}
private:
    std::shared_ptr<const qcMesh::SiteMultipoles> itsSites;
    RadialPoissonSolver itsSolver;
    rvec_t itsMetric;
};

} //namespace
