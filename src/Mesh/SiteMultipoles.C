// File: Mesh/SiteMultipoles.C  A field on an atom-centred mesh as radial channels per (site, l, m), and back at the points.
//
// An atom-centred mesh's weights already carry the partition, so a field weighted by them is split into per-site shares, and
// each share's moments shell by shell need no solve and no separate partition function:
//     Q_lm(i) = sum over shell i of q_g Y_lm(direction of g).
// Two clients share this: the radial Poisson solve (Becke & Dickson 1988) and the multipole fit basis.  Synthesize carries the
// channels of every site to every point, by a cubic spline in log r, at N_pts N_sites (l+1)^2 cost.
//
// PAST THE TABLE'S ENDS a channel continues as the exact power of a harmonic: r^l inside the first node, r^-(l+1) past the
// last.  That is exact for a potential.  A density or a v_xc channel is negligible at the last node, so the same rule serves.
//
// The bodies are in Imp/SiteMultipoles.C.
module;
#include <cstddef>
#include <utility>
#include <vector>
export module qchem.Mesh.SiteMultipoles;
export import qchem.Mesh;   // Mesh, RadialMesh, AngularMesh: the built meshes the transforms read
import qchem.Types;

export namespace qchem::qcMesh
{

//! \brief Cubic splines through \a nCol columns on shared nodes.  Node-major storage: one node's columns are contiguous.
//!
//! Two end conditions.  NATURAL (m = 0 at the ends) for a field read between nodes.  ZERO SLOPE (clamped, y' = 0 at the ends) for a quantity
//! the caller INTEGRATES: a clamped spline's exact integral is the trapezoid sum of its nodes plus \f$h^2(y'_a-y'_b)/12\f$, so with zero end
//! slopes it is EXACTLY the trapezoid sum, for any data.  The radial Poisson solve needs that: it re-integrates the charge per unit node index,
//! and the rule's own charge is that trapezoid sum, so each shell carries exactly its own charge, in K's columns as in its rows.  A natural
//! spline gave the innermost shell 1.13 of its charge (invisible forward, where \f$r^2\rho\to0\f$; the transpose read it, lane S 2026-11-23).
class ColumnSplines
{
public:
    enum class End {Natural, ZeroSlope};
    explicit ColumnSplines(std::vector<double> x, End end=End::Natural) : itsX(std::move(x)), itsEnd(end) {}
    //! Second derivatives of \a y, an (n x nCol) node-major table.  Thomas, once for every column.
    std::vector<double> SecondDerivatives(const std::vector<double>& y, size_t nCol) const;
    //! \brief The TRANSPOSE of \c SecondDerivatives: \f$\partial/\partial y\f$ of \f$\sum g\cdot m(y)\f$, for \a g an (n x nCol) table like \a m.
    //! \f$m=T^{-1}Ry\f$ on the solved nodes (a natural spline's ends are m = 0, so \a g there is not read), and the tridiagonal \f$T\f$ is
    //! symmetric, so the transpose is \f$R^TT^{-1}g\f$: the same solve, then the difference stencil scattered back.  CLIENT: the δ Coulomb
    //! operator's adjoint.
    std::vector<double> SecondDerivativesAdjoint(const std::vector<double>& g, size_t nCol) const;
    //! The interval holding \a t, clamped to the table.
    size_t Interval(double t) const;
    //! The four weights of the value at \a t in interval \a i: A y_i + B y_(i+1) + C m_i + D m_(i+1).
    void Weights(size_t i, double t, double& A, double& B, double& C, double& D) const;
    //! The four weights of the DERIVATIVE d/dt at \a t in interval \a i, the same combination as \c Weights.
    void DerivativeWeights(size_t i, double t, double& A, double& B, double& C, double& D) const;
private:
    //! \f$T^{-1}\f$ (Thomas, every column) on the solved rows of \a rhs: the interior for a natural spline (its end rows answer zero), every
    //! row for a clamped one.  The tridiagonal is symmetric, \f$T_{i,i+1}=T_{i+1,i}=x_{i+1}-x_i\f$.
    std::vector<double> Solve(const std::vector<double>& rhs, size_t nCol) const;
    std::vector<double> itsX;
    End itsEnd;
};

//! \brief The site/shell structure of an atom-centred mesh, and the transforms between point values and radial channels.
class SiteMultipoles
{
public:
    typedef std::vector<std::vector<double>> channels_t;   //!< per site, (radial node x lm) node-major, lm at index \f$l^2+l+m\f$

    //! \a centres are the mesh's site centres in block order, and \a radial / \a angular the rules every site was built from (a caller
    //! holding the recipe builds them with \c MakeRadial / \c MakeAngular: this class sees BUILT meshes only, V-MESH-FACTORY step 5).
    //! THROWS when the mesh has no site blocks, or a point sits on no shell: the transforms mean nothing on such a mesh.
    SiteMultipoles(const Mesh& mesh, const rvec3vec_t& centres, const RadialMesh& radial, const AngularMesh& angular, int lMax);

    size_t NumSites   () const {return itsCentres.size();}
    size_t NumLM      () const {return itsNLM;}
    int    LMax       () const {return itsLMax;}
    size_t NumChannels() const {return NumSites()*itsR.size()*itsNLM;}
    //! The nodes with r > 0.  An MHL rule's first node is r = 0 with weight 0: the point at the nucleus carries no charge.
    const std::vector<double>& Radii       () const {return itsR;}
    const std::vector<double>& ShellWeights() const {return itsW;}   //!< \f$r_i^2dr_i\f$ on those nodes

    //! Is point \a g on a shell of its site, so that a charge there has moments?  Only the rule's zero-weight node AT a nucleus is not.
    //! CLIENT: the δ operator's symmetrization, which has no column for such a point (\c RadialPoissonCoulomb::Potential).
    bool OnShell(size_t g) const {return itsOnShell[g]!=0;}

    //! \f$Q_{lm}(i)=\sum_{g\in\text{shell }i}q_gY_{lm}(\Omega_g)\f$ for each site, of \a q, one entry per mesh point.
    channels_t Moments(const rvec_t& q) const;

    //! \f$f(r_g)=\sum_{A,lm}V^A_{lm}(|r_g-R_A|)Y_{lm}\f$: the channels of every site at every mesh point.
    rvec_t Synthesize(const channels_t& V) const;

    //! \brief THE TRANSPOSE of \c Moments: \f$w_g=\sum_{lm}G^A_{lm}(i)Y_{lm}(\Omega_g)\f$ for g on shell i of site A, zero at a point on no shell.
    //! With \c SynthesizeAdjoint and the radial solve's adjoint it builds \f$K^T\f$ (the δ Coulomb operator's symmetrization, user ruling
    //! 2026-11-23: \f$H=\partial E_H/\partial D\f$ exactly).
    rvec_t MomentsAdjoint(const channels_t& G) const;

    //! \brief THE TRANSPOSE of \c Synthesize: \f$\partial/\partial V\f$ of \f$\sum_gw_gf(r_g)\f$, the channels' spline (values and second
    //! derivatives in \f$\log r\f$) and the exact powers past the table's ends read backwards.  Threaded over (site, point chunk) with a FIXED
    //! chunk count, each chunk's sums reduced in chunk order: the same bits at any thread count.
    channels_t SynthesizeAdjoint(const rvec_t& w) const;

    //! \brief THE CROSS-SITE PULLBACK of a synthesis: \f$\partial/\partial R_C\sum_h q_h\sum_{A\ne b(h)}V_A(r_h-R_A)\f$ for each site C, the
    //! channels \a V HELD and every point riding its owner b(h) (\f$r_h=R_{b(h)}+u_h\f$).
    //!
    //! WHY ONLY CROSS SITES.  Within its own site everything is rigid: its shells, the directions on them and so its moments and their
    //! radial solve move with the site, and a point's own-site channels are read at its own (fixed) offset.  Only a point read by ANOTHER
    //! site moves relative to that site, through \f$d=r_h-R_A\f$, by \f$dR_{b(h)}-dR_A\f$.  So each (A, h) pair adds \f$+q_h\nabla V_A(d)\f$
    //! to h's owner and \f$-q_h\nabla V_A(d)\f$ to A, and the sum over sites is zero (translation invariance) by construction.
    //! \f$\nabla[V_{lm}(r)Y_{lm}(\hat d)]=V'_{lm}Y_{lm}\hat d+\tfrac{V_{lm}}{r}(1-\hat d\hat d^T)\nabla_uY_{lm}\f$, the spline's own derivative in
    //! \f$\log r\f$ (and the exact powers past the table's ends).  CLIENT: \c RadialPoissonCoulomb's \c MakeSelfRepulsionPullback (lane S
    //! 2026-11-23, the δ Hartree force's operator leg).  Threaded over point blocks, each block's sums reduced in block order: the same
    //! bits at any thread count.  THROWS when a charged point belongs to no site block.
    rvec3vec_t CrossSitePullback(const channels_t& V, const rvec_t& q) const;

    //! \f$\langle f|g\rangle\f$ for the expansion \a F and the field whose channels are \a G, with \f$f\f$ evaluated ONLY on each site's
    //! own shells: \f$\sum_A\sum_iW_i\sum_da_df_A(r_i,\Omega_d)\,g(r_i,\Omega_d)\f$, on complete shells and without the partition.
    //!
    //! NEVER SYNTHESIZE \a F TO THE MESH INSTEAD.  A density's site shares \f$p_A\rho\f$ have a sharp edge at the Becke boundary, so
    //! their channels are not smooth in r, and a spline to the other sites' points loses charge: 30x the mesh's own error on TeF₆ PP
    //! 80/29, in the bonds, unmoved by the angular rule, and it stalls at 1e-4 with the radial count.  That made Eee 6.4 mHa wrong.
    //! Here \f$Y_{lm}\f$ is orthonormal, so \f$f\f$ keeps its charge exactly; \f$g\f$, the smooth one, is the one that travels.
    double OwnShellOverlap(const channels_t& F, const channels_t& G) const;

    //! The channels as one vector, sites in order: the coefficient axis of a fit basis over them.
    rvec_t     Flatten(const channels_t& Q) const;
    channels_t Split  (const rvec_t& f) const;   //!< the inverse of \c Flatten

private:
    static constexpr size_t theBlock=256;
    static constexpr size_t theChunks=16;   //!< \c SynthesizeAdjoint's (site x chunk) items, at least: fixed, so the sum order is too
    //! The entries of \a x at the nodes of \a r with r > 0.
    static std::vector<double> Positive(const rvec_t& r, const rvec_t& x);
    static std::vector<double> MakeLog(const std::vector<double>& r);

    //! Every channel of one site at displacement \a d from it.  Past the table's ends the channel continues as its exact power.
    double AtPoint(const ColumnSplines& logR, const std::vector<double>& V, const std::vector<double>& M, const rvec3_t& d, double* y) const;

    //! \c AtPoint read backwards: \a w times each weight it gives a value (\a gv) or a second derivative (\a gm), at displacement \a d.
    void AtPointAdjoint(const ColumnSplines& logR, double w, const rvec3_t& d, double* y, double* gv, double* gm) const;

    //! \f$\nabla_d\sum_{lm}V_{lm}(|d|)Y_{lm}(\hat d)\f$, \c AtPoint's gradient (see \c CrossSitePullback).  At \f$d=0\f$ it answers zero: only a site's
    //! own zero-weight nucleus node sits there, and it carries no charge.
    rvec3_t GradientAtPoint(const ColumnSplines& logR, const std::vector<double>& V, const std::vector<double>& M, const rvec3_t& d,
                            double* y, double* gu) const;

    void MakeYlmTables();
    //! All orthonormal real \f$Y_{lm}\f$ at unit vector \a u, at index \f$l^2+l+m\f$.  \f$\bar P_l^m/\sin^m\theta\f$ by the upward
    //! recursion, with \f$\sin^m\theta\cos m\phi\f$ and its sine partner as powers of \f$x+iy\f$, so a pole needs no division.
    void Ylm(const rvec3_t& u, double* y) const;

    //! \c Ylm and the gradient of the POLYNOMIAL it evaluates, \f$\bar P_l^m(z)/\sin^m\theta\f$ times the powers of \f$x+iy\f$, with x, y, z
    //! independent: \a gu holds \f$(\partial_x,\partial_y,\partial_z)\f$ at index \f$3(l^2+l+m)\f$.  Off the sphere the polynomial is not \f$Y_{lm}\f$, so
    //! a caller keeps only the tangential part (\c GradientAtPoint).
    void YlmGradient(const rvec3_t& u, double* y, double* gu) const;

    rvec3vec_t itsPoints, itsCentres;
    std::vector<double> itsR, itsW;
    rvec3vec_t itsDirs;                          //!< the angular rule's directions...
    rvec_t     itsDirW;                          //!< ...its weights, which sum to \f$4\pi\f$...
    std::vector<double> itsDirY;                 //!< ...and \f$Y_{lm}\f$ at each, direction-major
    int    itsLMax;
    size_t itsNLM;
    std::vector<double> itsLogR;
    std::vector<std::vector<size_t>> itsShell;   //!< point indices per (site, shell)
    std::vector<size_t> itsOwner;                //!< each point's site block (\c NumSites() for none)
    std::vector<char>   itsOnShell;              //!< 1 where a point is on a shell of its site (\c OnShell)
    std::vector<double> itsPmm, itsA, itsB;      //!< the \f$Y_{lm}\f$ recursion's constants
};

} //namespace
