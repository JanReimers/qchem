// File: Fitting/Internal/RadialPoissonCoulomb.C  The Coulomb potential of a mesh's charges, by a radial Poisson solve per (site, l, m).
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
module;
#include <algorithm>
#include <cmath>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <vector>
export module qchem.Fitting.Internal.RadialPoissonCoulomb;
export import qchem.Fitting.CoulombOperator;
import qchem.Mesh.SiteMultipoles;
import qchem.Parallel;   // OuterFor: every item writes only its own slots
import qchem.Types;

export namespace qchem::Fitting
{

//! \brief The radial Poisson solve of one site's channels, in place of a density's charge channels their potential's.
class RadialPoissonSolver
{
public:
    typedef qcMesh::SiteMultipoles::channels_t channels_t;
    explicit RadialPoissonSolver(std::shared_ptr<const qcMesh::SiteMultipoles> sites)
        : itsSites(std::move(sites)), itsR(itsSites->Radii()), itsLMax(itsSites->LMax()), itsNLM(itsSites->NumLM())
        , itsLogR(MakeLog(itsR)), itsIndex(MakeIndexNodes(itsR.size())), itsInnerPower(std::log(itsR[1]/itsR[0])/std::log(2.0))
        , itsNodes(std::vector<double>(itsIndex.begin()+1,itsIndex.end()-1)), itsSmooth(MakeSmoothLog()), itsSmoothM(itsNodes.SecondDerivatives(itsSmooth,1))
    {}

    //! Each (site, l) channel's radial solve: charge channels \a Q in, potential channels out.
    channels_t Solve(const channels_t& Q) const
    {
        const size_t nA=Q.size(), nR=itsR.size();
        channels_t V(nA,std::vector<double>(nR*itsNLM,0.0));
        OuterFor(nA*size_t(itsLMax+1), [&](size_t item)
        {
            const size_t a=item/(itsLMax+1); const int l=item%(itsLMax+1);
            Solve(Q[a],l,V[a]);
        });
        return V;
    }

private:
    static std::vector<double> MakeLog(const std::vector<double>& r)
    {
        std::vector<double> s(r.size());
        for (size_t i=0; i<r.size(); i++) s[i]=std::log(r[i]);
        return s;
    }
    //! The charge spline's nodes: index 0 is r = 0, nR+1 is r = infinity, where the charge per unit index is zero.
    static std::vector<double> MakeIndexNodes(size_t nR)
    {
        std::vector<double> t(nR+2);
        for (size_t i=0; i<nR+2; i++) t[i]=double(i);
        return t;
    }
    //! The rule's own smooth part: log r less its power at the origin, whose log t no polynomial in t follows.
    std::vector<double> MakeSmoothLog() const
    {
        std::vector<double> f(itsR.size());
        for (size_t i=0; i<f.size(); i++) f[i]=itsLogR[i]-itsInnerPower*std::log(double(i+1));
        return f;
    }
    //! r at index \a t in (0, nR]: a spline of that smooth part, and the power law inside the first node, where the rule starts.
    //! THIS MAP CARRIES THE ENERGY.  A quadratic in log r here put 95 μHa into TeF₆'s PP Hartree (2.7e-7 of E_H), at any rung.
    void RadialAt(double t, double& r) const
    {
        if (t<=1.0) {r=itsR[0]*std::pow(t,itsInnerPower); return;}
        double A,B,C,D;
        const size_t i=itsNodes.Interval(t);
        itsNodes.Weights(i,t,A,B,C,D);
        r=std::exp(itsInnerPower*std::log(t)+A*itsSmooth[i]+B*itsSmooth[i+1]+C*itsSmoothM[i]+D*itsSmoothM[i+1]);
    }

    //! One (site, l) block of channels: \a Q in, this l's columns of \a V out.  See the file header.
    void Solve(const std::vector<double>& Q, int l, std::vector<double>& V) const
    {
        static const double gx[4]={-0.861136311594053,-0.339981043584856,0.339981043584856,0.861136311594053};
        static const double gw[4]={ 0.347854845137454, 0.652145154862546,0.652145154862546,0.347854845137454};
        const size_t nR=itsR.size(), nm=2*l+1, c0=size_t(l)*l;
        std::vector<double> h((nR+2)*nm,0.0);
        for (size_t i=0; i<nR; i++) for (size_t m=0; m<nm; m++) h[(i+1)*nm+m]=Q[i*itsNLM+c0+m];
        const qcMesh::ColumnSplines index(itsIndex);
        const std::vector<double> hm=index.SecondDerivatives(h,nm);
        // Per interval i (index i to i+1, r_(i-1) to r_i with r_(-1) = 0): the charge b_q at each Gauss point, and its r_q.
        std::vector<double> b(nR*4*nm), rq(nR*4);
        for (size_t i=0; i<nR; i++)
            for (int k=0; k<4; k++)
            {
                const double t=double(i)+0.5*(1.0+gx[k]);
                RadialAt(t,rq[i*4+k]);
                double A,B,C,D; index.Weights(i,t,A,B,C,D);
                for (size_t m=0; m<nm; m++)
                    b[(i*4+k)*nm+m]=0.5*gw[k]*(A*h[i*nm+m]+B*h[(i+1)*nm+m]+C*hm[i*nm+m]+D*hm[(i+1)*nm+m]);
            }
        std::vector<double> U(nR*nm,0.0), T(nR*nm,0.0);
        for (size_t i=0; i<nR; i++)
        {
            const double ri=itsR[i], lo=i>0 ? itsR[i-1] : 0.0, damp=std::pow(lo/ri,l+1);
            for (size_t m=0; m<nm; m++)
            {
                double s=i>0 ? damp*U[(i-1)*nm+m] : 0.0;
                for (int k=0; k<4; k++) s+=std::pow(rq[i*4+k]/ri,l)/ri*b[(i*4+k)*nm+m];
                U[i*nm+m]=s;
            }
        }
        for (size_t i=nR-1; i>=1; i--)   // interval i runs r_(i-1) to r_i, so T at node i-1 takes it in
        {
            const double lo=itsR[i-1], damp=std::pow(lo/itsR[i],l);
            for (size_t m=0; m<nm; m++)
            {
                double s=damp*T[i*nm+m];
                for (int k=0; k<4; k++) s+=std::pow(lo/rq[i*4+k],l)/rq[i*4+k]*b[(i*4+k)*nm+m];
                T[(i-1)*nm+m]=s;
            }
        }
        const double f=4.0*std::numbers::pi/(2*l+1);
        for (size_t i=0; i<nR; i++) for (size_t m=0; m<nm; m++) V[i*itsNLM+c0+m]=f*(U[i*nm+m]+T[i*nm+m]);
    }

    std::shared_ptr<const qcMesh::SiteMultipoles> itsSites;
    const std::vector<double>& itsR;   //!< the sites' own nodes
    int    itsLMax;
    size_t itsNLM;
    std::vector<double> itsLogR, itsIndex;
    double itsInnerPower;
    qcMesh::ColumnSplines itsNodes;              //!< the radial nodes, by index
    std::vector<double> itsSmooth, itsSmoothM;   //!< the rule's smooth part, and its spline's second derivatives
};

//! \brief The potential at the mesh's points of charges at those points: the \f$\delta\f$ route's operator.
class RadialPoissonCoulomb : public virtual CoulombOperator
{
public:
    explicit RadialPoissonCoulomb(std::shared_ptr<const qcMesh::SiteMultipoles> sites) : itsSites(sites), itsSolver(sites) {}
    //! \copydoc CoulombOperator::Potential
    rvec_t Potential(const rvec_t& q) const override {return itsSites->Synthesize(itsSolver.Solve(itsSites->Moments(q)));}
private:
    std::shared_ptr<const qcMesh::SiteMultipoles> itsSites;
    RadialPoissonSolver itsSolver;
};

//! \brief The same solve on a multipole fit's own coefficient axis: charge channels in, potential channels out.
//!
//! THE SECOND LEG IS THE ONE THING THIS ROUTE ADDS OVER \f$\delta\f$ x Poisson.  \f$v\f$ and \f$J\f$ are the same numbers.  But here the
//! fit is a function, so \f$\langle\tilde\rho|\tilde\rho\rangle\f$ is available, taken on each site's own shells (see
//! \c SiteMultipoles::OwnShellOverlap), and FittedVee's \f$2E_{fit}-E_{fitfit}\f$ cancels the \f$l\f$ truncation to second order.  The
//! mesh's own charge error stays, first order, so at a given rung this matches \f$\delta\f$ and does not beat it.  What it buys is a lower
//! \f$l\f$: at fixed density, all-electron TeF₆ 150/47 is -6.6e-5 at l = 8, where \f$\delta\f$ needs l = 20 for -7.0e-5.
class MultipoleCoulomb : public virtual CoulombOperator
{
public:
    //! \a metric is the fit basis's diagonal.
    MultipoleCoulomb(std::shared_ptr<const qcMesh::SiteMultipoles> sites, rvec_t metric)
        : itsSites(sites), itsSolver(sites), itsMetric(std::move(metric)) {}
    //! \copydoc CoulombOperator::Potential
    rvec_t Potential(const rvec_t& q) const override {return itsSites->Flatten(itsSolver.Solve(itsSites->Split(q)));}
    //! \copydoc CoulombOperator::SelfRepulsion
    double SelfRepulsion(const rvec_t& q, const rvec_t& v) const override
    {
        rvec_t c(q.size());
        for (size_t a=0; a<c.size(); a++) c[a]=q[a]/itsMetric[a];
        return itsSites->OwnShellOverlap(itsSites->Split(c),itsSites->Split(v));
    }
private:
    std::shared_ptr<const qcMesh::SiteMultipoles> itsSites;
    RadialPoissonSolver itsSolver;
    rvec_t itsMetric;
};

} //namespace
