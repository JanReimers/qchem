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
module;
#include <algorithm>
#include <cmath>
#include <memory>
#include <numbers>
#include <stdexcept>
#include <vector>
export module qchem.Mesh.Internal.RadialPoissonCoulomb;
export import qchem.Mesh.CoulombOperator;
import qchem.Mesh.SiteMultipoles;
import qchem.Parallel;   // OuterFor: every item writes only its own slots
import qchem.Blaze;
import qchem.Types;

export namespace qchem::qcMesh
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
    {
        MakeGaussPoints();
    }

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
    //! The TRANSPOSE of \c Solve: \f$\partial/\partial Q\f$ of \f$\sum W\cdot V(Q)\f$.  Each (site, l) block is a linear map, read backwards
    //! (\c SolveAdjoint of one block).
    channels_t SolveAdjoint(const channels_t& W) const
    {
        const size_t nA=W.size(), nR=itsR.size();
        channels_t G(nA,std::vector<double>(nR*itsNLM,0.0));
        OuterFor(nA*size_t(itsLMax+1), [&](size_t item)
        {
            const size_t a=item/(itsLMax+1); const int l=item%(itsLMax+1);
            SolveAdjoint(W[a],l,G[a]);
        });
        return G;
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

    //! Per interval i (index i to i+1, r_(i-1) to r_i with r_(-1) = 0) and Gauss point k: its r_q, and its charge weights on the index
    //! spline, the Gauss weight folded in (A, B on h_i, h_(i+1); C, D on their second derivatives).  l-independent: built once.
    void MakeGaussPoints()
    {
        static const double gx[4]={-0.861136311594053,-0.339981043584856,0.339981043584856,0.861136311594053};
        static const double gw[4]={ 0.347854845137454, 0.652145154862546,0.652145154862546,0.347854845137454};
        const size_t nR=itsR.size();
        const qcMesh::ColumnSplines index(itsIndex);
        itsRq.resize(nR*4); itsGW.resize(nR*4*4);
        for (size_t i=0; i<nR; i++)
            for (int k=0; k<4; k++)
            {
                const double t=double(i)+0.5*(1.0+gx[k]);
                RadialAt(t,itsRq[i*4+k]);
                double* w=&itsGW[(i*4+k)*4];
                index.Weights(i,t,w[0],w[1],w[2],w[3]);
                for (int j=0; j<4; j++) w[j]*=0.5*gw[k];
            }
    }

    //! One (site, l) block of channels: \a Q in, this l's columns of \a V out.  See the file header.
    void Solve(const std::vector<double>& Q, int l, std::vector<double>& V) const
    {
        const size_t nR=itsR.size(), nm=2*l+1, c0=size_t(l)*l;
        std::vector<double> h((nR+2)*nm,0.0);
        for (size_t i=0; i<nR; i++) for (size_t m=0; m<nm; m++) h[(i+1)*nm+m]=Q[i*itsNLM+c0+m];
        const qcMesh::ColumnSplines index(itsIndex,qcMesh::ColumnSplines::End::ZeroSlope);   // see ColumnSplines: each shell keeps its charge
        const std::vector<double> hm=index.SecondDerivatives(h,nm);
        // The charge b_q at each Gauss point.
        const std::vector<double>& rq=itsRq;
        std::vector<double> b(nR*4*nm);
        for (size_t i=0; i<nR; i++)
            for (int k=0; k<4; k++)
            {
                const double* w=&itsGW[(i*4+k)*4];
                for (size_t m=0; m<nm; m++)
                    b[(i*4+k)*nm+m]=w[0]*h[i*nm+m]+w[1]*h[(i+1)*nm+m]+w[2]*hm[i*nm+m]+w[3]*hm[(i+1)*nm+m];
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
        // The last index interval (nR to nR+1: r_(nR-1) to infinity) holds the charge spline's tail past the last node: half the last
        // shell's charge, by the trapezoid's end weight.  It is placed AT the last node, so it is counted, and its potential is exact from
        // there in (the same r^-(l+1) Synthesize continues with past the table).  Dropping it was harmless forward (rho is ~0 there) but
        // halved the TRANSPOSE's potential at the outermost shells (0.086 Ha on a Gaussian, lane S 2026-11-23).  The cubic's exact integral.
        const double rLast=itsR[nR-1];
        for (size_t m=0; m<nm; m++) T[(nR-1)*nm+m]=(0.5*h[nR*nm+m]-(hm[nR*nm+m]+hm[(nR+1)*nm+m])/24.0)/rLast;
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

    //! \c Solve of one block read backwards, step by step: the potential's weights \a W in, this l's columns of the charges' \a G out.
    void SolveAdjoint(const std::vector<double>& W, int l, std::vector<double>& G) const
    {
        const size_t nR=itsR.size(), nm=2*l+1, c0=size_t(l)*l;
        const std::vector<double>& rq=itsRq;
        const double f=4.0*std::numbers::pi/(2*l+1);
        std::vector<double> gU(nR*nm), gT(nR*nm), gb(nR*4*nm,0.0);
        for (size_t i=0; i<nR; i++) for (size_t m=0; m<nm; m++) gU[i*nm+m]=gT[i*nm+m]=f*W[i*itsNLM+c0+m];
        for (size_t i=1; i<nR; i++)   // the inward recursion, backwards: T_(i-1) is final before T_i takes its share
        {
            const double lo=itsR[i-1], damp=std::pow(lo/itsR[i],l);
            for (int k=0; k<4; k++)
            {
                const double e=std::pow(lo/rq[i*4+k],l)/rq[i*4+k];
                for (size_t m=0; m<nm; m++) gb[(i*4+k)*nm+m]+=e*gT[(i-1)*nm+m];
            }
            for (size_t m=0; m<nm; m++) gT[i*nm+m]+=damp*gT[(i-1)*nm+m];
        }
        for (size_t i=nR; i-- >0;)   // the outward recursion, backwards
        {
            const double ri=itsR[i], lo=i>0 ? itsR[i-1] : 0.0, damp=std::pow(lo/ri,l+1);
            for (int k=0; k<4; k++)
            {
                const double e=std::pow(rq[i*4+k]/ri,l)/ri;
                for (size_t m=0; m<nm; m++) gb[(i*4+k)*nm+m]+=e*gU[i*nm+m];
            }
            if (i>0) for (size_t m=0; m<nm; m++) gU[(i-1)*nm+m]+=damp*gU[i*nm+m];
        }
        std::vector<double> gh((nR+2)*nm,0.0), ghm((nR+2)*nm,0.0);
        const double rLast=itsR[nR-1];   // the tail past the last node (see Solve)
        for (size_t m=0; m<nm; m++)
        {
            const double g=gT[(nR-1)*nm+m]/rLast;
            gh[nR*nm+m]+=0.5*g; ghm[nR*nm+m]-=g/24.0; ghm[(nR+1)*nm+m]-=g/24.0;
        }
        for (size_t i=0; i<nR; i++)
            for (int k=0; k<4; k++)
            {
                const double* w=&itsGW[(i*4+k)*4];
                for (size_t m=0; m<nm; m++)
                {
                    const double g=gb[(i*4+k)*nm+m];
                    gh[i*nm+m]+=w[0]*g; gh[(i+1)*nm+m]+=w[1]*g; ghm[i*nm+m]+=w[2]*g; ghm[(i+1)*nm+m]+=w[3]*g;
                }
            }
        const std::vector<double> fromM=qcMesh::ColumnSplines(itsIndex,qcMesh::ColumnSplines::End::ZeroSlope).SecondDerivativesAdjoint(ghm,nm);
        for (size_t i=0; i<nR; i++) for (size_t m=0; m<nm; m++) G[i*itsNLM+c0+m]=gh[(i+1)*nm+m]+fromM[(i+1)*nm+m];
    }

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
    rvec_t Potential(const rvec_t& q) const override
    {
        const rvec_t Kq =itsSites->Synthesize(itsSolver.Solve(itsSites->Moments(q)));
        const rvec_t KTq=itsSites->MomentsAdjoint(itsSolver.SolveAdjoint(itsSites->SynthesizeAdjoint(q)));
        rvec_t v(Kq.size());
        // A point on no shell (the rule's zero-weight node at a nucleus) carries no charge, so K has no column for it: its potential is Kq.
        for (size_t g=0; g<v.size(); g++) v[g]=itsSites->OnShell(g) ? 0.5*(Kq[g]+KTq[g]) : Kq[g];
        return v;
    }
    //! \copydoc CoulombOperator::MakeSelfRepulsionPullback
    //! Each site's moments and radial solve are rigid with it, so only the synthesis's cross-site reads move (\c SiteMultipoles::CrossSitePullback).
    //! \f$q^TdKq=q^Td\bar Kq\f$, so the forward channels serve.
    bool       CanPullbackSelfRepulsion() const override {return true;}
    rvec3vec_t MakeSelfRepulsionPullback(const rvec_t& q) const override {return itsSites->CrossSitePullback(itsSolver.Solve(itsSites->Moments(q)), q);}
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
    rvec_t Potential(const rvec_t& q) const override {return itsSites->Flatten(itsSolver.Solve(itsSites->Split(q)));}
    //! \copydoc CoulombOperator::SelfRepulsion
    double SelfRepulsion(const rvec_t& q, const rvec_t& v) const override
    {
        rvec_t c(q.size());
        for (size_t a=0; a<c.size(); a++) c[a]=q[a]/itsMetric[a];
        return itsSites->OwnShellOverlap(itsSites->Split(c),itsSites->Split(v));
    }
    //! \copydoc CoulombOperator::HartreeEnergy
    //! Here the legs differ, and the second is what cancels the \f$l\f$ truncation: first order in the mesh, second in \f$l\f$
    //! (accepted 2026-11-06).  \a trDV is NOT \f$q\cdot v\f$: each site's channels hold only that site's own potential, so \f$q\cdot v\f$
    //! misses every cross-site term (23.5 Ha on water).
    double HartreeEnergy(double trDV, const rvec_t& q, const rvec_t& v) const override {return trDV-0.5*SelfRepulsion(q,v);}
    bool UsesTrDV() const override {return true;}
private:
    std::shared_ptr<const qcMesh::SiteMultipoles> itsSites;
    RadialPoissonSolver itsSolver;
    rvec_t itsMetric;
};

} //namespace
