// File: Fitting/Internal/RadialPoissonCoulomb.C  The Coulomb potential of a mesh's charges, by a radial Poisson solve per (site, l, m).
//
// Becke & Dickson 1988.  An atom-centred mesh's weights already carry the partition, so the charges on site A's block ARE that
// site's share of rho, and their moments shell by shell need no solve and no separate partition function:
//     Q_lm(i) = sum over shell i of q_g Y_lm(direction of g).
// Each (A,l,m) channel then has the 1-D potential
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
#include <numbers>
#include <stdexcept>
#include <vector>
export module qchem.Fitting.Internal.RadialPoissonCoulomb;
export import qchem.Fitting.CoulombOperator;
import qchem.Parallel;   // OuterFor: every item writes only its own slots
import qchem.Types;

export namespace qchem::Fitting
{

//! \brief Natural cubic splines through \a nCol columns on shared nodes.  Node-major storage: one node's columns are contiguous.
class ColumnSplines
{
public:
    explicit ColumnSplines(std::vector<double> x) : itsX(std::move(x)) {}
    //! Second derivatives of \a y, an (n x nCol) node-major table.  Thomas, once for every column.
    std::vector<double> SecondDerivatives(const std::vector<double>& y, size_t nCol) const
    {
        const size_t n=itsX.size();
        std::vector<double> m(n*nCol,0.0), d(n,0.0);
        for (size_t i=1; i+1<n; i++)
        {
            const double a=itsX[i]-itsX[i-1], b=itsX[i+1]-itsX[i], piv=2.0*(a+b)-(i>1 ? a*d[i-1] : 0.0);
            d[i]=b/piv;
            for (size_t c=0; c<nCol; c++)
            {
                const double rhs=6.0*((y[(i+1)*nCol+c]-y[i*nCol+c])/b-(y[i*nCol+c]-y[(i-1)*nCol+c])/a);
                m[i*nCol+c]=(rhs-(i>1 ? a*m[(i-1)*nCol+c] : 0.0))/piv;
            }
        }
        for (size_t i=n-2; i>=1; i--)
        {
            for (size_t c=0; c<nCol; c++) m[i*nCol+c]-=d[i]*m[(i+1)*nCol+c];
            if (i==1) break;
        }
        return m;
    }
    //! The interval holding \a t, clamped to the table.
    size_t Interval(double t) const
    {
        const size_t i=std::upper_bound(itsX.begin(),itsX.end(),t)-itsX.begin();
        return std::min(std::max<size_t>(i,1),itsX.size()-1)-1;
    }
    //! The four weights of the value at \a t in interval \a i: A y_i + B y_(i+1) + C m_i + D m_(i+1).
    void Weights(size_t i, double t, double& A, double& B, double& C, double& D) const
    {
        const double h=itsX[i+1]-itsX[i];
        B=(t-itsX[i])/h; A=1.0-B;
        C=(A*A*A-A)*h*h/6.0; D=(B*B*B-B)*h*h/6.0;
    }
private:
    std::vector<double> itsX;
};

class RadialPoissonCoulomb : public virtual CoulombOperator
{
public:
    //! \a centres are the mesh's site centres in block order, \a radii the nodes of the radial rule every site was built from.
    //! THROWS when the mesh has no site blocks, or a point sits on no shell: this operator is meaningless on such a mesh.
    RadialPoissonCoulomb(const qcMesh::Mesh& mesh, const rvec3vec_t& centres, const rvec_t& radii, int lMax)
        : itsPoints(mesh.Points()), itsCentres(centres), itsR(ToStd(radii)), itsLMax(lMax), itsNLM((lMax+1)*(lMax+1))
        , itsLogR(MakeLog(itsR)), itsIndex(MakeIndexNodes(itsR.size()))
    {
        if (mesh.NSites()==0 || mesh.NSites()!=centres.size())
            throw std::logic_error("RadialPoissonCoulomb: needs one site block per centre; this mesh has "
                                   +std::to_string(mesh.NSites())+" blocks for "+std::to_string(centres.size())+" centres");
        if (lMax<0) throw std::invalid_argument("RadialPoissonCoulomb: lMax must not be negative");
        for (size_t i=1; i<itsR.size(); i++)
            if (!(itsR[i]>itsR[i-1])) throw std::logic_error("RadialPoissonCoulomb: the radial nodes must increase");
        MakeYlmTables();
        // Which shell of its own site each point is on, as lists per (site, shell).
        itsShell.assign(centres.size()*itsR.size(),{});
        for (size_t a=0; a<centres.size(); a++)
            for (size_t g=mesh.SiteBegin(a); g<mesh.SiteEnd(a); g++)
            {
                const double r=norm(itsPoints[g]-centres[a]);
                if (r==0.0) continue;   // the rule's zero-weight node at the nucleus: no charge, no direction
                const size_t k=std::min<size_t>(std::lower_bound(itsR.begin(),itsR.end(),r)-itsR.begin(),itsR.size()-1);
                const size_t i=(k>0 && r-itsR[k-1]<itsR[k]-r) ? k-1 : k;
                if (std::fabs(r-itsR[i])>1e-8*itsR[i])
                    throw std::logic_error("RadialPoissonCoulomb: a mesh point at r="+std::to_string(r)+" from its centre is on no shell of the radial rule");
                itsShell[a*itsR.size()+i].push_back(g);
            }
    }

    //! \copydoc CoulombOperator::Potential
    rvec_t Potential(const rvec_t& q) const override
    {
        const size_t nA=itsCentres.size(), nR=itsR.size();
        // Step 1: the shell moments.
        std::vector<std::vector<double>> Q(nA,std::vector<double>(nR*itsNLM,0.0));
        OuterFor(nA*nR, [&](size_t item)
        {
            const size_t a=item/nR, i=item%nR;
            std::vector<double> y(itsNLM);
            double* out=&Q[a][i*itsNLM];
            for (size_t g : itsShell[item])
            {
                const rvec3_t d=itsPoints[g]-itsCentres[a];
                Ylm(d/norm(d),y.data());
                for (size_t c=0; c<itsNLM; c++) out[c]+=q[g]*y[c];
            }
        });
        // Step 2: each (site, l) channel's radial solve, then the spline that carries it to the other sites' points.
        const ColumnSplines logR(itsLogR);
        std::vector<std::vector<double>> V(nA,std::vector<double>(nR*itsNLM,0.0)), M(nA);
        OuterFor(nA*size_t(itsLMax+1), [&](size_t item)
        {
            const size_t a=item/(itsLMax+1); const int l=item%(itsLMax+1);
            Solve(Q[a],l,V[a]);
        });
        OuterFor(nA, [&](size_t a) {M[a]=logR.SecondDerivatives(V[a],itsNLM);});
        // Step 3: the channels of every site at every point.
        const size_t n=itsPoints.size(), nBlock=(n+theBlock-1)/theBlock;
        rvec_t v(n);
        OuterFor(nBlock, [&](size_t b)
        {
            std::vector<double> y(itsNLM);
            for (size_t g=b*theBlock; g<std::min(n,(b+1)*theBlock); g++)
            {
                double s=0.0;
                for (size_t a=0; a<nA; a++) s+=AtPoint(logR,V[a],M[a],itsPoints[g]-itsCentres[a],y.data());
                v[g]=s;
            }
        });
        return v;
    }

private:
    static constexpr size_t theBlock=256;
    //! The nodes with r > 0.  An MHL rule's first node is r = 0 with weight 0: the point at the nucleus carries no charge, and
    //! it is the charge spline's own index-0 end anyway.
    static std::vector<double> ToStd(const rvec_t& r)
    {
        std::vector<double> s;
        for (size_t i=0; i<r.size(); i++) if (r[i]>0.0) s.push_back(r[i]);
        return s;
    }
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
    //! r at index \a t in (0, nR]: log r splined on the nodes, and a power law inside the first node, where the rule starts.
    void RadialAt(double t, double& r) const
    {
        if (t<=1.0) {r=itsR[0]*std::pow(t,itsInnerPower); return;}
        const size_t i=std::min<size_t>(size_t(t),itsR.size()-1)-1;   // nodes i+1, i+2 in index units
        const double u=t-double(i+1);
        r=std::exp((1.0-u)*itsLogR[i]+u*itsLogR[i+1]+(u*(u-1.0)/2.0)*Curvature(i));
    }
    double Curvature(size_t i) const   // second difference of log r about the interval: the quadratic term of RadialAt
    {
        const size_t k=std::min(std::max<size_t>(i,1),itsLogR.size()-2);
        return itsLogR[k+1]-2.0*itsLogR[k]+itsLogR[k-1];
    }

    //! One (site, l) block of channels: \a Q in, this l's columns of \a V out.  See the file header.
    void Solve(const std::vector<double>& Q, int l, std::vector<double>& V) const
    {
        static const double gx[4]={-0.861136311594053,-0.339981043584856,0.339981043584856,0.861136311594053};
        static const double gw[4]={ 0.347854845137454, 0.652145154862546,0.652145154862546,0.347854845137454};
        const size_t nR=itsR.size(), nm=2*l+1, c0=size_t(l)*l;
        std::vector<double> h((nR+2)*nm,0.0);
        for (size_t i=0; i<nR; i++) for (size_t m=0; m<nm; m++) h[(i+1)*nm+m]=Q[i*itsNLM+c0+m];
        const ColumnSplines index(itsIndex);
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

    //! Every channel of one site at displacement \a d from it.  Past the table's ends the channel continues as its exact power.
    double AtPoint(const ColumnSplines& logR, const std::vector<double>& V, const std::vector<double>& M, const rvec3_t& d, double* y) const
    {
        const double r=norm(d);
        if (r==0.0) return V[0]*itsPmm[0];   // only l = 0 survives at the centre
        Ylm(d/r,y);
        const size_t nR=itsR.size();
        double s=0.0;
        if (r<=itsR[0] || r>=itsR[nR-1])
        {
            const bool inner=r<=itsR[0];
            const double* v=&V[(inner ? 0 : nR-1)*itsNLM];
            for (int l=0; l<=itsLMax; l++)
            {
                const double scale=inner ? std::pow(r/itsR[0],l) : std::pow(itsR[nR-1]/r,l+1);
                for (int c=l*l; c<(l+1)*(l+1); c++) s+=scale*v[c]*y[c];
            }
            return s;
        }
        const double t=std::log(r);
        const size_t i=logR.Interval(t);
        double A,B,C,D; logR.Weights(i,t,A,B,C,D);
        const double *v0=&V[i*itsNLM], *v1=&V[(i+1)*itsNLM], *m0=&M[i*itsNLM], *m1=&M[(i+1)*itsNLM];
        for (size_t c=0; c<itsNLM; c++) s+=(A*v0[c]+B*v1[c]+C*m0[c]+D*m1[c])*y[c];
        return s;
    }

    void MakeYlmTables()
    {
        itsInnerPower=std::log(itsR[1]/itsR[0])/std::log(2.0);
        const int L=itsLMax+1;
        itsPmm.assign(L,0.0); itsA.assign(L*L,0.0); itsB.assign(L*L,0.0);
        itsPmm[0]=std::sqrt(1.0/(4.0*std::numbers::pi));
        for (int m=1; m<L; m++) itsPmm[m]=itsPmm[m-1]*std::sqrt((2.0*m+1.0)/(2.0*m));
        for (int m=0; m<L; m++)
            for (int l=m+1; l<L; l++)
            {
                const double l2=double(l)*l, m2=double(m)*m, k2=double(l-1)*(l-1);
                itsA[l*L+m]=std::sqrt((4.0*l2-1.0)/(l2-m2));
                itsB[l*L+m]=std::sqrt((k2-m2)/(4.0*k2-1.0));
            }
    }
    //! All orthonormal real \f$Y_{lm}\f$ at unit vector \a u, at index \f$l^2+l+m\f$.  \f$\bar P_l^m/\sin^m\theta\f$ by the upward
    //! recursion, with \f$\sin^m\theta\cos m\phi\f$ and its sine partner as powers of \f$x+iy\f$, so a pole needs no division.
    void Ylm(const rvec3_t& u, double* y) const
    {
        const int L=itsLMax+1;
        double c=1.0, s=0.0;
        for (int m=0; m<L; m++)
        {
            double p0=itsPmm[m], p1=0.0;
            for (int l=m; l<L; l++)
            {
                if (l>m) {const double p=itsA[l*L+m]*(u.z*p0-itsB[l*L+m]*p1); p1=p0; p0=p;}
                if (m==0) y[l*l+l]=p0;
                else {y[l*l+l+m]=std::numbers::sqrt2*p0*c; y[l*l+l-m]=std::numbers::sqrt2*p0*s;}
            }
            const double cn=c*u.x-s*u.y; s=c*u.y+s*u.x; c=cn;
        }
    }

    rvec3vec_t itsPoints, itsCentres;
    std::vector<double> itsR;
    int    itsLMax;
    size_t itsNLM;
    std::vector<double> itsLogR, itsIndex;
    std::vector<std::vector<size_t>> itsShell;   //!< the points of (site a, shell i), at a*nR+i
    double itsInnerPower=2.0;
    std::vector<double> itsPmm, itsA, itsB;
};

} //namespace
