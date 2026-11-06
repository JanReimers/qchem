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
module;
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>
export module qchem.Mesh.SiteMultipoles;
export import qchem.Mesh;
export import qchem.Mesh.Radial;
export import qchem.Mesh.Angular;
import qchem.Parallel;   // OuterFor: every item writes only its own slots
import qchem.Types;

export namespace qchem::qcMesh
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

//! \brief The site/shell structure of an atom-centred mesh, and the transforms between point values and radial channels.
class SiteMultipoles
{
public:
    typedef std::vector<std::vector<double>> channels_t;   //!< per site, (radial node x lm) node-major, lm at index \f$l^2+l+m\f$

    //! \a centres are the mesh's site centres in block order, and \a mp the rules every site was built from.  THROWS when the mesh
    //! has no site blocks, or a point sits on no shell: the transforms mean nothing on such a mesh.
    SiteMultipoles(const Mesh& mesh, const rvec3vec_t& centres, const MeshParams& mp, int lMax)
        : SiteMultipoles(mesh, centres, MakeRadial(mp), MakeAngular(mp), lMax) {}
    SiteMultipoles(const Mesh& mesh, const rvec3vec_t& centres, const RadialMesh& radial, const AngularMesh& angular, int lMax)
        : itsPoints(mesh.Points()), itsCentres(centres), itsR(Positive(radial.R(),radial.R())), itsW(Positive(radial.R(),radial.W()))
        , itsDirs(angular.Dirs()), itsDirW(angular.W()), itsLMax(lMax), itsNLM((lMax+1)*(lMax+1)), itsLogR(MakeLog(itsR))
    {
        if (mesh.NSites()==0 || mesh.NSites()!=centres.size())
            throw std::logic_error("SiteMultipoles: needs one site block per centre; this mesh has "
                                   +std::to_string(mesh.NSites())+" blocks for "+std::to_string(centres.size())+" centres");
        if (lMax<0) throw std::invalid_argument("SiteMultipoles: lMax must not be negative");
        for (size_t i=1; i<itsR.size(); i++)
            if (!(itsR[i]>itsR[i-1])) throw std::logic_error("SiteMultipoles: the radial nodes must increase");
        MakeYlmTables();
        itsDirY.resize(itsDirs.size()*itsNLM);
        for (size_t d=0; d<itsDirs.size(); d++) Ylm(itsDirs[d],&itsDirY[d*itsNLM]);
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
                    throw std::logic_error("SiteMultipoles: a mesh point at r="+std::to_string(r)+" from its centre is on no shell of the radial rule");
                itsShell[a*itsR.size()+i].push_back(g);
            }
    }

    size_t NumSites   () const {return itsCentres.size();}
    size_t NumLM      () const {return itsNLM;}
    int    LMax       () const {return itsLMax;}
    size_t NumChannels() const {return NumSites()*itsR.size()*itsNLM;}
    //! The nodes with r > 0.  An MHL rule's first node is r = 0 with weight 0: the point at the nucleus carries no charge.
    const std::vector<double>& Radii       () const {return itsR;}
    const std::vector<double>& ShellWeights() const {return itsW;}   //!< \f$r_i^2dr_i\f$ on those nodes

    //! \f$Q_{lm}(i)=\sum_{g\in\text{shell }i}q_gY_{lm}(\Omega_g)\f$ for each site, of \a q, one entry per mesh point.
    channels_t Moments(const rvec_t& q) const
    {
        const size_t nA=NumSites(), nR=itsR.size();
        channels_t Q(nA,std::vector<double>(nR*itsNLM,0.0));
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
        return Q;
    }

    //! \f$f(r_g)=\sum_{A,lm}V^A_{lm}(|r_g-R_A|)Y_{lm}\f$: the channels of every site at every mesh point.
    rvec_t Synthesize(const channels_t& V) const
    {
        const size_t nA=NumSites();
        const ColumnSplines logR(itsLogR);
        channels_t M(nA);
        OuterFor(nA, [&](size_t a) {M[a]=logR.SecondDerivatives(V[a],itsNLM);});
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

    //! \f$\langle f|g\rangle\f$ for the expansion \a F and the field whose channels are \a G, with \f$f\f$ evaluated ONLY on each site's
    //! own shells: \f$\sum_A\sum_iW_i\sum_da_df_A(r_i,\Omega_d)\,g(r_i,\Omega_d)\f$, on complete shells and without the partition.
    //!
    //! NEVER SYNTHESIZE \a F TO THE MESH INSTEAD.  A density's site shares \f$p_A\rho\f$ have a sharp edge at the Becke boundary, so
    //! their channels are not smooth in r, and a spline to the other sites' points loses charge: 30x the mesh's own error on TeF₆ PP
    //! 80/29, in the bonds, unmoved by the angular rule, and it stalls at 1e-4 with the radial count.  That made Eee 6.4 mHa wrong.
    //! Here \f$Y_{lm}\f$ is orthonormal, so \f$f\f$ keeps its charge exactly; \f$g\f$, the smooth one, is the one that travels.
    double OwnShellOverlap(const channels_t& F, const channels_t& G) const
    {
        const size_t nA=NumSites(), nR=itsR.size(), nD=itsDirs.size();
        const ColumnSplines logR(itsLogR);
        channels_t M(nA);
        OuterFor(nA, [&](size_t a) {M[a]=logR.SecondDerivatives(G[a],itsNLM);});
        std::vector<double> perShell(nA*nR);
        OuterFor(nA*nR, [&](size_t item)
        {
            const size_t a=item/nR, i=item%nR;
            std::vector<double> y(itsNLM);
            const double* f=&F[a][i*itsNLM];
            double s=0.0;
            for (size_t d=0; d<nD; d++)
            {
                const double* yd=&itsDirY[d*itsNLM];
                double fd=0.0;
                for (size_t c=0; c<itsNLM; c++) fd+=f[c]*yd[c];
                const rvec3_t p=itsCentres[a]+itsR[i]*itsDirs[d];
                double g=0.0;
                for (size_t b=0; b<nA; b++) g+=AtPoint(logR,G[b],M[b],p-itsCentres[b],y.data());
                s+=itsDirW[d]*fd*g;
            }
            perShell[item]=itsW[i]*s;
        });
        double e=0.0;
        for (double x : perShell) e+=x;
        return e;
    }

    //! The channels as one vector, sites in order: the coefficient axis of a fit basis over them.
    rvec_t Flatten(const channels_t& Q) const
    {
        rvec_t f(NumChannels());
        size_t k=0;
        for (const auto& site : Q) for (double x : site) f[k++]=x;
        return f;
    }
    channels_t Split(const rvec_t& f) const
    {
        channels_t Q(NumSites(),std::vector<double>(itsR.size()*itsNLM));
        size_t k=0;
        for (auto& site : Q) for (double& x : site) x=f[k++];
        return Q;
    }

private:
    static constexpr size_t theBlock=256;
    //! The entries of \a x at the nodes of \a r with r > 0.
    static std::vector<double> Positive(const rvec_t& r, const rvec_t& x)
    {
        std::vector<double> s;
        for (size_t i=0; i<r.size(); i++) if (r[i]>0.0) s.push_back(x[i]);
        return s;
    }
    static std::vector<double> MakeLog(const std::vector<double>& r)
    {
        std::vector<double> s(r.size());
        for (size_t i=0; i<r.size(); i++) s[i]=std::log(r[i]);
        return s;
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
    std::vector<double> itsR, itsW;
    rvec3vec_t itsDirs;                          //!< the angular rule's directions...
    rvec_t     itsDirW;                          //!< ...its weights, which sum to \f$4\pi\f$...
    std::vector<double> itsDirY;                 //!< ...and \f$Y_{lm}\f$ at each, direction-major
    int    itsLMax;
    size_t itsNLM;
    std::vector<double> itsLogR;
    std::vector<std::vector<size_t>> itsShell;   //!< point indices per (site, shell)
    std::vector<double> itsPmm, itsA, itsB;      //!< the \f$Y_{lm}\f$ recursion's constants
};

} //namespace
