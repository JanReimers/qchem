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
    std::vector<double> SecondDerivatives(const std::vector<double>& y, size_t nCol) const
    {
        const size_t n=itsX.size();
        std::vector<double> rhs(n*nCol,0.0);
        for (size_t i=1; i+1<n; i++)
        {
            const double a=itsX[i]-itsX[i-1], b=itsX[i+1]-itsX[i];
            for (size_t c=0; c<nCol; c++) rhs[i*nCol+c]=6.0*((y[(i+1)*nCol+c]-y[i*nCol+c])/b-(y[i*nCol+c]-y[(i-1)*nCol+c])/a);
        }
        if (itsEnd==End::ZeroSlope)
        {
            const double a=itsX[1]-itsX[0], b=itsX[n-1]-itsX[n-2];
            for (size_t c=0; c<nCol; c++)
            {
                rhs[c]=6.0*(y[nCol+c]-y[c])/a;
                rhs[(n-1)*nCol+c]=-6.0*(y[(n-1)*nCol+c]-y[(n-2)*nCol+c])/b;
            }
        }
        return Solve(rhs,nCol);
    }
    //! \brief The TRANSPOSE of \c SecondDerivatives: \f$\partial/\partial y\f$ of \f$\sum g\cdot m(y)\f$, for \a g an (n x nCol) table like \a m.
    //! \f$m=T^{-1}Ry\f$ on the solved nodes (a natural spline's ends are m = 0, so \a g there is not read), and the tridiagonal \f$T\f$ is
    //! symmetric, so the transpose is \f$R^TT^{-1}g\f$: the same solve, then the difference stencil scattered back.  CLIENT: the δ Coulomb
    //! operator's adjoint.
    std::vector<double> SecondDerivativesAdjoint(const std::vector<double>& g, size_t nCol) const
    {
        const size_t n=itsX.size();
        std::vector<double> rhs(g);
        if (itsEnd==End::Natural) for (size_t c=0; c<nCol; c++) {rhs[c]=0.0; rhs[(n-1)*nCol+c]=0.0;}
        const std::vector<double> z=Solve(rhs,nCol);
        std::vector<double> gy(n*nCol,0.0);
        for (size_t i=1; i+1<n; i++)
        {
            const double a=itsX[i]-itsX[i-1], b=itsX[i+1]-itsX[i];
            for (size_t c=0; c<nCol; c++)
            {
                const double s=6.0*z[i*nCol+c];
                gy[(i+1)*nCol+c]+=s/b; gy[i*nCol+c]-=s/b+s/a; gy[(i-1)*nCol+c]+=s/a;
            }
        }
        if (itsEnd==End::ZeroSlope)
        {
            const double a=itsX[1]-itsX[0], b=itsX[n-1]-itsX[n-2];
            for (size_t c=0; c<nCol; c++)
            {
                const double s0=6.0*z[c]/a, s1=6.0*z[(n-1)*nCol+c]/b;
                gy[nCol+c]+=s0; gy[c]-=s0;
                gy[(n-1)*nCol+c]-=s1; gy[(n-2)*nCol+c]+=s1;
            }
        }
        return gy;
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
    //! The four weights of the DERIVATIVE d/dt at \a t in interval \a i, the same combination as \c Weights.
    void DerivativeWeights(size_t i, double t, double& A, double& B, double& C, double& D) const
    {
        const double h=itsX[i+1]-itsX[i], b=(t-itsX[i])/h, a=1.0-b;
        A=-1.0/h; B=1.0/h;
        C=-(3.0*a*a-1.0)*h/6.0; D=(3.0*b*b-1.0)*h/6.0;
    }
private:
    //! \f$T^{-1}\f$ (Thomas, every column) on the solved rows of \a rhs: the interior for a natural spline (its end rows answer zero), every
    //! row for a clamped one.  The tridiagonal is symmetric, \f$T_{i,i+1}=T_{i+1,i}=x_{i+1}-x_i\f$.
    std::vector<double> Solve(const std::vector<double>& rhs, size_t nCol) const
    {
        const size_t n=itsX.size(), lo=itsEnd==End::Natural ? 1 : 0, hi=itsEnd==End::Natural ? n-2 : n-1;
        auto h=[&](size_t i) {return itsX[i+1]-itsX[i];};
        std::vector<double> m(n*nCol,0.0), d(n,0.0);
        for (size_t i=lo; i<=hi; i++)
        {
            const double a=i>0 ? h(i-1) : 0.0, b=i+1<n ? h(i) : 0.0, piv=2.0*(a+b)-(i>lo ? a*d[i-1] : 0.0);
            d[i]=b/piv;
            for (size_t c=0; c<nCol; c++) m[i*nCol+c]=(rhs[i*nCol+c]-(i>lo ? a*m[(i-1)*nCol+c] : 0.0))/piv;
        }
        for (size_t i=hi; i-- >lo;)
            for (size_t c=0; c<nCol; c++) m[i*nCol+c]-=d[i]*m[(i+1)*nCol+c];
        return m;
    }
    std::vector<double> itsX;
    End itsEnd;
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
        itsOnShell.assign(itsPoints.size(),0);
        for (const auto& shell : itsShell) for (size_t g : shell) itsOnShell[g]=1;
        itsOwner.assign(itsPoints.size(),centres.size());   // NumSites() = no owner
        for (size_t a=0; a<centres.size(); a++)
            for (size_t g=mesh.SiteBegin(a); g<mesh.SiteEnd(a); g++) itsOwner[g]=a;
    }

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

    //! \brief THE TRANSPOSE of \c Moments: \f$w_g=\sum_{lm}G^A_{lm}(i)Y_{lm}(\Omega_g)\f$ for g on shell i of site A, zero at a point on no shell.
    //! With \c SynthesizeAdjoint and the radial solve's adjoint it builds \f$K^T\f$ (the δ Coulomb operator's symmetrization, user ruling
    //! 2026-11-23: \f$H=\partial E_H/\partial D\f$ exactly).
    rvec_t MomentsAdjoint(const channels_t& G) const
    {
        const size_t nA=NumSites(), nR=itsR.size();
        rvec_t w(itsPoints.size());
        for (size_t g=0; g<w.size(); g++) w[g]=0.0;
        OuterFor(nA*nR, [&](size_t item)   // each point is on one shell: every item writes only its own points
        {
            const size_t a=item/nR, i=item%nR;
            std::vector<double> y(itsNLM);
            const double* in=&G[a][i*itsNLM];
            for (size_t g : itsShell[item])
            {
                const rvec3_t d=itsPoints[g]-itsCentres[a];
                Ylm(d/norm(d),y.data());
                double s=0.0;
                for (size_t c=0; c<itsNLM; c++) s+=in[c]*y[c];
                w[g]=s;
            }
        });
        return w;
    }

    //! \brief THE TRANSPOSE of \c Synthesize: \f$\partial/\partial V\f$ of \f$\sum_gw_gf(r_g)\f$, the channels' spline (values and second
    //! derivatives in \f$\log r\f$) and the exact powers past the table's ends read backwards.  Threaded over (site, point chunk) with a FIXED
    //! chunk count, each chunk's sums reduced in chunk order: the same bits at any thread count.
    channels_t SynthesizeAdjoint(const rvec_t& w) const
    {
        const size_t nA=NumSites(), nR=itsR.size(), n=itsPoints.size(), nCh=nR*itsNLM;
        if (w.size()!=n) throw std::logic_error("SiteMultipoles::SynthesizeAdjoint: the field is not on this mesh");
        const ColumnSplines logR(itsLogR);
        const size_t nChunk=std::max<size_t>(1,(theChunks+nA-1)/nA), chunk=(n+nChunk-1)/nChunk;
        std::vector<std::vector<double>> GV(nA*nChunk), GM(nA*nChunk);
        OuterFor(nA*nChunk, [&](size_t item)
        {
            const size_t a=item/nChunk, k=item%nChunk;
            std::vector<double>& gv=GV[item]; std::vector<double>& gm=GM[item];
            gv.assign(nCh,0.0); gm.assign(nCh,0.0);
            std::vector<double> y(itsNLM);
            for (size_t g=k*chunk; g<std::min(n,(k+1)*chunk); g++)
                if (w[g]!=0.0) AtPointAdjoint(logR,w[g],itsPoints[g]-itsCentres[a],y.data(),gv.data(),gm.data());
        });
        channels_t W(nA);
        OuterFor(nA, [&](size_t a)
        {
            std::vector<double> gv(nCh,0.0), gm(nCh,0.0);
            for (size_t k=0; k<nChunk; k++)
                for (size_t c=0; c<nCh; c++) {gv[c]+=GV[a*nChunk+k][c]; gm[c]+=GM[a*nChunk+k][c];}
            const std::vector<double> fromM=logR.SecondDerivativesAdjoint(gm,itsNLM);
            for (size_t c=0; c<nCh; c++) gv[c]+=fromM[c];
            W[a]=std::move(gv);
        });
        return W;
    }

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
    rvec3vec_t CrossSitePullback(const channels_t& V, const rvec_t& q) const
    {
        const size_t nA=NumSites(), n=itsPoints.size();
        if (q.size()!=n) throw std::logic_error("SiteMultipoles::CrossSitePullback: the charges are not on this mesh");
        const ColumnSplines logR(itsLogR);
        channels_t M(nA);
        OuterFor(nA, [&](size_t a) {M[a]=logR.SecondDerivatives(V[a],itsNLM);});
        const size_t nBlock=(n+theBlock-1)/theBlock;
        std::vector<std::vector<rvec3_t>> part(nBlock,std::vector<rvec3_t>(nA,rvec3_t(0,0,0)));
        OuterFor(nBlock, [&](size_t b)
        {
            std::vector<double> y(itsNLM), gu(3*itsNLM);
            std::vector<rvec3_t>& G=part[b];
            for (size_t g=b*theBlock; g<std::min(n,(b+1)*theBlock); g++)
            {
                if (q[g]==0.0) continue;   // an exactly uncharged point adds exactly nothing (no screen, pin 1)
                const size_t o=itsOwner[g];
                if (o>=nA) throw std::logic_error("SiteMultipoles::CrossSitePullback: a charged mesh point belongs to no site block");
                for (size_t a=0; a<nA; a++)
                {
                    if (a==o) continue;
                    const rvec3_t gr=q[g]*GradientAtPoint(logR,V[a],M[a],itsPoints[g]-itsCentres[a],y.data(),gu.data());
                    G[o]+=gr; G[a]-=gr;
                }
            }
        });
        rvec3vec_t F(nA);
        for (size_t a=0; a<nA; a++) F[a]=rvec3_t(0,0,0);
        for (const auto& G : part) for (size_t a=0; a<nA; a++) F[a]+=G[a];
        return F;
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
    static constexpr size_t theChunks=16;   //!< \c SynthesizeAdjoint's (site x chunk) items, at least: fixed, so the sum order is too
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

    //! \c AtPoint read backwards: \a w times each weight it gives a value (\a gv) or a second derivative (\a gm), at displacement \a d.
    void AtPointAdjoint(const ColumnSplines& logR, double w, const rvec3_t& d, double* y, double* gv, double* gm) const
    {
        const double r=norm(d);
        if (r==0.0) {gv[0]+=w*itsPmm[0]; return;}
        Ylm(d/r,y);
        const size_t nR=itsR.size();
        if (r<=itsR[0] || r>=itsR[nR-1])
        {
            const bool inner=r<=itsR[0];
            double* v=&gv[(inner ? 0 : nR-1)*itsNLM];
            for (int l=0; l<=itsLMax; l++)
            {
                const double scale=w*(inner ? std::pow(r/itsR[0],l) : std::pow(itsR[nR-1]/r,l+1));
                for (int c=l*l; c<(l+1)*(l+1); c++) v[c]+=scale*y[c];
            }
            return;
        }
        const double t=std::log(r);
        const size_t i=logR.Interval(t);
        double A,B,C,D; logR.Weights(i,t,A,B,C,D);
        A*=w; B*=w; C*=w; D*=w;
        double *v0=&gv[i*itsNLM], *v1=&gv[(i+1)*itsNLM], *m0=&gm[i*itsNLM], *m1=&gm[(i+1)*itsNLM];
        for (size_t c=0; c<itsNLM; c++) {const double yc=y[c]; v0[c]+=A*yc; v1[c]+=B*yc; m0[c]+=C*yc; m1[c]+=D*yc;}
    }

    //! \f$\nabla_d\sum_{lm}V_{lm}(|d|)Y_{lm}(\hat d)\f$, \c AtPoint's gradient (see \c CrossSitePullback).  At \f$d=0\f$ it answers zero: only a site's
    //! own zero-weight nucleus node sits there, and it carries no charge.
    rvec3_t GradientAtPoint(const ColumnSplines& logR, const std::vector<double>& V, const std::vector<double>& M, const rvec3_t& d,
                            double* y, double* gu) const
    {
        const double r=norm(d);
        if (r==0.0) return rvec3_t(0,0,0);
        const rvec3_t u=d/r;
        YlmGradient(u,y,gu);
        const size_t nR=itsR.size();
        double radial=0.0, tx=0.0, ty=0.0, tz=0.0;   // sum V' Y, and sum V grad_u Y
        auto add=[&](size_t c, double f, double fp)
        {
            radial+=fp*y[c];
            tx+=f*gu[3*c]; ty+=f*gu[3*c+1]; tz+=f*gu[3*c+2];
        };
        if (r<=itsR[0] || r>=itsR[nR-1])
        {
            const bool inner=r<=itsR[0];
            const double* v=&V[(inner ? 0 : nR-1)*itsNLM];
            for (int l=0; l<=itsLMax; l++)
            {
                const double scale=inner ? std::pow(r/itsR[0],l) : std::pow(itsR[nR-1]/r,l+1), power=inner ? double(l) : -double(l+1);
                for (int c=l*l; c<(l+1)*(l+1); c++) {const double f=scale*v[c]; add(c,f,power*f/r);}
            }
        }
        else
        {
            const double t=std::log(r);
            const size_t i=logR.Interval(t);
            double A,B,C,D;   logR.Weights(i,t,A,B,C,D);
            double Ad,Bd,Cd,Dd; logR.DerivativeWeights(i,t,Ad,Bd,Cd,Dd);
            const double *v0=&V[i*itsNLM], *v1=&V[(i+1)*itsNLM], *m0=&M[i*itsNLM], *m1=&M[(i+1)*itsNLM];
            for (size_t c=0; c<itsNLM; c++)
                add(c, A*v0[c]+B*v1[c]+C*m0[c]+D*m1[c], (Ad*v0[c]+Bd*v1[c]+Cd*m0[c]+Dd*m1[c])/r);   // d/dr = (1/r) d/d(log r)
        }
        const rvec3_t T(tx,ty,tz);
        return radial*u+(1.0/r)*(T-(u*T)*u);   // the tangential part of the polynomial's gradient: Y is read on the unit sphere
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

    //! \c Ylm and the gradient of the POLYNOMIAL it evaluates, \f$\bar P_l^m(z)/\sin^m\theta\f$ times the powers of \f$x+iy\f$, with x, y, z
    //! independent: \a gu holds \f$(\partial_x,\partial_y,\partial_z)\f$ at index \f$3(l^2+l+m)\f$.  Off the sphere the polynomial is not \f$Y_{lm}\f$, so
    //! a caller keeps only the tangential part (\c GradientAtPoint).
    void YlmGradient(const rvec3_t& u, double* y, double* gu) const
    {
        const int L=itsLMax+1;
        double c=1.0, s=0.0, cp=0.0, sp=0.0;   // (x+iy)^m and (x+iy)^(m-1)
        for (int m=0; m<L; m++)
        {
            double p0=itsPmm[m], p1=0.0, d0=0.0, d1=0.0;   // the z polynomial and its z derivative
            for (int l=m; l<L; l++)
            {
                if (l>m)
                {
                    const double A=itsA[l*L+m], B=itsB[l*L+m];
                    const double p=A*(u.z*p0-B*p1), d=A*(p0+u.z*d0-B*d1);
                    p1=p0; p0=p; d1=d0; d0=d;
                }
                if (m==0)
                {
                    const size_t k=l*l+l;
                    y[k]=p0; gu[3*k]=0.0; gu[3*k+1]=0.0; gu[3*k+2]=d0;
                }
                else
                {
                    const double r2=std::numbers::sqrt2;
                    const size_t kp=l*l+l+m, km=l*l+l-m;
                    y[kp]=r2*p0*c; y[km]=r2*p0*s;
                    // d(x+iy)^m/dx = m (x+iy)^(m-1), d/dy = i m (x+iy)^(m-1)
                    gu[3*kp]=r2*p0*m*cp;  gu[3*kp+1]=-r2*p0*m*sp; gu[3*kp+2]=r2*d0*c;
                    gu[3*km]=r2*p0*m*sp;  gu[3*km+1]= r2*p0*m*cp; gu[3*km+2]=r2*d0*s;
                }
            }
            cp=c; sp=s;
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
    std::vector<size_t> itsOwner;                //!< each point's site block (\c NumSites() for none)
    std::vector<char>   itsOnShell;              //!< 1 where a point is on a shell of its site (\c OnShell)
    std::vector<double> itsPmm, itsA, itsB;      //!< the \f$Y_{lm}\f$ recursion's constants
};

} //namespace
