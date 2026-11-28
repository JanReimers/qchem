// File: Mesh/Imp/SiteMultipoles.C  The bodies of qchem.Mesh.SiteMultipoles: column splines, and the transforms between point values and
// radial channels (the interface says what each one is and who asks).
module;
#include <algorithm>
#include <cmath>
#include <numbers>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.Mesh.SiteMultipoles;
import qchem.Parallel;   // OuterFor: every item writes only its own slots

namespace qchem::qcMesh
{

//------------------------------------------------------------------------------------------------------------------------------------
//
//  ColumnSplines
//
std::vector<double> ColumnSplines::SecondDerivatives(const std::vector<double>& y, size_t nCol) const
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

std::vector<double> ColumnSplines::SecondDerivativesAdjoint(const std::vector<double>& g, size_t nCol) const
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

size_t ColumnSplines::Interval(double t) const
{
    const size_t i=std::upper_bound(itsX.begin(),itsX.end(),t)-itsX.begin();
    return std::min(std::max<size_t>(i,1),itsX.size()-1)-1;
}

void ColumnSplines::Weights(size_t i, double t, double& A, double& B, double& C, double& D) const
{
    const double h=itsX[i+1]-itsX[i];
    B=(t-itsX[i])/h; A=1.0-B;
    C=(A*A*A-A)*h*h/6.0; D=(B*B*B-B)*h*h/6.0;
}

void ColumnSplines::DerivativeWeights(size_t i, double t, double& A, double& B, double& C, double& D) const
{
    const double h=itsX[i+1]-itsX[i], b=(t-itsX[i])/h, a=1.0-b;
    A=-1.0/h; B=1.0/h;
    C=-(3.0*a*a-1.0)*h/6.0; D=(3.0*b*b-1.0)*h/6.0;
}

std::vector<double> ColumnSplines::Solve(const std::vector<double>& rhs, size_t nCol) const
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

//------------------------------------------------------------------------------------------------------------------------------------
//
//  SiteMultipoles
//
SiteMultipoles::SiteMultipoles(const Mesh& mesh, const rvec3vec_t& centres, const RadialMesh& radial, const AngularMesh& angular, int lMax)
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
    itsItem.assign(itsPoints.size(),theNoShell);
    for (size_t item=0; item<itsShell.size(); item++) for (size_t g : itsShell[item]) itsItem[g]=item;
    itsOwner.assign(itsPoints.size(),centres.size());   // NumSites() = no owner
    for (size_t a=0; a<centres.size(); a++)
        for (size_t g=mesh.SiteBegin(a); g<mesh.SiteEnd(a); g++) itsOwner[g]=a;
}

std::vector<size_t> SiteMultipoles::All() const
{
    std::vector<size_t> idx(itsPoints.size());
    for (size_t g=0; g<idx.size(); g++) idx[g]=g;
    return idx;
}

SiteMultipoles::channels_t SiteMultipoles::Moments          (const rvec_t& q) const {return MomentsAt(All(),q);}
rvec_t                     SiteMultipoles::Synthesize       (const channels_t& V) const {return SynthesizeAt(V,All());}
rvec_t                     SiteMultipoles::MomentsAdjoint   (const channels_t& G) const {return MomentsAdjointAt(G,All());}
SiteMultipoles::channels_t SiteMultipoles::SynthesizeAdjoint(const rvec_t& w) const
{
    if (w.size()!=itsPoints.size()) throw std::logic_error("SiteMultipoles::SynthesizeAdjoint: the field is not on this mesh");
    return SynthesizeAdjointAt(All(),w);
}

SiteMultipoles::channels_t SiteMultipoles::MomentsAt(const std::vector<size_t>& idx, const rvec_t& q) const
{
    if (q.size()!=idx.size()) throw std::logic_error("SiteMultipoles::MomentsAt: one charge per listed point");
    const size_t nA=NumSites(), nR=itsR.size();
    // The listed points per (site, shell) item, in idx order: on the full list that is each shell's own (ascending) order.
    std::vector<std::vector<size_t>> onItem(nA*nR);
    for (size_t k=0; k<idx.size(); k++)
        if (itsItem[idx[k]]!=theNoShell) onItem[itsItem[idx[k]]].push_back(k);
    channels_t Q(nA,std::vector<double>(nR*itsNLM,0.0));
    OuterFor(nA*nR, [&](size_t item)
    {
        const size_t a=item/nR, i=item%nR;
        std::vector<double> y(itsNLM);
        double* out=&Q[a][i*itsNLM];
        for (size_t k : onItem[item])
        {
            const rvec3_t d=itsPoints[idx[k]]-itsCentres[a];
            Ylm(d/norm(d),y.data());
            for (size_t c=0; c<itsNLM; c++) out[c]+=q[k]*y[c];
        }
    });
    return Q;
}

rvec_t SiteMultipoles::SynthesizeAt(const channels_t& V, const std::vector<size_t>& idx) const
{
    const size_t nA=NumSites();
    const ColumnSplines logR(itsLogR);
    channels_t M(nA);
    OuterFor(nA, [&](size_t a) {M[a]=logR.SecondDerivatives(V[a],itsNLM);});
    const size_t n=idx.size(), nBlock=(n+theBlock-1)/theBlock;
    rvec_t v(n);
    OuterFor(nBlock, [&](size_t b)
    {
        std::vector<double> y(itsNLM);
        for (size_t k=b*theBlock; k<std::min(n,(b+1)*theBlock); k++)
        {
            double s=0.0;
            for (size_t a=0; a<nA; a++) s+=AtPoint(logR,V[a],M[a],itsPoints[idx[k]]-itsCentres[a],y.data());
            v[k]=s;
        }
    });
    return v;
}

rvec_t SiteMultipoles::MomentsAdjointAt(const channels_t& G, const std::vector<size_t>& idx) const
{
    const size_t nR=itsR.size(), n=idx.size(), nBlock=(n+theBlock-1)/theBlock;
    rvec_t w(n);
    OuterFor(nBlock, [&](size_t b)   // each listed point reads only its own shell's channels
    {
        std::vector<double> y(itsNLM);
        for (size_t k=b*theBlock; k<std::min(n,(b+1)*theBlock); k++)
        {
            const size_t g=idx[k], item=itsItem[g];
            if (item==theNoShell) {w[k]=0.0; continue;}
            const size_t a=item/nR, i=item%nR;
            const rvec3_t d=itsPoints[g]-itsCentres[a];
            Ylm(d/norm(d),y.data());
            const double* in=&G[a][i*itsNLM];
            double s=0.0;
            for (size_t c=0; c<itsNLM; c++) s+=in[c]*y[c];
            w[k]=s;
        }
    });
    return w;
}

SiteMultipoles::channels_t SiteMultipoles::SynthesizeAdjointAt(const std::vector<size_t>& idx, const rvec_t& w) const
{
    const size_t nA=NumSites(), nR=itsR.size(), n=idx.size(), nCh=nR*itsNLM;
    if (w.size()!=n) throw std::logic_error("SiteMultipoles::SynthesizeAdjointAt: one weight per listed point");
    const ColumnSplines logR(itsLogR);
    const size_t nChunk=std::max<size_t>(1,(theChunks+nA-1)/nA), chunk=(n+nChunk-1)/nChunk;
    std::vector<std::vector<double>> GV(nA*nChunk), GM(nA*nChunk);
    OuterFor(nA*nChunk, [&](size_t item)
    {
        const size_t a=item/nChunk, c=item%nChunk;
        std::vector<double>& gv=GV[item]; std::vector<double>& gm=GM[item];
        gv.assign(nCh,0.0); gm.assign(nCh,0.0);
        std::vector<double> y(itsNLM);
        for (size_t k=c*chunk; k<std::min(n,(c+1)*chunk); k++)
            if (w[k]!=0.0) AtPointAdjoint(logR,w[k],itsPoints[idx[k]]-itsCentres[a],y.data(),gv.data(),gm.data());
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

namespace
{
//! Solve the symmetric positive definite \a A (n x n, row-major) X = B for the n x m \a B in place, by Cholesky.  THROWS on a non-positive pivot.
void CholeskySolve(std::vector<double> A, std::vector<double>& B, size_t n, size_t m)
{
    for (size_t j=0; j<n; j++)
    {
        double d=A[j*n+j];
        for (size_t k=0; k<j; k++) d-=A[j*n+k]*A[j*n+k];
        if (!(d>0.0)) throw std::logic_error("SiteMultipoles::MakeChannelAction: the directions do not fix D^l (a singular normal matrix)");
        A[j*n+j]=std::sqrt(d);
        for (size_t i=j+1; i<n; i++)
        {
            double s=A[i*n+j];
            for (size_t k=0; k<j; k++) s-=A[i*n+k]*A[j*n+k];
            A[i*n+j]=s/A[j*n+j];
        }
    }
    for (size_t c=0; c<m; c++)
    {
        for (size_t i=0; i<n; i++)
        {
            double s=B[i*m+c];
            for (size_t k=0; k<i; k++) s-=A[i*n+k]*B[k*m+c];
            B[i*m+c]=s/A[i*n+i];
        }
        for (size_t i=n; i-- >0;)
        {
            double s=B[i*m+c];
            for (size_t k=i+1; k<n; k++) s-=A[k*n+i]*B[k*m+c];
            B[i*m+c]=s/A[i*n+i];
        }
    }
}
} //anonymous namespace

SiteMultipoles::ChannelAction SiteMultipoles::MakeChannelAction(const std::vector<Matrix3D<double>>& W,
                                                                const std::vector<std::vector<size_t>>& image) const
{
    if (W.size()!=image.size() || W.empty())
        throw std::invalid_argument("SiteMultipoles::MakeChannelAction: one site permutation per op, and at least one op");
    for (const auto& perm : image)
    {
        if (perm.size()!=NumSites()) throw std::invalid_argument("SiteMultipoles::MakeChannelAction: a site permutation of the wrong length");
        for (size_t b : perm) if (b>=NumSites()) throw std::invalid_argument("SiteMultipoles::MakeChannelAction: a site index out of range");
    }
    // Generic directions (a Fibonacci sphere), NOT the angular rule's: a low-degree rule is too small, and its own symmetry makes the high-l
    // harmonics linearly dependent on it (H₂O at the test rung: a singular normal matrix).  2(lMax+1)^2 of them over-determine every block.
    const size_t nD=2*itsNLM;
    rvec3vec_t dirs(nD);
    std::vector<double> yU(nD*itsNLM), yW(nD*itsNLM);
    for (size_t d=0; d<nD; d++)
    {
        const double z=1.0-(2.0*d+1.0)/double(nD), rho=std::sqrt(1.0-z*z), phi=d*std::numbers::pi*(3.0-std::sqrt(5.0));
        dirs[d]=rvec3_t(rho*std::cos(phi),rho*std::sin(phi),z);
        Ylm(dirs[d],&yU[d*itsNLM]);
    }
    ChannelAction act{image,std::vector<std::vector<std::vector<double>>>(W.size(),std::vector<std::vector<double>>(itsLMax+1))};
    for (size_t k=0; k<W.size(); k++)
    {
        for (size_t d=0; d<nD; d++) Ylm(W[k]*dirs[d],&yW[d*itsNLM]);
        for (int l=0; l<=itsLMax; l++)
        {
            const size_t n=2*l+1, o=size_t(l)*l;
            // Rows d: A_d = Y_l(u_d), B_d = Y_l(W u_d) = A_d D^T.  Normal equations (A^T A) D^T = A^T B.
            std::vector<double> AtA(n*n,0.0), AtB(n*n,0.0);
            for (size_t d=0; d<nD; d++)
            {
                const double* a=&yU[d*itsNLM+o];
                const double* b=&yW[d*itsNLM+o];
                for (size_t i=0; i<n; i++)
                    for (size_t j=0; j<n; j++) {AtA[i*n+j]+=a[i]*a[j]; AtB[i*n+j]+=a[i]*b[j];}
            }
            CholeskySolve(AtA,AtB,n,n);   // AtB now holds D^T
            std::vector<double>& D=act.D[k][l];
            D.assign(n*n,0.0);
            for (size_t i=0; i<n; i++) for (size_t j=0; j<n; j++) D[i*n+j]=AtB[j*n+i];
            double res=0.0;
            for (size_t d=0; d<nD; d++)
            {
                const double* a=&yU[d*itsNLM+o];
                const double* b=&yW[d*itsNLM+o];
                for (size_t i=0; i<n; i++)
                {
                    double s=0.0;
                    for (size_t j=0; j<n; j++) s+=D[i*n+j]*a[j];
                    res=std::max(res,std::fabs(s-b[i]));
                }
            }
            if (!(res<1e-9))
                throw std::logic_error("SiteMultipoles::MakeChannelAction: Y_l(W u) is not D^l Y_l(u) (residual "+std::to_string(res)
                                       +" at l="+std::to_string(l)+"): W is not a rotation, or the rule's directions do not fix D^l");
        }
    }
    return act;
}

void SiteMultipoles::SymmetrizeChannels(channels_t& G, const ChannelAction& act) const
{
    const size_t nA=NumSites(), nR=itsR.size();
    if (G.size()!=nA || act.image.empty() || act.image.front().size()!=nA)
        throw std::invalid_argument("SiteMultipoles::SymmetrizeChannels: the channels or the action do not match this mesh's sites");
    const double inv=1.0/double(act.image.size());
    channels_t out(nA,std::vector<double>(nR*itsNLM,0.0));
    OuterFor(nA, [&](size_t a)   // each site writes only its own channels
    {
        for (size_t k=0; k<act.image.size(); k++)
        {
            const std::vector<double>& src=G[act.image[k][a]];
            for (int l=0; l<=itsLMax; l++)
            {
                const size_t n=2*l+1, o=size_t(l)*l;
                const std::vector<double>& D=act.D[k][l];
                for (size_t i=0; i<nR; i++)
                {
                    const double* in=&src[i*itsNLM+o];
                    double* to=&out[a][i*itsNLM+o];
                    for (size_t m=0; m<n; m++)
                        for (size_t mp=0; mp<n; mp++) to[mp]+=D[m*n+mp]*in[m];   // D^T in
                }
            }
        }
        for (double& x : out[a]) x*=inv;
    });
    G=std::move(out);
}

rvec3vec_t SiteMultipoles::CrossSitePullback(const channels_t& V, const rvec_t& q) const
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

double SiteMultipoles::OwnShellOverlap(const channels_t& F, const channels_t& G) const
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

rvec_t SiteMultipoles::Flatten(const channels_t& Q) const
{
    rvec_t f(NumChannels());
    size_t k=0;
    for (const auto& site : Q) for (double x : site) f[k++]=x;
    return f;
}

SiteMultipoles::channels_t SiteMultipoles::Split(const rvec_t& f) const
{
    channels_t Q(NumSites(),std::vector<double>(itsR.size()*itsNLM));
    size_t k=0;
    for (auto& site : Q) for (double& x : site) x=f[k++];
    return Q;
}

std::vector<double> SiteMultipoles::Positive(const rvec_t& r, const rvec_t& x)
{
    std::vector<double> s;
    for (size_t i=0; i<r.size(); i++) if (r[i]>0.0) s.push_back(x[i]);
    return s;
}

std::vector<double> SiteMultipoles::MakeLog(const std::vector<double>& r)
{
    std::vector<double> s(r.size());
    for (size_t i=0; i<r.size(); i++) s[i]=std::log(r[i]);
    return s;
}

double SiteMultipoles::AtPoint(const ColumnSplines& logR, const std::vector<double>& V, const std::vector<double>& M, const rvec3_t& d, double* y) const
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

void SiteMultipoles::AtPointAdjoint(const ColumnSplines& logR, double w, const rvec3_t& d, double* y, double* gv, double* gm) const
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

rvec3_t SiteMultipoles::GradientAtPoint(const ColumnSplines& logR, const std::vector<double>& V, const std::vector<double>& M, const rvec3_t& d,
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

void SiteMultipoles::MakeYlmTables()
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

void SiteMultipoles::Ylm(const rvec3_t& u, double* y) const
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

void SiteMultipoles::YlmGradient(const rvec3_t& u, double* y, double* gu) const
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

} //namespace
