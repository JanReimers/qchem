// File: Mesh/Internal/Imp/RadialPoissonCoulomb.C  The bodies of qchem.Mesh.Internal.RadialPoissonCoulomb: the radial Poisson solve, its
// transpose, and the two operators on it (the interface carries the method and its measurements).
module;
#include <cmath>
#include <memory>
#include <numbers>
#include <vector>
module qchem.Mesh.Internal.RadialPoissonCoulomb;
import qchem.Parallel;   // OuterFor: every item writes only its own slots
import qchem.Blaze;

namespace qchem::qcMesh
{

//------------------------------------------------------------------------------------------------------------------------------------
//
//  RadialPoissonSolver
//
RadialPoissonSolver::RadialPoissonSolver(std::shared_ptr<const qcMesh::SiteMultipoles> sites)
    : itsSites(std::move(sites)), itsR(itsSites->Radii()), itsLMax(itsSites->LMax()), itsNLM(itsSites->NumLM())
    , itsLogR(MakeLog(itsR)), itsIndex(MakeIndexNodes(itsR.size())), itsInnerPower(std::log(itsR[1]/itsR[0])/std::log(2.0))
    , itsNodes(std::vector<double>(itsIndex.begin()+1,itsIndex.end()-1)), itsSmooth(MakeSmoothLog()), itsSmoothM(itsNodes.SecondDerivatives(itsSmooth,1))
{
    MakeGaussPoints();
}

RadialPoissonSolver::channels_t RadialPoissonSolver::Solve(const channels_t& Q) const
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

RadialPoissonSolver::channels_t RadialPoissonSolver::SolveAdjoint(const channels_t& W) const
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

std::vector<double> RadialPoissonSolver::MakeLog(const std::vector<double>& r)
{
    std::vector<double> s(r.size());
    for (size_t i=0; i<r.size(); i++) s[i]=std::log(r[i]);
    return s;
}

std::vector<double> RadialPoissonSolver::MakeIndexNodes(size_t nR)
{
    std::vector<double> t(nR+2);
    for (size_t i=0; i<nR+2; i++) t[i]=double(i);
    return t;
}

std::vector<double> RadialPoissonSolver::MakeSmoothLog() const
{
    std::vector<double> f(itsR.size());
    for (size_t i=0; i<f.size(); i++) f[i]=itsLogR[i]-itsInnerPower*std::log(double(i+1));
    return f;
}

void RadialPoissonSolver::RadialAt(double t, double& r) const
{
    if (t<=1.0) {r=itsR[0]*std::pow(t,itsInnerPower); return;}
    double A,B,C,D;
    const size_t i=itsNodes.Interval(t);
    itsNodes.Weights(i,t,A,B,C,D);
    r=std::exp(itsInnerPower*std::log(t)+A*itsSmooth[i]+B*itsSmooth[i+1]+C*itsSmoothM[i]+D*itsSmoothM[i+1]);
}

void RadialPoissonSolver::MakeGaussPoints()
{
    constexpr double gx[4]={-0.861136311594053,-0.339981043584856,0.339981043584856,0.861136311594053};
    constexpr double gw[4]={ 0.347854845137454, 0.652145154862546,0.652145154862546,0.347854845137454};
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

void RadialPoissonSolver::Solve(const std::vector<double>& Q, int l, std::vector<double>& V) const
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

void RadialPoissonSolver::SolveAdjoint(const std::vector<double>& W, int l, std::vector<double>& G) const
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

//------------------------------------------------------------------------------------------------------------------------------------
//
//  RadialPoissonCoulomb
//
rvec_t RadialPoissonCoulomb::Potential(const rvec_t& q) const
{
    const rvec_t Kq =itsSites->Synthesize(itsSolver.Solve(itsSites->Moments(q)));
    const rvec_t KTq=itsSites->MomentsAdjoint(itsSolver.SolveAdjoint(itsSites->SynthesizeAdjoint(q)));
    rvec_t v(Kq.size());
    // A point on no shell (the rule's zero-weight node at a nucleus) carries no charge, so K has no column for it: its potential is Kq.
    for (size_t g=0; g<v.size(); g++) v[g]=itsSites->OnShell(g) ? 0.5*(Kq[g]+KTq[g]) : Kq[g];
    return v;
}

rvec3vec_t RadialPoissonCoulomb::MakeSelfRepulsionPullback(const rvec_t& q) const
{
    return itsSites->CrossSitePullback(itsSolver.Solve(itsSites->Moments(q)), q);
}

//------------------------------------------------------------------------------------------------------------------------------------
//
//  MultipoleCoulomb
//
rvec_t MultipoleCoulomb::Potential(const rvec_t& q) const
{
    return itsSites->Flatten(itsSolver.Solve(itsSites->Split(q)));
}

double MultipoleCoulomb::SelfRepulsion(const rvec_t& q, const rvec_t& v) const
{
    rvec_t c(q.size());
    for (size_t a=0; a<c.size(); a++) c[a]=q[a]/itsMetric[a];
    return itsSites->OwnShellOverlap(itsSites->Split(c),itsSites->Split(v));
}

double MultipoleCoulomb::HartreeEnergy(double trDV, const rvec_t& q, const rvec_t& v) const
{
    return trDV-0.5*SelfRepulsion(q,v);
}

} //namespace
