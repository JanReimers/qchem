// File: Fitting/Internal/DirectSumCoulomb.C  The plain Coulomb sum over a mesh's points.
//
// The simplest operator that must exist.  It is not the one to use: a Becke mesh is dense at nuclei and sparse between them,
// which is a poor rule for a singular kernel near each point (TeF₆: J good to about 1e-2 where radial Poisson reaches the
// mesh's own 5e-4, and no better at 80 radial points than at 50).  A density on it is still a correct, if coarse, answer.
//
// THE SELF TERM.  1/|r_g-r_g| is singular, so the diagonal is replaced by the integral of 1/r over the sphere of the point's
// own volume w_g, charge q_g spread uniformly through it:
//     (q_g/w_g) * int_0^R 4 pi r dr = (q_g/w_g) 2 pi R^2 = 3 q_g / (2 R),     R = (3 w_g / 4 pi)^(1/3).
// Evaluated at the sphere's centre, which is where the point is.
module;
#include <algorithm>
#include <cmath>
#include <numbers>
#include <vector>
export module qchem.Fitting.Internal.DirectSumCoulomb;
export import qchem.Fitting.CoulombOperator;
import qchem.Parallel;   // OuterFor: each block of points writes only its own slots
import qchem.Types;

export namespace qchem::Fitting
{

class DirectSumCoulomb : public virtual CoulombOperator
{
public:
    explicit DirectSumCoulomb(const qcMesh::Mesh& mesh)
        : itsX(mesh.size()), itsY(mesh.size()), itsZ(mesh.size()), itsW(mesh.Weights())
    {
        const rvec3vec_t& r=mesh.Points();
        for (size_t g=0; g<r.size(); g++) {itsX[g]=r[g].x; itsY[g]=r[g].y; itsZ[g]=r[g].z;}
    }

    //! \copydoc CoulombOperator::Potential
    //! A point with exactly zero charge adds exactly nothing, so it leaves the source list; that is no screen (pin 1).  Coincident
    //! points drop each other's pair: an MHL rule repeats its r = 0 node once per direction, all at zero weight, so nothing is lost.
    rvec_t Potential(const rvec_t& q) const override
    {
        const size_t n=itsW.size();
        std::vector<double> sx, sy, sz, sq;
        for (size_t h=0; h<n; h++)
            if (q[h]!=0.0) {sx.push_back(itsX[h]); sy.push_back(itsY[h]); sz.push_back(itsZ[h]); sq.push_back(q[h]);}
        rvec_t v(n);
        const size_t nBlock=(n+theBlock-1)/theBlock, ns=sq.size();
        OuterFor(nBlock, [&](size_t b)
        {
            for (size_t g=b*theBlock; g<std::min(n,(b+1)*theBlock); g++)
            {
                const double x=itsX[g], y=itsY[g], z=itsZ[g];
                double s=0.0;
                for (size_t h=0; h<ns; h++)
                {
                    const double dx=x-sx[h], dy=y-sy[h], dz=z-sz[h], r2=dx*dx+dy*dy+dz*dz;
                    s+= r2>0.0 ? sq[h]/std::sqrt(r2) : 0.0;
                }
                v[g]=s+SelfTerm(q[g],itsW[g]);
            }
        });
        return v;
    }

private:
    static constexpr size_t theBlock=256;   //!< points per OuterFor item: far more items than threads, so it balances
    static double SelfTerm(double q, double w)
    {
        if (w<=0.0) return 0.0;   // the Becke partition gives some points exactly zero weight (see OrthogonalFit)
        const double R=std::cbrt(3.0*w/(4.0*std::numbers::pi));
        return 1.5*q/R;
    }
    rvec_t itsX, itsY, itsZ, itsW;
};

} //namespace
