// File: Hamiltonian/Internal/Imp/PP_NonLocal.C  Separable (Kleinman-Bylander) non-local pseudopotential term.
//
// Like PP_Local, the term owns only the MODEL (the KB projectors) and asks the STRUCTURE for its integration
// mesh (CreateIntegrationMesh -- the polymorphic geometry capability), then builds each projector's overlap
// vector b_p = <chi_i | beta_p(|r-R|) Y_lm> with the generic qcMesh::Overlap and accumulates the rank-1
// D|b><b|.  Geometry-neutral: no basis-specific dispatch, mesh type follows the geometry.
module;
#include <cassert>
#include <cmath>
#include <iostream>

module qchem.Hamiltonian.Internal.Terms;
import qchem.Energy;
import qchem.Mesh.Quadrature;           // qcMesh::Overlap(mesh, VectorFunction, ScalarFunction) -> projection vector
import qchem.VectorFunction;            // the orbital basis IS-A VectorFunction (passed direct)
import qchem.Blaze;                     // rmat_t, rvec_t, trans, outer product
import qchem.Math;                      // norm, Pi, sqrt

namespace qchem::Hamiltonian
{

namespace {

// Real (tesseral) spherical harmonic Y_lm(rhat), UNIT-NORMALISED on the sphere (integral |Y|^2 dOmega = 1),
// evaluated from the UNIT vector n=(x,y,z).  Only Sum_m |Y_lm><Y_lm| (the degree-l projector) enters the KB
// term, so any orthonormal real degree-l set is equivalent; these are the standard Cartesian forms.  m runs
// -l..+l.  GTH/HGH projectors reach l=3 (f) at most.
double RealYlm(int l, int m, double x, double y, double z)
{
    using std::sqrt;
    const double pi=Pi;
    switch (l)
    {
    case 0: return 0.5/sqrt(pi);                                            // Y_00
    case 1:
        switch (m)
        {
        case -1: return sqrt(0.75/pi)*y;
        case  0: return sqrt(0.75/pi)*z;
        case  1: return sqrt(0.75/pi)*x;
        }
        break;
    case 2:
        switch (m)
        {
        case -2: return sqrt(15.0/(4*pi))*x*y;
        case -1: return sqrt(15.0/(4*pi))*y*z;
        case  0: return sqrt( 5.0/(16*pi))*(3*z*z-1.0);                     // x^2+y^2+z^2=1
        case  1: return sqrt(15.0/(4*pi))*x*z;
        case  2: return sqrt(15.0/(16*pi))*(x*x-y*y);
        }
        break;
    case 3:
        switch (m)
        {
        case -3: return sqrt( 35.0/(32*pi))*y*(3*x*x-y*y);
        case -2: return sqrt(105.0/(4 *pi))*x*y*z;
        case -1: return sqrt( 21.0/(32*pi))*y*(5*z*z-1.0);
        case  0: return sqrt(  7.0/(16*pi))*z*(5*z*z-3.0);
        case  1: return sqrt( 21.0/(32*pi))*x*(5*z*z-1.0);
        case  2: return sqrt(105.0/(16*pi))*z*(x*x-y*y);
        case  3: return sqrt( 35.0/(32*pi))*x*(x*x-3*y*y);
        }
        break;
    }
    assert(false && "RealYlm: unsupported (l,m)");
    return 0.0;
}

// One KB projector channel as a scalar field: beta_p(|r-R|) * Y_lm((r-R)^).  At r=R the angular factor is
// singular but beta_p(0)=0 for l>=1 (the r^l prefactor) and finite for l=0, so the product is well-defined;
// guard the unit-vector division.
class BetaYlmField : public ScalarFunction<double>
{
    rvec3_t R; const Pseudopotential::SeparablePotential_R& v; int Z; size_t p; int l, m;
public:
    BetaYlmField(const rvec3_t& R_, const Pseudopotential::SeparablePotential_R& v_, int Z_, size_t p_, int l_, int m_)
        : R(R_), v(v_), Z(Z_), p(p_), l(l_), m(m_) {}
    double operator()(const rvec3_t& r) const override
    {
        rvec3_t d=r-R; double rr=norm(d);
        if (rr<1e-12) return l==0 ? v.BetaR(Z,p,0.0)*RealYlm(0,0,0,0,0) : 0.0;
        return v.BetaR(Z,p,rr) * RealYlm(l,m, d.x/rr, d.y/rr, d.z/rr);
    }
    rvec3_t Gradient(const rvec3_t&) const override {return rvec3_t(0,0,0);}   // unused by Overlap
};

} //anon

PP_NonLocal::PP_NonLocal(const st_t& st, sep_t sep, const qcMesh::MeshParams& mp)
    : theStructure(st), itsSep(std::move(sep)), itsMeshParams(mp)
{
    assert(theStructure);
    assert(itsSep);
}

rsmat_t PP_NonLocal::CalculateMatrix(const robs_t* bs, const Spin&) const
{
    qcMesh::Mesh mesh = theStructure->CreateIntegrationMesh(itsMeshParams);   // the geometry's own mesh
    size_t n=bs->GetVectorSize();
    rmat_t V(n,n,0.0);
    for (size_t a=0; a<theStructure->GetNumAtoms(); a++)
    {
        const Atom* at=(*theStructure)[a];
        int Z=at->itsZ;
        for (size_t p=0; p<itsSep->NumProjectors(Z); p++)
        {
            int    l=itsSep->AngularMomentum(Z,p);
            double D=itsSep->Coefficient    (Z,p);
            for (int m=-l; m<=l; m++)
            {
                rvec_t b=qcMesh::Overlap(mesh, *bs, BetaYlmField(at->itsR,*itsSep,Z,p,l,m));
                for (size_t i=0;i<n;i++)         // rank-1 update D|b><b| (real symmetric, bit-exact)
                    for (size_t j=0;j<n;j++) V(i,j)+=D*b[i]*b[j];
            }
        }
    }
    return rsmat_t(V);
}

void PP_NonLocal::GetEnergy(EnergyBreakdown& te, const rDM_CD* cd) const
{
    te.Een += cd->DM_Contract(this);   // electron-ion (KB nonlocal) energy = Tr(D V_NL)
}

std::ostream& PP_NonLocal::Write(std::ostream& os) const
{
    return os << "    Separable (KB) non-local pseudopotential (mesh quadrature), " << theStructure->GetNumAtoms()
              << " atom(s)." << std::endl;
}

} //namespace
