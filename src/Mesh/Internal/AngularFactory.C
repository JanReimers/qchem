// File: Internal/AngularFactory.C  MakeAngular -- dispatch to the per-scheme builder functions.
module;
#include <cmath>
#include <stdexcept>
module qchem.Mesh.Angular;

namespace qchem::qcMesh
{

//! p.angRot: rotate the whole direction set rigidly (Rodrigues, about the fixed generic axis
//! (1,2,3)/sqrt(14)).  Exactness degree is rotation-invariant; the point is steering a rule's
//! special orbits OFF structure axes (the Lebedev-on-bond-axes fix, free runs only -- an IMPOSED
//! run's site-adapted grid must stay G_s-invariant and never routes through this factory's grid).
static AngularMesh Rotate(AngularMesh m, double angle)
{
    const double n = 1.0/std::sqrt(14.0);
    const rvec3_t k(1.0*n, 2.0*n, 3.0*n);
    const double c = std::cos(angle), s = std::sin(angle);
    rvec3vec_t d(m.size());
    for (size_t i = 0; i < m.size(); ++i)
    {
        const rvec3_t& v = m.Dirs()[i];
        d[i] = v*c + Cross(k,v)*s + k*((k*v)*(1.0-c));
    }
    return AngularMesh(std::move(d), rvec_t(m.W()));
}

AngularMesh MakeAngular(const MeshParams& p)
{
    AngularMesh m;
    switch (p.angular)
    {
    case AngularKind::Lebedev:       m = LebedevAngular(p.nAngular); break;
    case AngularKind::GaussLegendre: m = GaussLegendreAngular(p.nAngular); break;
    case AngularKind::EulerMaclaren: m = EulerMaclarenAngular(p.nAngular, p.em_m); break;
    default: throw std::runtime_error("MakeAngular: unknown AngularKind");
    }
    return (p.angRot != 0.0) ? Rotate(std::move(m), p.angRot) : m;
}

} //namespace qchem::qcMesh
