// File: Internal/AngularTypes.C  Concrete angular mesh classes (transplanted numerics).
module;
export module qchem.Mesh1.Angular.Internal;
export import qchem.Mesh1.Angular;

//! Gauss product-free angular schemes for numDir in {1,2,6,8,12,24,30,32,50} (exact to L=0..11).
export class GaussAngularMesh : public AngularMesh
{
public:
    explicit GaussAngularMesh(int numDir);
};

//! Gauss-Legendre in theta (x cos-nodes) times uniform phi.  Exact for spherical harmonics up to L.
export class GaussLegendreAngularMesh : public AngularMesh
{
public:
    explicit GaussLegendreAngularMesh(int L);
};

//! Euler-Maclaren in theta times uniform phi.  m in {1,2,3} controls the theta clustering.
export class EulerMaclarenAngularMesh : public AngularMesh
{
public:
    EulerMaclarenAngularMesh(int L, int m);
};
