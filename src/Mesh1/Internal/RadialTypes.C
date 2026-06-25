// File: Internal/RadialTypes.C  Concrete radial mesh classes (transplanted numerics).
module;
export module qchem.Mesh1.Radial.Internal;
export import qchem.Mesh1.Radial;

//! Murray-Handy-Laming radial mesh (MOLECULAR PHYSICS 1993, 78, 997).
export class MHLRadialMesh : public RadialMesh
{
public:
    MHLRadialMesh(int NumPoints, int m, double alpha);
};

//! Logarithmic radial mesh: r_i = start * q^i.
export class LogRadialMesh : public RadialMesh
{
public:
    LogRadialMesh(double start, double stop, int NumPoints);
};

//! Uniform radial mesh: r_i evenly spaced, trapezoidal r^2 weights.  Elementary quadrature
//! (no transplanted kernel) -- present for completeness; MHL/Log are the workhorses.
export class LinearRadialMesh : public RadialMesh
{
public:
    LinearRadialMesh(double start, double stop, int NumPoints);
};
