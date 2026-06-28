// File: Common/Math.C  The qchem.Math umbrella -- a convenience aggregator for outside consumers.
//
// Re-exports the STABLE, CHEAP numeric core: the <cmath> leaf (qchem.CMath, which itself brings the small
// Constants/IntPow/Factorials leaves) plus the fixed-size geometry types (Vector3D/Matrix3D).  Deliberately
// LEAN: it does NOT pull in the heavy / still-in-flux modules (Blaze, Types, FFT, FourierMap,
// ScalarFunction, VectorFunction) -- those stay explicit imports so a change to one of them doesn't force a
// rebuild of everything that just wanted sqrt, and so the (large) Blaze BMI is only loaded where it is used.
// As individual interfaces stabilise this umbrella can be widened with little rebuild-fan-out cost.
//
// NOTE: code INSIDE the math library must import qchem.CMath (the leaf), never this umbrella -- the umbrella
// re-exports Vector3D/Matrix3D, so importing it from those would be a cycle.
export module qchem.Math;
export import qchem.CMath;        // std::sqrt/… + Constants + IntPow + Factorials
