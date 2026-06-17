// File: BasisSet/Atom/Evaluators/Internal/RelAngularIntegrals.C
// Angular coefficients Ak for <LL|LL> relativistic ERIs.
// The large component of spinor (κ,mj) is Σ_{ms} CG(la,mj-ms;½,ms|ja,mj) Y_la^{mj-ms} χ_ms.
// Spin integration reduces each integral to a sum over ms of squared CG weights
// times the corresponding NR angular integral.
module;
export module qchem.BasisSet.Atom.Evaluators.Internal.RelAngularIntegrals;
export import qchem.Types;

export namespace RelAngularIntegrals
{
    // Ak coefficients for direct  <κa mja, κc mjc | 1/r12 | κa mja, κc mjc>
    rvec11_t Direct  (int κa, int κc, double mja, double mjc);
    // Ak coefficients for exchange <κa mja, κb mjb | 1/r12 | κb mjb, κa mja>
    rvec11_t Exchange(int κa, int κb, double mja, double mjb);
    // mj-summed versions (averaged over all mj, useful for unpolarized calculations)
    rvec11_t Direct  (int κa, int κc);
    rvec11_t Exchange(int κa, int κb);
}
