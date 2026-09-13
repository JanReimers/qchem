// File: Math/Gaussian.C  One radial Gaussian term -- a function, with no context.
module;
export module qchem.Math.Gaussian;

export namespace qchem::Math
{

//! \brief One term \f$c\,r^{2n}\,e^{-\alpha r^2}\f$: an even polynomial times a Gaussian, in a radial variable.
//! Pure data, no dependencies.  Anything that expands a radial function in closed Gaussian form -- a
//! pseudopotential's short-range core, a Kleinman-Bylander projector, a fitted density -- speaks in these; a
//! consumer holding a Gaussian basis can then integrate against them analytically instead of by quadrature.
//! (V1.2: this ONE type replaced two byte-identical structs that lived in two modules purely so neither would
//! name the other.  A channel's \f$r^l\f$, where there is one, is carried by the owner, not the term.)
struct Gaussian
{
    double c;       //!< coefficient
    int    n;       //!< polynomial degree, as \f$r^{2n}\f$
    double alpha;   //!< exponent
};

} //namespace
