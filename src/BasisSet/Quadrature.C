// File: BasisSet/Quadrature.C  The UNIVERSAL half of a fit basis's grid contract: its quadrature mesh.
//
// A FIT BASIS IS A FAMILY OF WEIGHT VECTORS OVER SHARED POINTS (user, 2026-08-22): projecting a field onto
// fit function c means integrating it against the weights W_c[g] = w_g c*(r_g), so a fit basis is not even
// DEFINED without a point set -- the mesh is CONSTITUTIVE of the basis, not a sibling of it.  This is the
// face that says so, and it is deliberately the NARROWEST one possible: points and weights, i.e. exactly
// qcMesh::Mesh, which every point set in the project already speaks.
//
// WHY IT IS ITS OWN FACE (doc/OpenWork.md "separation of concerns in the XC terms", step 1).  G_Quadrature
// used to BUNDLE two unrelated capabilities:
//
//     GridPoints(), Integral(f)                       -- any point set, any weights   (universal)
//     RhoOnGrid, ForwardFFT, GridCoeff, FieldCoeffs   -- RASTER ONLY
//
// so a fit basis over an atom-centred (Becke) mesh could provide the first row and not the second, and --
// the two being one interface -- could therefore answer NOTHING.  That is what made "delta fit on the fit
// raster" and "plane-wave fit on a Becke mesh" inexpressible, not any mathematics.  Split, the raster
// transforms keep their own face (G_RasterTransform, qchem.BasisSet.G_FieldEvaluator) and a consumer asks
// for exactly the half it consumes.
module;
export module qchem.BasisSet.Quadrature;
export import qchem.Mesh;   // qcMesh::Mesh -- the points+weights value type this face hands out

export namespace qchem::BasisSet
{

//! \brief The capability "I CARRY THE QUADRATURE I am defined over" -- points and weights, nothing else.
//!
//! Reached by an "I want more" cross-cast from a fit face a consumer already holds (never a cast into a
//! concrete basis), exactly as the \c G_* seams are.  A plane-wave fit basis answers with its FFT raster
//! (fractional corners \f$i/N\f$, weight \f$\Omega/N_{pts}\f$); a \f$\delta\f$ fit basis answers with the
//! mesh it was built over -- a uniform cell mesh, or an atom-centred Becke mesh with its site blocks.
//! Consumers integrate \f$\int f = \sum_g w_g f(r_g)\f$ (\c qcMesh::Integrate) and sample fields at
//! \c Mesh().Points(); neither needs to know which kind of mesh it got.
//!
//! \note The weights need NOT be quadrature weights in the narrow sense -- they are whatever the basis
//! contracts a field against.  That is what makes "a fit basis is a family of weight vectors" a TYPE and
//! not a metaphor.
class Quadrature
{
public:
    virtual ~Quadrature() = default;
    //! The points \f$r_g\f$ and weights \f$w_g\f$ this basis's fits and integrals run over.
    virtual const qcMesh::Mesh& Mesh() const=0;
};

} //namespace
