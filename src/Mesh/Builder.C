// File: Builder.C  Incremental Mesh accumulation, the efficient way.
//
// When a builder appends points one at a time and the final count is not known up front (the Becke
// molecular mesh: each atom contributes a variable number of surviving points), Mesh::Append's
// resize-by-one is O(N^2).  MeshBuilder accumulates into the project's doubling blazem::VecBuilder
// (O(N) total) and hands over a finished Mesh via the from-arrays constructor.  Opt-in module so
// the core qchem.Mesh stays free of the qchem.Blaze dependency.
module;
#include <vector>            // the site-block starts handed to the finished Mesh
export module qchem.Mesh.Builder;
export import qchem.Mesh;
import qchem.Blaze;          // blazem::VecBuilder

export namespace qchem::qcMesh
{

class MeshBuilder
{
public:
    void   Append(const rvec3_t& r, double w) {itsR.Append(r); itsW.Append(w);}
    size_t size  () const                     {return itsW.size();}
    //! \brief Open a new SITE BLOCK: every point appended from here on belongs to the next centre.
    //! An atom-centred builder calls this once per centre BEFORE appending that centre's surviving
    //! points, and the finished mesh can then answer \c SiteIntegrals -- \f$\int w_A f\f$ per site, which
    //! is what an ATOMIC moment or charge actually is (doc/OpenWork.md Step 0a).  A builder that never
    //! calls it produces a mesh with no site structure, exactly as before.
    void   BeginSite()                        {itsSiteStart.push_back(itsW.size());}
    Mesh   take  ()                           {return Mesh(itsR.take(), itsW.take(), std::move(itsSiteStart));}
private:
    blazem::VecBuilder<rvec3_t> itsR;
    blazem::VecBuilder<double>  itsW;
    std::vector<size_t>         itsSiteStart;
};

} //export namespace qchem::qcMesh
