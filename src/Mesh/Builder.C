// File: Builder.C  Incremental Mesh accumulation, the efficient way.
//
// When a builder appends points one at a time and the final count is not known up front (the Becke
// molecular mesh: each atom contributes a variable number of surviving points), Mesh::Append's
// resize-by-one is O(N^2).  MeshBuilder accumulates into the project's doubling blazem::VecBuilder
// (O(N) total) and hands over a finished Mesh via the from-arrays constructor.  Opt-in module so
// the core qchem.Mesh stays free of the qchem.Blaze dependency.
module;
#include <stdexcept>
#include <string>
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
    //! The atom-centred form: \a u is the point's offset from its owner's centre (\c Mesh::OwnerOffsets).  A builder appends every point
    //! one way or the other; mixing the two is a defect \c take refuses.
    void   Append(const rvec3_t& r, double w, const rvec3_t& u) {itsR.Append(r); itsW.Append(w); itsU.Append(u);}
    size_t size  () const                     {return itsW.size();}
    //! \brief Open a new SITE BLOCK: every point appended from here on belongs to the next centre.
    //! An atom-centred builder calls this once per centre BEFORE appending that centre's surviving
    //! points, and the finished mesh can then answer \c SiteIntegrals -- \f$\int w_A f\f$ per site, which
    //! is what an ATOMIC moment or charge actually is (doc/OpenWork.md Step 0a).  A builder that never
    //! calls it produces a mesh with no site structure, exactly as before.
    void   BeginSite()                        {itsSiteStart.push_back(itsW.size());}
    Mesh   take  ()
    {
        const size_t n=itsW.size(), nu=itsU.size();
        if (nu!=0 && nu!=n) throw std::logic_error("MeshBuilder::take: owner offsets for "+std::to_string(nu)+" of "+std::to_string(n)+" points");
        return Mesh(itsR.take(), itsW.take(), std::move(itsSiteStart), nu ? itsU.take() : rvec3vec_t());
    }
private:
    blazem::VecBuilder<rvec3_t> itsU;   //!< the owner offsets (empty: not recorded)
    blazem::VecBuilder<rvec3_t> itsR;
    blazem::VecBuilder<double>  itsW;
    std::vector<size_t>         itsSiteStart;
};

} //export namespace qchem::qcMesh
