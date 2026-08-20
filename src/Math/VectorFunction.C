// File: VectorFunction.C  Mixin interface for real-space vector (basis) functions.
//
// When evaluated at a point r it returns a vector of values
// [phi_i(r)].  
module;
export module qchem.VectorFunction;
export import qchem.Types;

namespace qchem {

export template <class T> class VectorFunction
{
public:
    virtual ~VectorFunction()  {};

    virtual size_t       GetVectorSize()             const=0;
    virtual vec_t<T>     operator() (const rvec3_t&) const=0;
    virtual vec3vec_t<T> Gradient   (const rvec3_t&) const=0;

    //! \brief Batch evaluation: \f$[\phi_i(r_q)]\f$ at every point at once, as an
    //! (\a rs.size() \f$\times\f$ \c GetVectorSize()) table -- the MATRIX mirror of
    //! \c ScalarFunction::operator()(const rvec3vec_t&), which returns a vector because its pointwise
    //! form returns a scalar where this one returns a vector.  The default just loops the pointwise
    //! \c op(); an overrider that can beat \f$O(N_\text{pts})\f$ pointwise cost does so here.
    //! THE CASE THAT MOTIVATED IT: a periodic (Bloch) basis sums over lattice images per point, so a
    //! TRANSFORMED view of one (the spherical lattice view) applies its cart->sphere transform once per
    //! IMAGE where the transform is linear and once per POINT would do -- unfixable pointwise, because
    //! the caller holds an abstract basis and cannot know a transform is inside it (doc/OpenWork.md
    //! Step 3).  It is also the seam the \f$\Phi\f$-SPARSITY increment extends.
    virtual mat_t<T> operator()(const rvec3vec_t& rs) const
    {
        mat_t<T> P(rs.size(), GetVectorSize());
        for (size_t q=0;q<rs.size();q++)
        {
            const vec_t<T> v=(*this)(rs[q]);
            for (size_t i=0;i<v.size();i++) P(q,i)=v[i];
        }
        return P;
    }
};

} // namespace qchem