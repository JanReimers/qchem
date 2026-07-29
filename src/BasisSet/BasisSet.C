// File: BasisSet.C Quantum Chemistry basis set expressed as a sequence of Irrep basis sets.
module;
#include <vector>
#include <memory>

export module qchem.BasisSet;
export import qchem.BasisSet.Orbital_1E_IBS;
export import qchem.BasisSet.Fit_IBS;
export import qchem.Structure;
export import qchem.Symmetry;
export import qchem.ElectronConfiguration;
export import qchem.Matrix3D;   // Matrix3D -- the crystal point-group ops a periodic basis exposes (IBZ symmetrization)

import qchem.Iterators;
export import qchem.Streamable;

export namespace qchem::BasisSet
{
typedef std::vector<Irrep> irrepv_t; 

//----------------------------------------------------------------------------
//
//  Interface for a BasisSet which is assumed to a list of Irrep Basis Sets.
//
template <class T> class tBasisSet
    : public virtual Streamable
{
public:
    typedef Orbital_1E_IBS<T> obs_t;

    virtual ~tBasisSet() {};
    virtual size_t   GetNumFunctions() const=0;
    virtual irrepv_t GetIrreps(const Spin& ms) const=0;

    virtual FIT_CD_ABS<T>* CreateCDFitBasisSet(const Structure* cl, const qcMesh::MeshParams&) const;
    virtual FIT_SF_ABS<T>* CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams&) const;

    //! The crystal RECIPROCAL point group (τ=0 symmorphic ops) for IBZ density symmetrization.  Default {} =
    //! trivial {E} = no-op -- molecules / Γ / unfolded bases.  A periodic GPW basis returns the ops WHEN it
    //! folds the mesh (reduceBZ), so the composite density ctor-injects them and the reduced density is star-
    //! averaged (doc/GPWPlan1.md item 3).  It is a basis property (the basis computes the space group), so the
    //! density/WF read it here rather than being told via a setter.
    virtual std::vector<Matrix3D<double>> GetReciprocalPointOps() const {return {};}

    // Iterate() with no type argument yields the base obs_t* directly (no cast);
    // Iterate<D>() dynamic_cast's each IBS to the requested derived type D.
    // Built on the two primitives below, so storage stays private to the
    // concrete BasisSet.
    auto Iterate() const {return IndexProxy<const tBasisSet>(this, GetNumIBS());}
    template <class D> auto Iterate() const
    {
        return D_IndexProxy<const D, const tBasisSet>(this, GetNumIBS());
    }
    const obs_t* operator[](size_t i) const {return GetIBS(i);}

protected:
    virtual size_t      GetNumIBS()    const=0; //Number of Irrep basis sets.
    virtual const obs_t* GetIBS(size_t) const=0; //The only storage-specific primitive.
};

typedef tBasisSet<double>    Real_BS;
typedef tBasisSet<dcmplx> Complex_BS;

}//namespace