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
export import qchem.Symmetry.Lattice_3D.SpaceGroup;   // ReciprocalOp {U|τ} -- the density-symmetrization ops (glide phase)

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

    //! The crystal RECIPROCAL point group as \f${U|\tau}\f$ ops for IBZ density symmetrization.  Default {} =
    //! trivial {E} = no-op -- molecules / Γ / unfolded bases.  A periodic GPW basis returns the ops WHEN it
    //! folds the mesh (reduceBZ), so the composite density ctor-injects them and the reduced density is star-
    //! averaged with the glide phase \f$e^{+2\pi i(Um)\cdot\tau}\f$ (doc/GPWPlan1.md items 3 + 5).  It is a basis
    //! property (the basis computes the space group), so the density/WF read it here rather than via a setter.
    virtual std::vector<Symmetry::Lattice_3D::ReciprocalOp> GetReciprocalPointOps() const {return {};}

    //! The FULL DETECTED crystal reciprocal point group \f${U|\tau}\f$ -- populated on every periodic run,
    //! IMPOSED or not (unlike \c GetReciprocalPointOps, which is the ops-to-impose set and stays empty on a
    //! free run).  This is what the §3 order-parameter diagnostic (doc/SymmetryUpgradePlan.md) measures the
    //! converged density against: a free run reports the symmetry it actually found, per op
    //! (\c SymmetryDefects).  Default {} = no lattice / nothing detected.
    virtual std::vector<Symmetry::Lattice_3D::ReciprocalOp> GetDetectedReciprocalOps() const {return {};}

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