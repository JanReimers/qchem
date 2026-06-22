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

import Common.Iterators;
export import qchem.Streamable;

export namespace BasisSet
{
typedef std::vector<Irrep> irrepv_t; 

//----------------------------------------------------------------------------
//
//  Interface for a BasisSet which is assumed to a list of Irrep Basis Sets.
//
template <class T> class BasisSet
    : public virtual Streamable
{
public:
    typedef Orbital_1E_IBS<T> bs_t;

    virtual ~BasisSet() {};
    virtual size_t   GetNumFunctions() const=0;
    virtual irrepv_t GetIrreps(const Spin& ms) const=0;

    virtual Fit_IBS* CreateCDFitBasisSet(const Structure* cl) const;
    virtual Fit_IBS* CreateVxcFitBasisSet(const Structure* cl) const;

    // Iterate returns IBS pointers dynamic_cast'ed to the requested derived
    // type D.  It is built on the two primitives below, so how a concrete
    // BasisSet stores its IBSs stays private to that class.
    template <class D> auto Iterate() const
    {
        return D_IndexProxy<const D, BasisSet>(this, GetNumIBS());
    }
    const bs_t* operator[](size_t i) const {return GetIBS(i);}

protected:
    virtual size_t      GetNumIBS()    const=0; //Number of Irrep basis sets.
    virtual const bs_t* GetIBS(size_t) const=0; //The only storage-specific primitive.
};

typedef BasisSet<double>    Real_BS;
typedef BasisSet<dcmplx> Complex_BS;

}//namespace