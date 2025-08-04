// File: Irrep_BS.C  Interface for an Irrep Basis Set
module;
#include <vector>
#include <memory>

export module qchem.Irrep_BS;
export import qchem.Symmetry;
export import qchem.Symmetry.ElectronConfiguration;
export import qchem.VectorFunction;

import qchem.BasisSet.Internal.Integrals;
import qchem.LASolver;
import Common.UniqueID; 
import Common.Iterators;
import qchem.Streamable;
import oml;

//----------------------------------------------------------------------------
//
//  Interface for an irreducible representation basis sets.  H is block diagonal with one
//  block for each IrrepBasisSet,  For atoms each L gets an IrrepBasisSet and an H  block. 
//  The Symmetry could be spherical (l,m QNs) for atoms, point group for molecules, or
//  transational (with wave vector k) for solids.
//
export template <class T> class IrrepBasisSet
    : public virtual UniqueID
    , public virtual Streamable
    , public virtual VectorFunction<T>
{
public:
    typedef std::shared_ptr<const Symmetry> sym_t;
    virtual sym_t  GetSymmetry() const=0;
    
    virtual size_t GetNumFunctions() const=0;
    size_t GetVectorSize() const {return GetNumFunctions();}
};

export typedef IrrepBasisSet<double>    Real_IBS;
export typedef IrrepBasisSet<dcmplx> Complex_IBS;

//
// Define an orbital irrep basis set which supports integrals for SCF orbital calculations.
// Mix-in the integral interfaces required for an orbital basis. 
//
export template <class T> class Orbital_IBS
    : public virtual IrrepBasisSet<T>
    , public virtual Integrals_Overlap<T> 
    , public virtual Integrals_Kinetic<T> 
    , public virtual Integrals_Nuclear<T> 
{
public:    
    virtual void         Set(const LAParams&)=0;
    virtual LASolver<T>* CreateSolver() const=0;
};

export typedef Orbital_IBS<double>    Real_OIBS;
export typedef Orbital_IBS<dcmplx> Complex_OIBS;
