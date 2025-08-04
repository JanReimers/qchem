// File: Irrep_BS.C  Interface for an Irrep Basis Set
module;
#include <vector>
#include <memory>

export module qchem.Irrep_BS;
export import qchem.BasisFunction;
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
//  block for  IrrepBasisSet,  For atoms each L get.s an IrrepBasisSet and an H  block. 
//  
//  The quantum number could be L for atoms, Irreducable rep for molecules, or
//  the wave vector k for solids.
//

export class IrrepBasisSet
    : public virtual UniqueID
    , public virtual Streamable
{
public:
    typedef Real_BF bf_t;
    typedef std::vector<std::shared_ptr<bf_t>> bfv_t;
    typedef bfv_t::      iterator       iterator;
    typedef bfv_t::const_iterator const_iterator;
    typedef std::shared_ptr<const Symmetry> sym_t;

    virtual size_t GetNumFunctions() const=0;
    virtual sym_t  GetSymmetry() const=0;
//
//  Streamable stuff.
//
    virtual IrrepBasisSet* Clone  (const RVec3&) const=0;

};

//----------------------------------------------------------------------------
//
//  Extend basis to be a set of real or complex valued functions
//
export template <class T> class TIrrepBasisSet
    : public virtual IrrepBasisSet
    , public virtual VectorFunction<T>
{
public:
    size_t GetVectorSize() const {return GetNumFunctions();}
};

//
// Define an orbital irrep basis set which supports integrals for SCF orbital calculations.
// Mix-in the integral interfaces required for an orbital basis. 
//
export template <class T> class TOrbital_IBS
    : public virtual TIrrepBasisSet<T>
    , public virtual Integrals_Overlap<T> 
    , public virtual Integrals_Kinetic<T> 
    , public virtual Integrals_Nuclear<T> 
{
public:    
    virtual void   Set(const LAParams&)=0;
    virtual LASolver<T>* CreateSolver() const=0;
};
