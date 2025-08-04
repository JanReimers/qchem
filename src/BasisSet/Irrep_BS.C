// File: Irrep_BS.C  Interface for an Irrep Basis Set
module;
#include <memory>

export module qchem.Irrep_BS;
export import qchem.Symmetry;
export import qchem.VectorFunction;

import qchem.BasisSet.Internal.Integrals;
import qchem.LASolver;
import Common.UniqueID; 
import qchem.Streamable;

//! \brief Interface for overlap integrals.
export template <class T> class Integrals_Overlap
{
public:
    //! Single basis set Overlap \f$ \left\langle a\left|1\right|b\right\rangle =\int d^{3}\vec{r}\:g_{a}\left(\vec{r}\right)g_{b}\left(\vec{r}\right) \f$ 
    virtual const SMatrix<T>& Overlap() const=0;
};
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
    , public virtual Integrals_Overlap<T>
{
public:
    typedef std::shared_ptr<const Symmetry> sym_t;
    virtual sym_t  GetSymmetry() const=0;
    
    virtual size_t GetNumFunctions() const=0;
    size_t GetVectorSize() const {return GetNumFunctions();}
};

export typedef IrrepBasisSet<double>    Real_IBS;
export typedef IrrepBasisSet<dcmplx> Complex_IBS;

