// File: Atom/radial/Slater/BS_Common.H  l/ml/kappa/mj independent part of BasisSet for atom Slater Basis Sets.
module;
#include <vector>
#include <iosfwd>
export module qchem.BasisSet.Atom.Internal.radial.SlaterBS;
import qchem.BasisSet.Internal.IrrepBasisSet;
import qchem.BasisSet.Atom.IEClient;
import qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Internal.Common;
export import qchem.BasisSet;
export import qchem.Types;

export namespace Slater
{
// Common implementation for orbital and fit basis sets.
class IrrepBasisSet
    : public virtual Real_IBS
    , public         IrrepBasisSet_Common<double>
{
    typedef typename VectorFunction<double>::Vec     Vec;  //Vector of scalars.
    typedef typename VectorFunction<double>::Vec3Vec Vec3Vec;//vector of 3 space vectors.
public:
    IrrepBasisSet(IBS_Evaluator*, Symmetry*);
    virtual size_t  GetNumFunctions() const {return itsEval->size();}
    virtual size_t  size() const {return itsEval->size();}
    virtual Vec     operator() (const RVec3& r) const {return itsEval->operator()(r);}
    virtual Vec3Vec Gradient   (const RVec3& r) const {return itsEval->Gradient(r);}
    virtual std::ostream&  Write(std::ostream&    ) const;
private:
    IBS_Evaluator* itsEval;
};

// Creates the Rk tool for HF ERIs
class BS_Common
: public ::BS_Common
, public ::AtomIE_BS_2E<double>
{
protected:
    BS_Common(BS_Evaluator* bse) : AtomIE_BS_2E<double>(bse) {};
    virtual void Insert(bs_t* bs);
};

}

