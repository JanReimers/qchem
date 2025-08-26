// File: Atom/radial/Gaussian/BS_Common.H  l/ml/kappa/mj independent part of BasisSet for atom Gaussians.
module;
#include <vector>
#include <iosfwd>
export module qchem.BasisSet.Atom.Internal.radial.GaussianBS;
import qchem.BasisSet.Internal.IrrepBasisSet;
import qchem.BasisSet.Atom.IEClient;
import qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Internal.Common;

export namespace Gaussian
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
    virtual size_t  GetNumFunctions() const {return size();}
    virtual size_t  size() const {return itsEval->size();}
    virtual Vec     operator() (const RVec3& r) const {return itsEval->operator()(r);}
    virtual Vec3Vec Gradient   (const RVec3& r) const {return itsEval->Gradient(r);}
    virtual std::ostream&  Write(std::ostream&    ) const;
private:
    IBS_Evaluator* itsEval;
};
// Common base handles all the radial aspects.
class BS_Common
: public ::BS_Common
, public ::AtomIE_BS_2E<double>
{
protected:
    BS_Common(BS_Evaluator* bse) : AtomIE_BS_2E<double>(bse) {};
    virtual void Insert(bs_t* bs);
};

}

