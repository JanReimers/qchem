// File: Atom/radial/Slater/BS_Common.H  l/ml/kappa/mj independent part of BasisSet for atom Slater Basis Sets.
module;
#include <vector>
#include <iosfwd>
export module qchem.BasisSet.Atom.Internal.radial.SlaterBS;
import qchem.BasisSet.Internal.IBS_Common;
import qchem.BasisSet.Atom.IEClient;
import qchem.BasisSet.Atom.IE;
import qchem.BasisSet.Internal.Common;
export import qchem.Types;

export namespace Slater
{
// Common implementation for orbital and fit basis sets.
class IrrepBasisSet
    : public virtual ::IrrepBasisSet
    , public         TIBS_Common1<double>
    , public         AtomIrrepIEClient
{
    typedef typename VectorFunction<double>::Vec     Vec;  //Vector of scalars.
    typedef typename VectorFunction<double>::Vec3Vec Vec3Vec;//vector of 3 space vectors.
public:
    IrrepBasisSet(const Vector<double>& exponents, Symmetry*,size_t L);
    IrrepBasisSet(const Vector<double>& exponents, Symmetry*,size_t L, const std::vector<int>& ml);
    virtual size_t  GetNumFunctions() const {return size();}
    virtual Vec     operator() (const RVec3&) const;
    virtual Vec3Vec Gradient   (const RVec3&) const;
    virtual std::ostream&  Write(std::ostream&    ) const;
private:
    Vector<double> Norms(const Vector<double>& exponents, size_t l) const;
};

// Common base handles all the radial aspects.
class BS_Common
: public ::BS_Common
, public ::AtomIE_BS_2E<double>
{
protected:
    virtual void Insert(bs_t* bs);
private:
    virtual const Cacheable* Create(size_t ia,size_t ic,size_t ib,size_t id) const;
    virtual Vector<double> loop_4_direct  (size_t id, size_t la, size_t lc)  const;
    virtual Vector<double> loop_4_exchange(size_t id, size_t la, size_t lc)  const;
};

}

