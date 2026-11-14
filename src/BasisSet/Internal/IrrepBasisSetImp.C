// File: BasisSet/Internal/IrrepBasisSetImp.C Implement a generic IrrepBasisSet
module;
#include <cassert>
#include <iostream>
#include <memory>
#include <stdexcept>
export module qchem.BasisSet.Internal.IrrepBasisSetImp;
import qchem.BasisSet.IrrepBasisSet;
import qchem.Symmetry.Irrep;
export import qchem.BasisSet.BasisCache;   // every concrete block's constructor names it

export namespace qchem::BasisSet
{

template <class T> class IrrepBasisSetImp
    : public virtual IrrepBasisSet<T>
{
public:
    //! \a cache is the PARENT's cache (doc/IntegralCacheOwnershipPlan.md S3a): the parent BasisSet's factory makes it and passes it to
    //! every block it builds; a unit test holds one in its harness.  Never null.
    IrrepBasisSetImp(const sym_t& sym, const std::shared_ptr<BasisCache>& cache)
        : itsSymmetry(sym), itsCache(cache), itsCacheIndex(Enrol(cache)) {assert(itsSymmetry);}
    //! A copy is a NEW block of the same parent: same symmetry and cache, its OWN index (and its own instance token), so it never reads
    //! the original's slots (a pruned subset has fewer functions).
    IrrepBasisSetImp(const IrrepBasisSetImp& o) : itsSymmetry(o.itsSymmetry), itsCache(o.itsCache), itsCacheIndex(Enrol(o.itsCache)) {}
    IrrepBasisSetImp& operator=(const IrrepBasisSetImp&) = delete;
    //! Readonly ref to the polymorphic Symmetry object.
    virtual const Symmetry::Symmetry& GetSymmetry() const {return *itsSymmetry;}
    virtual const           sym_t   & GetSymt    () const {return  itsSymmetry;}
    //! Irrep basis sets are spin agnostic, so caller must specify the spin in order to a full set of QNs.
    virtual Irrep GetIrrep(const Spin& s) const
    {
        return Irrep(s,itsSymmetry);
    }
    
    virtual const InstanceToken& Instance() const override {return itsInstance;}
protected:
    virtual const BasisCache& Cache     () const override {return *itsCache;}
    virtual size_t      CacheIndex() const override {return itsCacheIndex;}
private:
    static size_t Enrol(const std::shared_ptr<BasisCache>& cache)
    {
        if (!cache) throw std::logic_error("IrrepBasisSetImp: a block was built without its parent's cache (the factory must pass one)");
        return cache->Enrol();
    }
    sym_t                       itsSymmetry;
    std::shared_ptr<BasisCache> itsCache;
    size_t                      itsCacheIndex;
    InstanceToken               itsInstance;
};

} //namespace