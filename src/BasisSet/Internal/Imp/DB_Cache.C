// File: BasisSet/Internal/Imp/DB_Cache.C
// Definition of the construct-on-first-use accessor theCache<T>() declared in the DB_Cache module
// interface.  It lives here (an implementation unit of the SAME module) rather than in the interface
// so the single function-static instance is emitted in exactly one TU; importing DB_Cache_RAM gives
// us the complete IntegralsCache_RAM<T> needed to build it.
module;
module qchem.BasisSet.Internal.DB_Cache;
import qchem.BasisSet.Internal.DB_Cache_RAM;

namespace qchem::BasisSet
{

template <class T> IntegralsCache<T>& theCache()
{
    static IntegralsCache_RAM<T> instance(/*makelog=*/false);
    return instance;
}

template IntegralsCache<double>& theCache<double>();
// Complex (plane-wave / Bloch) cache: Hermitian (chmat_t) I2C Overlap etc.
template IntegralsCache<dcmplx>& theCache<dcmplx>();

} //namespace
