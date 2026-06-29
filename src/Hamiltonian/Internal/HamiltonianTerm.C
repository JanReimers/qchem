// File: HamiltonianTerm.C  General implementation of a HamiltonianTerm term in the Hamiltonian.
//
// Templated on the matrix element type T (double for atoms/molecules; dcmplx for the plane-wave
// lattice lineage).  hmat_t<double> IS rsmat_t and tobs_t<double> IS obs_t, so the <double> aliases
// at the bottom leave existing real code unchanged.  The cache-lookup bodies are inline here (not in a
// separate Imp unit) because the dcmplx terms instantiate these templates in other translation units
// and so need the definitions visible.
module;
#include <map>
#include <cassert>
#include <cstddef>
export module qchem.Hamiltonian.Internal.Term;
export import qchem.Hamiltonian;
import qchem.Hamiltonian.Types;

import qchem.Symmetry.Irrep;

export namespace qchem::Hamiltonian
{

template <class T> class tHT_Common
{
protected:
    typedef std::map<Irrep,hmat_t<T>> CacheMap;
    mutable CacheMap   itsCache;       //Cache the H matrices for total energy calculations.
};


template <class T> class tStatic_HT_Imp
    : public virtual tStatic_HT<T>
    , protected tHT_Common<T>
{
public:
    virtual const hmat_t<T>& GetMatrix(const tobs_t<T>* bs,const Spin& s) const
    {
        assert(bs);
        Irrep qns(bs->GetIrrep(s));
        auto i=this->itsCache.find(qns);
        if (i==this->itsCache.end())
            return this->itsCache[qns]=CalculateMatrix(bs,s);
        else
            return i->second;
    }

protected:
    // Unconditional calculation, does not use cache.
    virtual hmat_t<T> CalculateMatrix(const tobs_t<T>*,const Spin&) const=0;
};

template <class T> class tDynamic_HT_Imp
    : public virtual tDynamic_HT<T>
    , protected tHT_Common<T>
{
public:
    tDynamic_HT_Imp() : itsCacheVersion(size_t(-1)), itsFitVersion(size_t(-1)) {}
    virtual const hmat_t<T>& GetMatrix(const tobs_t<T>* bs,const Spin& s,const tChargeDensity<T>* cd) const
    {
        assert(bs);
        assert(cd);
        // The cache keys on Irrep, not the density.  If the density has changed (its logical-clock serial
        // differs from the one the cache was built for), the cached matrices are stale -- drop them and
        // recompute for the new cd.  This makes GetMatrix self-correcting for ANY caller (the serial is
        // intrinsic to the density, stamped at construction/mutation), not just densities the SCFIterator
        // hands us -- which is why a plain density-pointer key, or an SCFIterator-only stamp, both failed.
        // This is now the SOLE owner of cache-freshness; newCD() only signals the concrete term to refit.
        if (cd->Version()!=itsCacheVersion)
        {
            this->itsCache.clear();
            itsCacheVersion=cd->Version();
        }
        Irrep qns(bs->GetIrrep(s));
        if (auto i=this->itsCache.find(qns);i==this->itsCache.end())
            return this->itsCache[qns]=CalcMatrix(bs,s,cd);
        else
            return i->second; //Cache hit (same density serial, already computed for this Irrep)
    }

protected:
    // Unconditional calculation, does not use cache.
    virtual hmat_t<T> CalcMatrix(const tobs_t<T>*,const Spin&,const tChargeDensity<T>* cd) const=0;
    // Refit trigger only (the OTHER half of the former itsCD double-duty): true exactly once per new
    // density serial, so a concrete term (FittedVee/FittedVxc) refits its fitted potential exactly once.
    // Cache-freshness is GetMatrix's job (itsCacheVersion) -- newCD no longer touches the cache.
    bool newCD(const tChargeDensity<T>* cd) const
    {
        assert(cd);
        if (cd->Version()==itsFitVersion) return false;
        itsFitVersion=cd->Version();
        return true;
    }

    mutable size_t itsCacheVersion;   //!< density serial the Irrep cache was built for
    mutable size_t itsFitVersion;     //!< density serial the concrete term last refit for
};

// Used for polarized potentials (Vxc) which each polarization will handle its own cache.
template <class T> class tDynamic_HT_Imp_NoCache
    : public virtual tDynamic_HT<T>
{
public:
    virtual const hmat_t<T>& GetMatrix(const tobs_t<T>* bs,const Spin& s,const tChargeDensity<T>* cd) const
    {
        return itsMat=CalcMatrix(bs,s,cd);
    }

protected:
    virtual hmat_t<T> CalcMatrix(const tobs_t<T>*,const Spin&,const tChargeDensity<T>* cd) const=0;
    mutable hmat_t<T> itsMat;
};

// A density-dependent Hamiltonian potential term.  (No longer a Fitting::ScalarFFClient -- the
// "I can answer the fitter's questions" role belongs on the concrete LDAVxc, not this term interface.)
template <class T> class tFittablePotential
    : public virtual tDynamic_HT<T>
{
public:
    virtual void UseChargeDensity(const tChargeDensity<T>*)       =0;
};

// r* = <double>, c* = <dcmplx> (mirrors rsmat_t/chmat_t); bare names transitional (= r*), rename pinned.
using rHT_Common              = tHT_Common<double>;              using cHT_Common              = tHT_Common<dcmplx>;
using rStatic_HT_Imp          = tStatic_HT_Imp<double>;          using cStatic_HT_Imp          = tStatic_HT_Imp<dcmplx>;
using rDynamic_HT_Imp         = tDynamic_HT_Imp<double>;         using cDynamic_HT_Imp         = tDynamic_HT_Imp<dcmplx>;
using rDynamic_HT_Imp_NoCache = tDynamic_HT_Imp_NoCache<double>; using cDynamic_HT_Imp_NoCache = tDynamic_HT_Imp_NoCache<dcmplx>;
using rFittablePotential      = tFittablePotential<double>;      using cFittablePotential      = tFittablePotential<dcmplx>;
using HT_Common              = rHT_Common;
using Static_HT_Imp          = rStatic_HT_Imp;
using Dynamic_HT_Imp         = rDynamic_HT_Imp;
using Dynamic_HT_Imp_NoCache = rDynamic_HT_Imp_NoCache;
using FittablePotential      = rFittablePotential;

} //namespace
