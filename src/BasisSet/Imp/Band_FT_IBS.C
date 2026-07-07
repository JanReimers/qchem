// File: BasisSet/Imp/Band_FT_IBS.C  Cached D-free reciprocal-space 3-centre accessors.
//
// The plane-wave mirror of Imp/Orbital_DFT_IBS.C: the public Repulsion3C/Overlap3C route through the
// process-wide integrals cache (theCache<dcmplx>(), keyed by BasisSetID), and the one-time build is the
// concrete basis's MakeRepulsion3C/MakeOverlap3C.  The tensor is intrinsic to the orbital's {G} grid, so it
// is valid for the basis's lifetime regardless of the density -- exactly the ERI3 caching contract.
module;
module qchem.BasisSet.Band_FT_IBS;
import qchem.BasisSet.Internal.DB_Cache;
import qchem.Types;   // dcmplx

namespace qchem::BasisSet
{

const G_ERI3& Band_FT_IBS::Repulsion3C(const cFIT_CD_ABS& c) const
{
    return theCache<dcmplx>().Get(IntegralsCache_Base::I3C::Repulsion,this,&c,
        [this,&c]{ return MakeRepulsion3C(c); });
}

const G_ERI3& Band_FT_IBS::Overlap3C(const cFIT_SF_ABS& c) const
{
    return theCache<dcmplx>().Get(IntegralsCache_Base::I3C::Overlap,this,&c,
        [this,&c]{ return MakeOverlap3C(c); });
}

} //namespace
