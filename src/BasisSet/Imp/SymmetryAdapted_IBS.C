// File: BasisSet/Imp/SymmetryAdapted_IBS.C
module;
#include <string>
#include <iostream>
#include <memory>
#include <cassert>
module qchem.BasisSet.SymmetryAdapted_IBS;
import qchem.Blaze;          // trans, submatrix, matrix/vector ops

namespace BasisSet
{

// Symmetrize a (numerically near-symmetric) square matrix into an rsmat_t.
static rsmat_t SymCopy(const rmat_t& P)
{
    size_t n = P.rows();
    rsmat_t S(n);
    for (size_t i=0;i<n;i++) for (size_t j=0;j<=i;j++) S(i,j) = 0.5*(P(i,j)+P(j,i));
    return S;
}

// Build the AO Coulomb (exch=false) or exchange (true) from a cd-irrep density block: back-
// transform Ocd Dcd Ocd^T to AO, then the raw whole-molecule AccumulateDirect/Exchange (reuses
// the cached AO ERIs -- no 4-index transform).
static rsmat_t BuildAOFock(bool exch, const Orbital_HF_IBS<double>* raw,
                           const rmat_t& Ocd, const rsmat_t& Dcd)
{
    rsmat_t Dao = SymCopy(Ocd * Dcd * blazem::trans(Ocd));
    rsmat_t M   = blazem::zero<double>(Ocd.rows());
    if (blazem::max(blazem::abs(Dao)) > 0.0)              // pre-screen (raw asserts non-zero)
    {
        if (exch) raw->AccumulateExchange(M, Dao, raw);
        else      raw->AccumulateDirect  (M, Dao, raw);
    }
    return M;
}
// Fab += O^T Mao O (symmetrized): slice an AO matrix down to this irrep's block.
static void AddSlice(rsmat_t& Fab, const rmat_t& O, const rsmat_t& Mao)
{
    rmat_t B = blazem::trans(O) * Mao * O;
    for (size_t i=0;i<B.rows();++i) for (size_t j=0;j<=i;++j) Fab(i,j) += 0.5*(B(i,j)+B(j,i));
}

// Memoize the AO J/K for a cd-irrep, rebuilding only when that block's density changes.
const rsmat_t& SymFockCache::Direct(const Orbital_HF_IBS<double>* raw, const void* cd,
                                    const rmat_t& Ocd, const rsmat_t& Dcd)
{
    Entry& e = itsJ[cd];
    if (!e.valid || blazem::max(blazem::abs(e.D - Dcd)) > 0.0)
    {
        e.M = BuildAOFock(false, raw, Ocd, Dcd);
        e.D = Dcd; e.valid = true;
    }
    return e.M;
}
const rsmat_t& SymFockCache::Exchange(const Orbital_HF_IBS<double>* raw, const void* cd,
                                      const rmat_t& Ocd, const rsmat_t& Dcd)
{
    Entry& e = itsK[cd];
    if (!e.valid || blazem::max(blazem::abs(e.D - Dcd)) > 0.0)
    {
        e.M = BuildAOFock(true, raw, Ocd, Dcd);
        e.D = Dcd; e.valid = true;
    }
    return e.M;
}

SymmetryAdapted_IBS::SymmetryAdapted_IBS(const Orbital_1E_IBS<double>* raw, const rmat_t& Oblock,
                                         const std::string& label, const sym_t& sym,
                                         std::shared_ptr<SymFockCache> cache)
    : IrrepBasisSetImp<double>(sym), itsRaw(raw)
    , itsRawHF(dynamic_cast<const Orbital_HF_IBS<double>*>(raw))   // same object, HF face
    , itsO(Oblock), itsLabel(label), itsCache(cache)
{}

// O^T Mraw O, explicitly symmetrized (the product is symmetric only up to roundoff).
rsmat_t SymmetryAdapted_IBS::Transform(const rsmat_t& Mraw) const
{
    return SymCopy(blazem::trans(itsO) * Mraw * itsO);   // dGamma x dGamma
}

// Coulomb / exchange are linear in the density, so build the AO matrix from the cd-irrep's
// density block and slice to this irrep (no 4-index ERI transform).  The AO build depends only
// on the cd-irrep, so a shared cache builds it once per iteration (N instead of N^2 builds);
// without a cache it is built directly each call.
void SymmetryAdapted_IBS::AccumulateDirect(rsmat_t& Jab, const rsmat_t& Dcd,
                                           const Orbital_HF_IBS<double>* bs_cd) const
{
    const SymmetryAdapted_IBS* cd = dynamic_cast<const SymmetryAdapted_IBS*>(bs_cd);
    assert(cd && itsRawHF);
    if (itsCache) AddSlice(Jab, itsO, itsCache->Direct(itsRawHF, bs_cd, cd->itsO, Dcd));
    else          AddSlice(Jab, itsO, BuildAOFock(false, itsRawHF, cd->itsO, Dcd));
}

void SymmetryAdapted_IBS::AccumulateExchange(rsmat_t& Kab, const rsmat_t& Dcd,
                                             const Orbital_HF_IBS<double>* bs_cd) const
{
    const SymmetryAdapted_IBS* cd = dynamic_cast<const SymmetryAdapted_IBS*>(bs_cd);
    assert(cd && itsRawHF);
    if (itsCache) AddSlice(Kab, itsO, itsCache->Exchange(itsRawHF, bs_cd, cd->itsO, Dcd));
    else          AddSlice(Kab, itsO, BuildAOFock(true, itsRawHF, cd->itsO, Dcd));
}

// Use the raw COMPUTE (MakeX), not the cached accessor (X()): the decorator's own cached
// accessor is mid Has/Set when MakeOverlap runs, and a nested cache Has/Set on the raw would
// clobber the cache's stateful "last key" and mis-file the transformed block.  (So the raw
// 1-e matrix is recomputed per irrep; cheap, and the transformed block is still cached.)
rsmat_t SymmetryAdapted_IBS::MakeOverlap()                 const {return Transform(itsRaw->MakeOverlap());}
rsmat_t SymmetryAdapted_IBS::MakeKinetic()                 const {return Transform(itsRaw->MakeKinetic());}
rsmat_t SymmetryAdapted_IBS::MakeNuclear(const Cluster* cl) const {return Transform(itsRaw->MakeNuclear(cl));}

std::string SymmetryAdapted_IBS::RadialID()  const {return itsRaw->RadialID();}
std::string SymmetryAdapted_IBS::AngularID() const {return itsRaw->AngularID() + "_" + itsLabel;}
std::string SymmetryAdapted_IBS::Name()      const {return itsRaw->Name() + "[" + itsLabel + "]";}

rvec_t SymmetryAdapted_IBS::operator()(const rvec3_t& r) const
{
    return blazem::trans(itsO) * (*itsRaw)(r);     // dGamma SALC values from nAO AO values
}

rvec3vec_t SymmetryAdapted_IBS::Gradient(const rvec3_t& r) const
{
    rvec3vec_t g = itsRaw->Gradient(r);            // nAO gradients
    size_t nAO = itsO.rows(), dG = itsO.columns();
    rvec3vec_t out(dG);
    for (size_t a=0;a<dG;a++)
    {
        rvec3_t s(0,0,0);
        for (size_t i=0;i<nAO;i++) s += itsO(i,a)*g[i];
        out[a] = s;
    }
    return out;
}

std::ostream& SymmetryAdapted_IBS::Write(std::ostream& os) const
{
    return os << "SymmetryAdapted_IBS[" << itsLabel << "] dim=" << GetNumFunctions();
}

} //namespace
