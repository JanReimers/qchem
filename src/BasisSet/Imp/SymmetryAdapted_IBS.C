// File: BasisSet/Imp/SymmetryAdapted_IBS.C
module;
#include <string>
#include <iostream>
#include <memory>
#include <cassert>
module qchem.BasisSet.SymmetryAdapted_IBS;
import qchem.Blaze;          // trans, submatrix, matrix/vector ops

namespace qchem::BasisSet
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
const rsmat_t& SymFockCache::Direct(const Orbital_HF_IBS<double>* raw, const Orbital_HF_IBS<double>* cd,
                                    const rmat_t& Ocd, const rsmat_t& Dcd)
{
    Entry& e = itsJ[cd->BasisSetID()];
    if (!e.valid || blazem::max(blazem::abs(e.D - Dcd)) > 0.0)
    {
        e.M = BuildAOFock(false, raw, Ocd, Dcd);
        e.D = Dcd; e.valid = true;
    }
    return e.M;
}
const rsmat_t& SymFockCache::Exchange(const Orbital_HF_IBS<double>* raw, const Orbital_HF_IBS<double>* cd,
                                      const rmat_t& Ocd, const rsmat_t& Dcd)
{
    Entry& e = itsK[cd->BasisSetID()];
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
    , itsRawHF (dynamic_cast<const Orbital_HF_IBS <double>*>(raw))  // same object, HF  face
    , itsRawDFT(dynamic_cast<const Orbital_DFT_IBS<double>*>(raw))  // same object, DFT face
    , itsO(Oblock), itsLabel(label), itsCache(cache)
{}

// O^T Mraw O, explicitly symmetrized (the product is symmetric only up to roundoff).
rsmat_t SymmetryAdapted_IBS::Transform(const rsmat_t& Mraw) const
{
    return SymCopy(blazem::trans(itsO) * Mraw * itsO);   // dGamma x dGamma
}

// Transform each fit-function matrix of a 3C tensor: O^T m O (linear in the AO 3C, like J/K).
ERI3<double> SymmetryAdapted_IBS::TransformERI3(const ERI3<double>& raw) const
{
    ERI3<double> out; out.reserve(raw.size());
    for (const auto& m : raw) out.push_back(Transform(m));         // dGamma x dGamma per fit fn
    return out;
}

// DFT 3-centre integrals in the irrep basis.  These are the cache-miss compute hooks invoked by
// the inherited Overlap3C/Repulsion3C accessors; they transform the raw basis's *cached* 3C (the
// integral cache is re-entrant now, so the nested cached access is safe).  The raw 3C is therefore
// computed once and shared across all irreps; only the cheap O^T (.) O transform is per irrep.
ERI3<double> SymmetryAdapted_IBS::MakeOverlap3C  (const FIT_SF_ABS& c) const { return TransformERI3(itsRawDFT->Overlap3C(c));   }
ERI3<double> SymmetryAdapted_IBS::MakeRepulsion3C(const FIT_CD_ABS& c) const { return TransformERI3(itsRawDFT->Repulsion3C(c)); }

// Fit bases are atom-centred (geometry, not symmetry), so reuse the raw basis's unchanged.
FIT_CD_ABS* SymmetryAdapted_IBS::CreateCDFitBasisSet (const Structure* cl, const qcMesh::MeshParams& mp) const { return itsRawDFT->CreateCDFitBasisSet(cl,mp);  }
FIT_SF_ABS* SymmetryAdapted_IBS::CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const { return itsRawDFT->CreateVxcFitBasisSet(cl,mp); }

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

// Transform the raw basis's *cached* 1-e matrix (Overlap/Kinetic/Nuclear): the raw matrix is
// computed once and shared by every irrep.  The nested cached access is safe now that the integral
// cache is re-entrant (these MakeXxx run as the cache-miss hook of this irrep's own cached block).
rsmat_t SymmetryAdapted_IBS::MakeOverlap()                 const {return Transform(itsRaw->Overlap());}
// Kinetic() here is the \f$\langle p^2\rangle=\langle-\nabla^2\rangle\f$ block (no 1/2), just transformed
// into the SALC basis -- the 1/2 / energy meaning is applied later (see BasisSet/Orbital_1E_IBS.C).
rsmat_t SymmetryAdapted_IBS::MakeKinetic()                 const {return Transform(itsRaw->Kinetic());}
rsmat_t SymmetryAdapted_IBS::MakeNuclear(const Structure* cl) const {return Transform(itsRaw->Nuclear(cl));}

std::string SymmetryAdapted_IBS::RadialID()   const {return itsRaw->RadialID();}
std::string SymmetryAdapted_IBS::AngularID()  const {return itsRaw->AngularID();}
std::string SymmetryAdapted_IBS::BasisSetID() const {return itsRaw->BasisSetID() + "[" + itsLabel + "]";}
std::string SymmetryAdapted_IBS::Name()       const {return itsRaw->Name() + "[" + itsLabel + "]";}

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
