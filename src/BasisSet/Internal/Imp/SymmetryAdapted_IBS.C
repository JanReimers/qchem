// File: BasisSet/Internal/Imp/SymmetryAdapted_IBS.C
module;
#include <vector>
#include <memory>
#include <string>
#include <iostream>
#include <memory>
#include <cassert>
#include <stdexcept>
module qchem.BasisSet.Internal.SymmetryAdapted_IBS;
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

// --- WholeSystemFock_IBS (V1.31) -----------------------------------------------------------------
// The three primitives the composite density drives instead of the per-irrep-pair loop.  Between them they
// say the whole physics of this class's 2-electron path: fold each block's density up to AO, build ONCE,
// slice back down.  J and K are both LINEAR in D, which is what makes the single build exact rather than an
// approximation -- sum(J_AO(O_C D_C O_C^T)) == J_AO(sum(O_C D_C O_C^T)).
void SymmetryAdapted_IBS::AddAODensity(rsmat_t& Dao, const rsmat_t& D) const
{
    Dao += SymCopy(itsO * D * blazem::trans(itsO));
}

rsmat_t SymmetryAdapted_IBS::MakeAOFock(const rsmat_t& Dao, bool exchange) const
{
    assert(itsRawHF);
    rsmat_t M = blazem::zero<double>(itsO.rows());
    if (blazem::max(blazem::abs(Dao)) > 0.0)          // pre-screen (the raw basis asserts non-zero)
    {
        if (exchange) itsRawHF->AccumulateExchange(M, Dao, itsRawHF);
        else          itsRawHF->AccumulateDirect  (M, Dao, itsRawHF);
    }
    return M;
}

void SymmetryAdapted_IBS::SliceAOFock(rsmat_t& Fab, const rsmat_t& Fao) const
{
    AddSlice(Fab, itsO, Fao);
}

SymmetryAdapted_IBS::SymmetryAdapted_IBS(const Orbital_1E_IBS<double>* raw, const rmat_t& Oblock,
                                         const std::string& label, const sym_t& sym,
                                         std::shared_ptr<const std::vector<Symmetry::SymOp>> group)
    : IrrepBasisSetImp<double>(sym)
    , Symmetry::IdentityRepresentation<double>(group ? Symmetry::Stabiliser::OpList(std::move(group)) : Symmetry::Stabiliser())
    , itsRaw(raw)
    , itsRawHF (dynamic_cast<const Orbital_HF_IBS <double>*>(raw))  // same object, HF  face
    , itsRawDFT(dynamic_cast<const Orbital_DFT_IBS<double>*>(raw))  // same object, DFT face
    , itsO(Oblock), itsLabel(label)
{}

// O^T Mraw O, explicitly symmetrized (the product is symmetric only up to roundoff).
rsmat_t SymmetryAdapted_IBS::Transform(const rsmat_t& Mraw) const
{
    return SymCopy(blazem::trans(itsO) * Mraw * itsO);   // dGamma x dGamma
}

// Transform each fit-function matrix of a 3C tensor: O^T m O (linear in the AO 3C, like J/K).
PairProjection<double> SymmetryAdapted_IBS::Transform3C(const PairProjection<double>& raw) const
{
    PairProjection<double> out; out.dense.reserve(raw.dense.size());
    for (const auto& m : raw.dense) out.dense.push_back(Transform(m));         // dGamma x dGamma per fit fn
    return out;
}

// DFT 3-centre integrals in the irrep basis.  These are the cache-miss compute hooks invoked by
// the inherited Overlap3C/Repulsion3C accessors; they transform the raw basis's *cached* 3C (the
// integral cache is re-entrant now, so the nested cached access is safe).  The raw 3C is therefore
// computed once and shared across all irreps; only the cheap O^T (.) O transform is per irrep.
PairProjection<double> SymmetryAdapted_IBS::MakeOverlap3C  (const rFIT_SF_ABS& c) const { return Transform3C(itsRawDFT->Overlap3C(c));   }
PairProjection<double> SymmetryAdapted_IBS::MakeRepulsion3C(const rFIT_CD_ABS& c) const { return Transform3C(itsRawDFT->Repulsion3C(c)); }

// Fit bases are atom-centred (geometry, not symmetry), so reuse the raw basis's unchanged.
rFIT_CD_ABS* SymmetryAdapted_IBS::CreateCDFitBasisSet (const Structure* cl, const qcMesh::MeshParams& mp) const { return itsRawDFT->CreateCDFitBasisSet(cl,mp);  }
rFIT_SF_ABS* SymmetryAdapted_IBS::CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const { return itsRawDFT->CreateVxcFitBasisSet(cl,mp); }

// Coulomb / exchange are linear in the density, so build the AO matrix from the cd-irrep's
// density block and slice to this irrep (no 4-index ERI transform).  The AO build depends only
// on the cd-irrep, so a shared cache builds it once per iteration (N instead of N^2 builds);
// without a cache it is built directly each call.
// V1.10: no cast to the partner's CONCRETE type.  These used to reach into a sibling for its private itsO
// (its SALC columns) purely to fold the partner's density up to AO.  That operation already exists as an
// abstract question -- WholeSystemFock_IBS::AddAODensity, added by V1.31 -- so the three steps are now three
// calls on the face: the PARTNER folds its own block up, we build once, we slice into ours.  Nobody hands
// out their columns, and the cast is abstract->abstract (the sanctioned direction).
rsmat_t SymmetryAdapted_IBS::PairAOFock(const rsmat_t& Dcd, const Orbital_HF_IBS<double>* bs_cd,
                                        bool exchange) const
{
    auto* ws = dynamic_cast<const WholeSystemFock_IBS<double>*>(bs_cd);
    if (!ws) throw std::runtime_error("SymmetryAdapted_IBS: the HF partner basis has no whole-system Fock "
                                      "face, so its density block cannot be folded up to the AO space.");
    assert(itsRawHF);
    rsmat_t Dao = blazem::zero<double>(AODimension());
    ws->AddAODensity(Dao, Dcd);              // Dao = O_cd D_cd O_cd^T -- the partner uses its OWN columns
    return MakeAOFock(Dao, exchange);
}

void SymmetryAdapted_IBS::AccumulateDirect(rsmat_t& Jab, const rsmat_t& Dcd,
                                           const Orbital_HF_IBS<double>* bs_cd) const
{
    SliceAOFock(Jab, PairAOFock(Dcd, bs_cd, false));
}

void SymmetryAdapted_IBS::AccumulateExchange(rsmat_t& Kab, const rsmat_t& Dcd,
                                             const Orbital_HF_IBS<double>* bs_cd) const
{
    SliceAOFock(Kab, PairAOFock(Dcd, bs_cd, true));
}

// No per-irrep-pair ERI4 here (see the header): serve the bra-ket partner as two independent AO slices.
// Ji=irrep i's block (this=irrep i's SALC basis); Jj=irrep j's block (cd=irrep j's SALC basis).  The
// per-side zero guard mirrors the density's pre-screen and skips the (cheap) O^T(.)O transform on an
// empty block.
void SymmetryAdapted_IBS::AccumulateDirectBoth(rsmat_t& Ji, rsmat_t& Jj, const rsmat_t& Di, const rsmat_t& Dj,
                                               const Orbital_HF_IBS<double>* cd) const
{
    if (blazem::max(blazem::abs(Dj))>0.0)     AccumulateDirect(Ji, Dj, cd);   // Ji += J(i,j)·Dj
    if (blazem::max(blazem::abs(Di))>0.0) cd->AccumulateDirect(Jj, Di, this); // Jj += J(j,i)·Di
}

// Exchange counterpart: same AO-slice fallback (no per-irrep-pair ERI4 in the SALC path).
void SymmetryAdapted_IBS::AccumulateExchangeBoth(rsmat_t& Ki, rsmat_t& Kj, const rsmat_t& Di, const rsmat_t& Dj,
                                                 const Orbital_HF_IBS<double>* cd) const
{
    if (blazem::max(blazem::abs(Dj))>0.0)     AccumulateExchange(Ki, Dj, cd);   // Ki += K(i,j)·Dj
    if (blazem::max(blazem::abs(Di))>0.0) cd->AccumulateExchange(Kj, Di, this); // Kj += K(j,i)·Di
}

// Transform the raw basis's *cached* 1-e matrix (Overlap/Kinetic/Nuclear): the raw matrix is
// computed once and shared by every irrep.  The nested cached access is safe now that the integral
// cache is re-entrant (these MakeXxx run as the cache-miss hook of this irrep's own cached block).
rsmat_t SymmetryAdapted_IBS::MakeOverlap()                 const {return Transform(itsRaw->Overlap());}
// Kinetic() here is the \f$\langle p^2\rangle=\langle-\nabla^2\rangle\f$ block (no 1/2), just transformed
// into the SALC basis -- the 1/2 / energy meaning is applied later (see BasisSet/Orbital_1E_IBS.C).
rsmat_t SymmetryAdapted_IBS::MakeKinetic()                 const {return Transform(itsRaw->Kinetic());}
rsmat_t SymmetryAdapted_IBS::MakeNuclear(const Structure* cl) const {return Transform(itsRaw->Nuclear(cl));}

bool SymmetryAdapted_IBS::CanAssembleField(const SpeciesRadialField& f) const
{
    const auto* pp=dynamic_cast<const Orbital_PP_IBS<double>*>(itsRaw);
    return pp && pp->CanAssembleField(f);
}
bool SymmetryAdapted_IBS::CanAssembleProjectors(const SpeciesProjectorSet& ps) const
{
    const auto* pp=dynamic_cast<const Orbital_PP_IBS<double>*>(itsRaw);
    return pp && pp->CanAssembleProjectors(ps);
}
rsmat_t SymmetryAdapted_IBS::MakeSpeciesFieldMatrix(const Structure* st, const SpeciesRadialField& f, FieldRange r) const
{
    return Transform(dynamic_cast<const Orbital_PP_IBS<double>&>(*itsRaw).MakeSpeciesFieldMatrix(st,f,r));
}
rsmat_t SymmetryAdapted_IBS::MakeProjectorMatrix(const Structure* st, const SpeciesProjectorSet& ps) const
{
    return Transform(dynamic_cast<const Orbital_PP_IBS<double>&>(*itsRaw).MakeProjectorMatrix(st,ps));
}

// --- d/dR (GeometryOptimizationPlan §13.2; lane B 2026-11-09) --------------------------------------------------
// Every pullback is the RAW basis's on this block's density folded up to AO (see the header): Tr(D ∂(Oᵀ M O)) = Tr(O D Oᵀ ∂M), O fixed.
rsmat_t SymmetryAdapted_IBS::FoldUp(const rsmat_t& D) const
{
    rsmat_t Dao=blazem::zero<double>(AODimension());
    AddAODensity(Dao, D);
    return Dao;
}
namespace
{
const Orbital_PP_IBS<double>* RawPP(const Orbital_1E_IBS<double>* raw) {return dynamic_cast<const Orbital_PP_IBS<double>*>(raw);}  // abstract -> abstract
Outcome<rvec3vec_t,AD::NoGradient> NoRaw(const char* face, const char* what)
{return NoGeometryDerivativeYet<rvec3vec_t>(std::string("a symmetry-adapted (SALC) block whose raw basis has no ")+face, what);}
}

bool SymmetryAdapted_IBS::CanPullbackOverlap() const {return itsRaw->CanPullbackOverlap();}
bool SymmetryAdapted_IBS::CanPullbackKinetic() const {return itsRaw->CanPullbackKinetic();}
bool SymmetryAdapted_IBS::CanPullbackNuclear() const {return itsRaw->CanPullbackNuclear();}
Outcome<rvec3vec_t,AD::NoGradient> SymmetryAdapted_IBS::MakeOverlapPullback(const hmat_t<double>& D, const Structure* st) const
{return itsRaw->MakeOverlapPullback(FoldUp(D), st);}
Outcome<rvec3vec_t,AD::NoGradient> SymmetryAdapted_IBS::MakeKineticPullback(const hmat_t<double>& D, const Structure* st) const
{return itsRaw->MakeKineticPullback(FoldUp(D), st);}
bool SymmetryAdapted_IBS::CanVirialOverlap() const {return itsRaw->CanVirialOverlap();}
bool SymmetryAdapted_IBS::CanVirialKinetic() const {return itsRaw->CanVirialKinetic();}
Outcome<rmat3d_t,AD::NoGradient> SymmetryAdapted_IBS::MakeOverlapVirial(const hmat_t<double>& D, const Structure* st) const
{return itsRaw->MakeOverlapVirial(FoldUp(D), st);}
Outcome<rmat3d_t,AD::NoGradient> SymmetryAdapted_IBS::MakeKineticVirial(const hmat_t<double>& D, const Structure* st) const
{return itsRaw->MakeKineticVirial(FoldUp(D), st);}
Outcome<rvec3vec_t,AD::NoGradient> SymmetryAdapted_IBS::MakeNuclearPullback(const hmat_t<double>& D, const Structure* st) const
{return itsRaw->MakeNuclearPullback(FoldUp(D), st);}

bool SymmetryAdapted_IBS::CanPullbackSpeciesField(const SpeciesRadialField& f) const
{
    const auto* pp=RawPP(itsRaw);
    return pp && pp->CanPullbackSpeciesField(f);
}
bool SymmetryAdapted_IBS::CanPullbackProjectors(const SpeciesProjectorSet& ps) const
{
    const auto* pp=RawPP(itsRaw);
    return pp && pp->CanPullbackProjectors(ps);
}
Outcome<rvec3vec_t,AD::NoGradient> SymmetryAdapted_IBS::MakeSpeciesFieldPullback(const Structure* st, const SpeciesRadialField& f, FieldRange r,
                                                                                const hmat_t<double>& D) const
{
    const auto* pp=RawPP(itsRaw);
    if (!pp) return NoRaw("species-field face", "species-field pullback");
    return pp->MakeSpeciesFieldPullback(st, f, r, FoldUp(D));
}
Outcome<rvec3vec_t,AD::NoGradient> SymmetryAdapted_IBS::MakeProjectorPullback(const Structure* st, const SpeciesProjectorSet& ps, const hmat_t<double>& D) const
{
    const auto* pp=RawPP(itsRaw);
    if (!pp) return NoRaw("species-field face", "projector pullback");
    return pp->MakeProjectorPullback(st, ps, FoldUp(D));
}

bool SymmetryAdapted_IBS::CanPullbackRepulsion3C(const rFIT_CD_ABS& c) const {return itsRawDFT && itsRawDFT->CanPullbackRepulsion3C(c);}
Outcome<rvec3vec_t,AD::NoGradient> SymmetryAdapted_IBS::MakeRepulsion3CPullback(const rFIT_CD_ABS& c, const hmat_t<double>& D, const rvec_t& a,
                                                                               const Structure* st) const
{
    if (!itsRawDFT) return NoRaw("DFT face", "3-centre repulsion pullback");
    return itsRawDFT->MakeRepulsion3CPullback(c, FoldUp(D), a, st);
}

// The partner must fold its own block up (the whole-system face) into the SAME AO space; the raw basis must pull back against itself.
bool SymmetryAdapted_IBS::CanPullbackERI4(const Orbital_HF_IBS<double>* other) const
{
    const auto* ws=dynamic_cast<const WholeSystemFock_IBS<double>*>(other);
    return ws && ws->AODimension()==AODimension() && itsRawHF && itsRawHF->CanPullbackERI4(itsRawHF);
}
static Outcome<rvec3vec_t,AD::NoGradient> PairERI4(const SymmetryAdapted_IBS& self, const Orbital_HF_IBS<double>* rawHF, const rsmat_t& Dao_a,
                                                    const rsmat_t& Db, const Orbital_HF_IBS<double>* other, const Structure* st, bool exchange)
{
    const char* what = exchange ? "exchange pullback" : "Coulomb pullback";
    if (!self.CanPullbackERI4(other))
        return NoGeometryDerivativeYet<rvec3vec_t>("a symmetry-adapted (SALC) block against a partner it cannot fold into its AO space, or "
                                                   "over a raw basis with no ERI pullback", what);
    rsmat_t Dao_b=blazem::zero<double>(self.AODimension());
    dynamic_cast<const WholeSystemFock_IBS<double>&>(*other).AddAODensity(Dao_b, Db);   // the partner uses its OWN columns
    return exchange ? rawHF->MakeExchangePullback(Dao_a, Dao_b, rawHF, st) : rawHF->MakeDirectPullback(Dao_a, Dao_b, rawHF, st);
}
Outcome<rvec3vec_t,AD::NoGradient> SymmetryAdapted_IBS::MakeDirectPullback(const rsmat_t& Da, const rsmat_t& Db, const Orbital_HF_IBS<double>* other,
                                                                          const Structure* st) const
{return PairERI4(*this, itsRawHF, FoldUp(Da), Db, other, st, false);}
Outcome<rvec3vec_t,AD::NoGradient> SymmetryAdapted_IBS::MakeExchangePullback(const rsmat_t& Da, const rsmat_t& Db, const Orbital_HF_IBS<double>* other,
                                                                            const Structure* st) const
{return PairERI4(*this, itsRawHF, FoldUp(Da), Db, other, st, true);}

// The whole-system twin of MakeAOFock: ONE raw pullback with the assembled AO density on both sides.
Outcome<rvec3vec_t,AD::NoGradient> SymmetryAdapted_IBS::MakeAOPullback(const rsmat_t& Dao, bool exchange, const Structure* st) const
{
    const char* what = exchange ? "exchange pullback" : "Coulomb pullback";
    if (!itsRawHF || !itsRawHF->CanPullbackERI4(itsRawHF)) return NoRaw("ERI pullback", what);
    if (!st) throw std::invalid_argument("SymmetryAdapted_IBS::MakeAOPullback: a geometry pullback needs the structure");
    if (Dao.rows()!=AODimension()) throw std::logic_error("SymmetryAdapted_IBS::MakeAOPullback: the AO density is not AO-sized");
    if (blazem::max(blazem::abs(Dao))==0.0) return Outcome<rvec3vec_t,AD::NoGradient>::Ok(rvec3vec_t(st->GetNumAtoms(), rvec3_t(0,0,0)));   // bilinear: an empty density gives zero
    return exchange ? itsRawHF->MakeExchangePullback(Dao, Dao, itsRawHF, st) : itsRawHF->MakeDirectPullback(Dao, Dao, itsRawHF, st);
}

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
