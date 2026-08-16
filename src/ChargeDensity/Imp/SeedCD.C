// File: ChargeDensity/Imp/SeedCD.C  Plane-wave SAD seed (G-space form-factor assembly).
module;
#include <cassert>
#include <cmath>
#include <complex>   // std::operator*(double, complex) -- the seed's k*rho-tilde
#include <map>
#include <memory>
#include <string>
#include <vector>
#include <cstddef>
#include <utility>

module qchem.ChargeDensity.SeedCD;
import qchem.ReciprocalLattice;        // ReciprocalLattice + UnitCell::MakeReciprocalCell (the seed's own Poisson metric)
import qchem.BasisSet.G_FieldEvaluator; // G_StructureFactor: the fit basis's analytic MakeFourierDensity
import qchem.Matrix3D;                 // Invert/Transpose (the periodic image window: (A^T A)^-1 diagonals)
import qchem.Parallel;                 // WorkerThreads (GPW_OMP_THREADS -- the batched seed sampling)

namespace qchem::ChargeDensity
{

// The seed is periodic (a plane-wave density), so its Structure IS a UnitCell -- derive the reciprocal
// lattice ONCE at construction (the Poisson metric B), rather than casting on every GetRepulsion3C.
// GetReciprocalLattice is the shared pry-out helper (qchem.ReciprocalLattice); it throws std::bad_cast
// when the Structure is not periodic, which is exactly this class's precondition.

SeedCD::SeedCD(std::shared_ptr<const BasisSet::cFIT_CD_ABS> fitBasis, const Structure* st,
                             const std::string& functional, const std::map<size_t,int>& ionicNvalByZ,
                             const Spin& channel)
    : itsFitBasis(fitBasis), itsStructure(st), itsRecip(GetReciprocalLattice(st)), itsCell(&GetUnitCell(st))
    , itsChannel(channel), itsCharge(0.0)
    , itsVersion(NextDensityVersion())   // shared global clock (no cross-kind collisions)
{
    assert(fitBasis);
    assert(st);
    const bool chan = IsPolarized(channel);   // Up/Down = one channel of the polarized seed
    const std::string db = "atomic_valence_densities.json";
    for (size_t i=0;i<st->GetNumAtoms();i++)
    {
        const Atom* a = (*st)[i];
        const size_t Z = a->itsZ;
        if (itsRadByZ.find(Z)==itsRadByZ.end())
        {
            // The NEUTRAL density fixes this species' neutral valence count; the IonicSAD target (Nval-q) is
            // ionicNvalByZ (default: neutral, i.e. no charge transfer).
            RadialDensity neutral = GetAtomicDensity((int)Z, functional, db);       // Nval<0 => neutral
            const int neutralNval = (int)std::lround(neutral.Charge());
            auto ti = ionicNvalByZ.find(Z);
            const int target = (ti!=ionicNvalByZ.end()) ? ti->second : neutralNval;

            std::shared_ptr<const RadialDensity> rad, radFlip;
            double scale;
            if (target<=0)                                                          // stripped cation -> nothing
                { rad=radFlip=std::make_shared<const RadialDensity>(std::move(neutral)); scale=0.0; }
            else if (chan && HasAtomicSpinPair((int)Z, functional, db, target))     // the Hund (majority, minority) pair
            {
                auto [maj,min] = GetAtomicSpinPair((int)Z, functional, db, target);
                auto majp = std::make_shared<const RadialDensity>(std::move(maj));
                auto minp = std::make_shared<const RadialDensity>(std::move(min));
                rad     = (channel==Spin::Up) ? majp : minp;                        // unflipped: up = majority
                radFlip = (channel==Spin::Up) ? minp : majp;                        // flipped: the pair swapped (-m site)
                scale   = 1.0;
            }
            else if (HasAtomicDensity((int)Z, functional, db, target))              // library has the charge state
                { rad=radFlip=std::make_shared<const RadialDensity>(GetAtomicDensity((int)Z,functional,db,target));
                  scale = chan ? 0.5 : 1.0; }                                       // pairless channel: rho/2
            else                                                                    // fall back: amplitude-scale neutral
                { scale = (chan ? 0.5 : 1.0)*double(target)/double(neutralNval);
                  rad=radFlip=std::make_shared<const RadialDensity>(std::move(neutral)); }

            itsRadByZ[Z]     = rad;
            itsRadFlipByZ[Z] = radFlip;
            itsScaleByZ[Z]   = scale;
        }
        const bool flip = chan && a->itsSpinFlip;                     // None: the total is flip-blind
        const auto& rc = flip ? itsRadFlipByZ.at(Z) : itsRadByZ.at(Z);
        itsCharge += itsScaleByZ.at(Z) * rc->Charge();                // this atom's (scaled) electron count
        itsRecentred.emplace_back(rc, a->itsR);                       // for real-space op(r), flip-aware
    }
    // Flip-partitioned sub-cells for the G assembly: the basis's MakeFourierDensity is SPECIES-keyed
    // (formFactor(Z, g2)), so flipped and unflipped sites of one species -- which carry DIFFERENT channel
    // radials -- must be summed in separate calls.  Only a channel seed with at least one flipped atom
    // needs the split; everything else keeps the single full-structure call.
    if (chan)
    {
        bool anyFlip=false;
        for (size_t i=0;i<st->GetNumAtoms();i++) if ((*st)[i]->itsSpinFlip) anyFlip=true;
        if (anyFlip)
        {
            const UnitCell& cell=GetUnitCell(st);   // checked precondition (periodic seed)
            itsGroupA=std::make_shared<UnitCell>(cell.GetCellMatrix());
            itsGroupB=std::make_shared<UnitCell>(cell.GetCellMatrix());
            for (size_t i=0;i<st->GetNumAtoms();i++)
            {
                const Atom* a=(*st)[i];
                (a->itsSpinFlip ? itsGroupB : itsGroupA)->Insert(new Atom(*a));
            }
        }
    }
    assert(itsCharge >= 0);
    assert(chan || itsCharge > 0);   // a lone CHANNEL may be legitimately empty (a fully-polarized minority)
}

// One flip-group's rho-tilde(G) = (1/Omega) Sum_{atoms in group} F(Z,|B.G|) e^{-i(B.G).R}: the seed's OWN
// density-fit basis (its grid engine) does the analytic structure-factor assembly over its {G}; we supply
// the per-species form factor (the group's radial 1-D Fourier transform).  Memoize the FT per (Z,g2) --
// many G share |G|.
ΔG_Map SeedCD::GroupDensity(const Structure* group,
                            const std::map<size_t,std::shared_ptr<const RadialDensity>>& radByZ) const
{
    if (group->GetNumAtoms()==0) return {};
    auto memo = std::make_shared<std::map<std::pair<int,double>,double>>();
    auto formFactor = [this,&radByZ,memo](int Z, double g2)->double
    {
        // The memo holds the RAW species form factor; the uniform ReScale (itsScale) and the per-species
        // multiplier (itsScaleByZ) are applied on the way out.
        double s = itsScale * itsScaleByZ.at(Z);
        auto key = std::make_pair(Z,g2);
        auto it = memo->find(key);
        if (it!=memo->end()) return s*it->second;
        double ff = radByZ.at(Z)->FormFactor(std::sqrt(g2));
        (*memo)[key]=ff;
        return s*ff;
    };
    auto* ge=dynamic_cast<const BasisSet::G_StructureFactor*>(itsFitBasis.get());
    assert(ge && "SeedCD's density-fit basis must be a G_StructureFactor (plane-wave grid engine)");
    return ge->MakeFourierDensity(group, formFactor);
}

// The whole seed's rho-tilde: one full-structure sweep for the spin-agnostic total (and for a channel with
// no flipped sites); a flip-staggered channel sums its two group sweeps (unflipped x this channel's radial
// + flipped x the swapped one).
ΔG_Map SeedCD::StructureFactorDensity() const
{
    if (!itsGroupA)
        return GroupDensity(itsStructure, itsRadByZ);
    ΔG_Map ro = GroupDensity(itsGroupA.get(), itsRadByZ);
    for (const auto& [dm,rt] : GroupDensity(itsGroupB.get(), itsRadFlipByZ)) ro[dm]+=rt;
    return ro;
}

// The seed's metric-free rho-tilde: no D to contract, so it IS the structure-factor density (the Vxc fit
// basis argument is ignored -- the seed's support is the orbital difference set, like every rho-tilde).
ΔG_Map SeedCD::GetFourierDensity(const BasisSet::cFIT_SF_ABS&) const
{
    return StructureFactorDensity();
}

// The seed's Coulomb projection V_H = 4pi rho-tilde/|G|^2: apply the diagonal Poisson kernel (reciprocal-
// LATTICE physics -- needs only B, not the basis's {G} set) to the structure-factor rho-tilde.  The seed
// owns its reciprocal lattice (itsRecip, built at construction), so there is NO basis round-trip.
ΔG_Map SeedCD::GetRepulsion3C(const BasisSet::cFIT_CD_ABS&) const
{
    ΔG_Map VH;
    for (const auto& [dm,rt] : StructureFactorDensity())
        if (double k=itsRecip.CoulombKernel(dm); k!=0.0) VH[dm]=k*rt;   // skip dm=0 (k==0)
    return VH;
}

// PERIODIC real-space evaluation (2026-08-11).  A plane-wave seed is cell-periodic, and its mesh consumers
// sample it at points WRAPPED into the home cell (MakePeriodicBeckeMesh stores kpt = r - A*n0; the uniform
// raster's cell-corner regions likewise) -- where the nearest atom is often a lattice IMAGE.  The image-less
// sum this used to be reads ~0 there: on MnO's AFM cell that erased the CORNER Mn's majority-density hump
// from the v_xc channel rasters and broke the sublattice mirror in the very first Fock build (the
// iteration-1 m_stag collapse -- doc/SymmetryUpgradePlan.md; gate GPW_SCF.MnOSeedVxcMirrorOnBeckeMesh).
// Each radial table has compact support [0, Rmax] (RadialDensity clamps to 0 beyond its last node), so the
// lattice-image sum is FINITE and EXACT: the per-axis image window is |f_j - n_j| <= Rmax*sqrt((A^T A)^-1_jj)
// (the interplanar-spacing relation MakePeriodicBeckeMesh's shell floor uses), and the only truncation is
// the table's own support -- no distance cut.
//
// The per-image visitor: every lattice image of atom \a i whose support ball reaches \a r.
template <class F> static void ForEachImage(const UnitCell& cell, const RecentredAtomicDensity& rc,
                                            const Matrix3D<double>& Minv, const rvec3_t& r, const F& f)
{
    const double  rmax=rc.Rmax(), r2max=rmax*rmax;
    const rvec3_t fr=cell.ToFractional(r-rc.Center());
    const rvec3_t b(rmax*std::sqrt(Minv(1,1)), rmax*std::sqrt(Minv(2,2)), rmax*std::sqrt(Minv(3,3)));
    for (int nx=int(std::ceil(fr.x-b.x)); nx<=int(std::floor(fr.x+b.x)); nx++)
        for (int ny=int(std::ceil(fr.y-b.y)); ny<=int(std::floor(fr.y+b.y)); ny++)
            for (int nz=int(std::ceil(fr.z-b.z)); nz<=int(std::floor(fr.z+b.z)); nz++)
            {
                const rvec3_t L=cell.ToCartesian(rvec3_t(nx,ny,nz));
                const rvec3_t d=r-rc.Center()-L;
                if (d.x*d.x+d.y*d.y+d.z*d.z<r2max) f(L);
            }
}

double SeedCD::operator()(const rvec3_t& r) const
{
    const Matrix3D<double>& A=itsCell->GetCellMatrix();
    const Matrix3D<double>  Minv=Invert(Transpose(A)*A);
    double rho=0;   // itsRecentred is parallel to the structure's atoms (flip-aware) -> per-atom scale by Z
    for (size_t i=0;i<itsRecentred.size();i++)
    {
        const double s=itsScaleByZ.at((*itsStructure)[i]->itsZ);
        if (s==0.0) continue;
        const RecentredAtomicDensity& rc=itsRecentred[i];
        ForEachImage(*itsCell, rc, Minv, r, [&](const rvec3_t& L){ rho += s*rc(r-L); });
    }
    return itsScale*rho;
}

// Batch: the ScalarFunction default loop, threaded (see the interface doc) -- and with the cell metric
// (A^T A)^-1 lifted OUT of the point loop, where op() recomputes a 3x3 inverse per point.
rvec_t SeedCD::operator()(const rvec3vec_t& rs) const
{
    const Matrix3D<double>& A=itsCell->GetCellMatrix();
    const Matrix3D<double>  Minv=Invert(Transpose(A)*A);
    rvec_t out(rs.size());
    auto at=[&](size_t q)
    {
        double rho=0;
        for (size_t i=0;i<itsRecentred.size();i++)
        {
            const double s=itsScaleByZ.at((*itsStructure)[i]->itsZ);
            if (s==0.0) continue;
            const RecentredAtomicDensity& rc=itsRecentred[i];
            ForEachImage(*itsCell, rc, Minv, rs[q], [&](const rvec3_t& L){ rho += s*rc(rs[q]-L); });
        }
        out[q]=itsScale*rho;
    };
#ifdef QCHEM_OPENMP
    if (const int nthreads=qchem::WorkerThreads(); nthreads>1)
    {
        #pragma omp parallel for schedule(static) num_threads(nthreads)
        for (size_t q=0;q<rs.size();q++) at(q);      // pure per-point work: no throw path, no reduction
        return out;
    }
#endif
    for (size_t q=0;q<rs.size();q++) at(q);
    return out;
}

rvec3_t SeedCD::Gradient(const rvec3_t& r) const
{
    const Matrix3D<double>& A=itsCell->GetCellMatrix();
    const Matrix3D<double>  Minv=Invert(Transpose(A)*A);
    rvec3_t g(0,0,0);
    for (size_t i=0;i<itsRecentred.size();i++)
    {
        const double s=itsScaleByZ.at((*itsStructure)[i]->itsZ);
        if (s==0.0) continue;
        const RecentredAtomicDensity& rc=itsRecentred[i];
        ForEachImage(*itsCell, rc, Minv, r, [&](const rvec3_t& L){ g += s*rc.Gradient(r-L); });
    }
    return itsScale*g;
}

void SeedCD::ReScale(double factor)
{
    itsScale *= factor;
    itsVersion = NextDensityVersion();
}

//------------------------------------------------------------------------------- PolarizedSeedCD

PolarizedSeedCD::PolarizedSeedCD(std::shared_ptr<const BasisSet::cFIT_CD_ABS> fitBasis, const Structure* st,
                                 const std::string& functional, const std::map<size_t,int>& ionicNvalByZ)
    : itsUp(fitBasis, st, functional, ionicNvalByZ, Spin::Up  )
    , itsDn(fitBasis, st, functional, ionicNvalByZ, Spin::Down)
    , itsVersion(NextDensityVersion())   // one serial for the pair (the channels mutate together)
{
    assert(GetTotalCharge()>0);
}

const cChargeDensity* PolarizedSeedCD::GetChannel(const Spin& s) const
{
    assert(IsPolarized(s) && "PolarizedSeedCD::GetChannel: ask for Up or Down, not the total");
    return s==Spin::Up ? static_cast<const cChargeDensity*>(&itsUp) : &itsDn;
}

ΔG_Map PolarizedSeedCD::GetFourierDensity(const BasisSet::cFIT_SF_ABS& c) const
{
    ΔG_Map ro = itsUp.GetFourierDensity(c);
    for (const auto& [dm,rt] : itsDn.GetFourierDensity(c)) ro[dm]+=rt;
    return ro;
}

ΔG_Map PolarizedSeedCD::GetRepulsion3C(const BasisSet::cFIT_CD_ABS& c) const
{
    ΔG_Map VH = itsUp.GetRepulsion3C(c);
    for (const auto& [dm,rt] : itsDn.GetRepulsion3C(c)) VH[dm]+=rt;
    return VH;
}

double  PolarizedSeedCD::operator()(const rvec3_t& r) const {return itsUp(r)+itsDn(r);}
rvec3_t PolarizedSeedCD::Gradient  (const rvec3_t& r) const {return itsUp.Gradient(r)+itsDn.Gradient(r);}

double PolarizedSeedCD::GetTotalCharge() const {return itsUp.GetTotalCharge()+itsDn.GetTotalCharge();}
double PolarizedSeedCD::GetTotalSpin  () const {return itsUp.GetTotalCharge()-itsDn.GetTotalCharge();}

void PolarizedSeedCD::ReScale(double factor)
{
    itsUp.ReScale(factor);
    itsDn.ReScale(factor);
    itsVersion = NextDensityVersion();
}

} //namespace
