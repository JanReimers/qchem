// File: BasisSet.C  Polarized Gaussian basis set, for MO calculations.
module;
#include <cassert>
#include <stdexcept>
#include <complex>   // dcmplx operators (phase*chi in BlochPointValues); same as GPW's Evaluator.C
#include <algorithm> //Need std::max
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <vector>
#include <functional>

module qchem.BasisSet.Gaussian.PG_Cart;
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.BasisSet.Gaussian.Point.Readers.Gaussian94;
import qchem.BasisSet.Gaussian.Point.BasisFiles;   // the auto fit-basis files (path owned by BasisFiles)
import qchem.BasisSet;
import qchem.Structure;
import qchem.UnitCell;   // UnitCell (CollocateDensity grid<->cell map)
import qchem.Symmetry.Unit;
import qchem.BasisSet.Gaussian.PG_Cart.Symmetry;   // ExtractAoShells(const PGData&) -- for GetAoShells()
import qchem.stl_io;
import qchem.Streamable;
import qchem.Blaze;
import qchem.LatticePoints;   // LatticePointsInSphere (V-STRUCTURE-ISP F2: images off the cell)
import qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.Polarization;
import qchem.BasisSet.Gaussian.Lattice.LibcintPullbackWorkaround;   // ⚠ TEMPORARY (#3747): the lattice S/T derivatives, the ONE MnD/libcint seam

namespace qchem::BasisSet::Gaussian::PG_Cart
{
using namespace ::qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD;  // Cartesian glue moved out to PG_Cart_MnD

template <class T> T Max(const std::vector<T>& v)
{
    return *std::max_element(v.begin(), v.end());
}

std::vector<Polarization> MakePolarizations(const std::vector<int>& Ls)
{
    std::vector<Polarization> ret;
    for (std::vector<int>::const_iterator bl(Ls.begin()); bl!=Ls.end(); bl++)
        for(int m=0; m<=*bl; m++)
            for(int l=0; l<=*bl-m; l++)
                ret.push_back(Polarization(*bl-m-l,l,m));

    return ret;
}
std::vector<Polarization> MakePolarizations(int L)
{
    std::vector<Polarization> ret;
        for(int m=0; m<=L; m++)
            for(int l=0; l<=L-m; l++)
                ret.push_back(Polarization(L-m-l,l,m));

    return ret;
}

//----------------------------------------------------------------
//
//  Common implementation for orbital and fit basis sets.
//
IrrepBasisSet::IrrepBasisSet(Reader* bsr, const Structure* cl)
    : IrrepBasisSetImp<double>(sym_t(new UnitQN))
{
    //
    //  Read in all the radial functions (contracted or single Gaussians) and their l-sets (ReadShells: a shared-exponent
    //  shell of another l is its own block, for orbital and fit bases alike).
    //
    std::vector<ReadShell> shells=ReadShells(*bsr, *cl);
    std::vector<std::unique_ptr<GaussianRF>> radials;
    std::vector<std::vector<int> >    Ls;
    for (auto& sh : shells) { radials.push_back(std::move(sh.rf)); Ls.push_back(sh.Ls); }

//
//  Automatically build the basis set from a list of atoms and a basis function reader.
//
    int i=0;
    for (auto& r:radials)
    {
        Block* bfb=new Block(r.release());   // hand ownership of the radial to the Block (Block deletes it)
        for (auto& p:MakePolarizations(Ls[i]))
        {
            bfb->Add(p);
        }
        itsBlocks.push_back(std::unique_ptr<Block>(bfb));
        i++;
    }
    
    std::vector<const Block*> bls;
    for (auto& bl:itsBlocks) bls.push_back(bl.get());
    PGData::Init(bls);
};

IrrepBasisSet::IrrepBasisSet(const rvec_t& es, size_t LMax, const Structure* cl)
    : IrrepBasisSetImp<double>(sym_t(new UnitQN))
{
    for (const Site& atom : cl->Sites())
    {
        for (size_t L=0;L<=LMax;L++)
        {
            std::vector<Polarization> Ps=MakePolarizations(L);
            for (auto e:es)
            {
                GaussianRF* r=new GaussianRF(e,atom.R,L);
                Block* bfb=new Block(r);
                for (auto& p:Ps)
                {
                    bfb->Add(p);
                }
                itsBlocks.push_back(std::unique_ptr<Block>(bfb));
            }
            
        }
    }
    std::vector<const Block*> bls;
    for (auto& bl:itsBlocks) bls.push_back(bl.get());
    PGData::Init(bls);
}
// Single atom version
IrrepBasisSet::IrrepBasisSet(const rvec_t& es, size_t L)
    : IrrepBasisSetImp<double>(sym_t(new UnitQN))
{
    std::vector<Polarization> Ps=MakePolarizations(L);
    for (auto e:es)
    {
        GaussianRF* r=new GaussianRF(e,rvec3_t(0,0,0),L);
        Block* bfb=new Block(r);
        for (auto& p:Ps)
        {
            bfb->Add(p);
        }
        itsBlocks.push_back(std::unique_ptr<Block>(bfb));
    }
    std::vector<const Block*> bls;
    for (auto& bl:itsBlocks) bls.push_back(bl.get());
    PGData::Init(bls);
}
//----------------------------------------------------------------
//
//  SUBSET / prune constructor: copy \a bs's symmetry state, but take a SUBSET of its shells (\a theBlocks).
//  Deep-copies (Clone) the passed blocks so the new basis owns them, then re-flattens via PGData::Init -- so
//  the new basis has fewer functions (a vetted basis; doc/GPWPlan1.md §4a).  The blocks are Cloned, not moved,
//  so \a theBlocks may alias \a bs's own itsBlocks.
//
IrrepBasisSet::IrrepBasisSet(const IrrepBasisSet* bs, const bv_t& theBlocks)
    : IrrepBasisSetImp<double>(*bs)
{
    for (const auto& bl : theBlocks) itsBlocks.push_back(std::unique_ptr<Block>(bl->Clone()));
    std::vector<const Block*> bls;
    for (auto& bl : itsBlocks) bls.push_back(bl.get());
    PGData::Init(bls);
}

IrrepBasisSet::~IrrepBasisSet() {};


rvec_t IrrepBasisSet::operator() (const rvec3_t& r) const
{
    rvec_t ret(size(), 0.0);   // ZERO-FILLED: a screened-out function is left at its exact value, 0
    // SHELL HOIST: the flattened (radial x polarization) layout stores a shell's components
    // CONSECUTIVELY over the SAME GaussianRF (PGData::Init), so the contracted radial -- the exp
    // calls, the whole cost of this loop -- is evaluated ONCE per shell instead of once per
    // component (6x on a d shell, 10x on an f).  Pure hoist: identical value, identical order.
    // This is THE pointwise basis sweep behind the atom-centred XC mesh's Phi tables (and every
    // molecular Becke-grid evaluation), so it is on the hot path of every mesh-quadrature run.
    //
    // MAGNITUDE SCREEN (2026-08-20): skip a function whose bound |chi_i| < 1e-10 at this distance
    // (PGData::Reaches -- one-time, geometry-fixed).  It matters because the PERIODIC caller
    // (GPW_Evaluator::Eval) runs this sweep ONCE PER LATTICE IMAGE per mesh point: the image list
    // reaches as far as the MOST diffuse function, so an alpha=36 Mn d shell was paying a contracted
    // exp() at 20+ bohr, where its value is ~e^-14000.  The test is a squared distance against a
    // cached radius; the skipped work is exp() and intpow().  This is the project's own discipline --
    // the magnitude screen is the only truncation, and 1e-10 sits far below the ~1e-4 quadrature error
    // of any mesh this feeds.
    // size() is VIRTUAL and lives across a .so boundary, so it is neither inlined nor hoisted by the
    // compiler -- as a loop CONDITION it cost 2.7% of run cycles (PGData::size + IrrepBasisSet::size +
    // the PLT thunk, measured 2026-08-20), comparable to the whole cart->spherical transform layer.
    // Hoisting is bit-identical: the basis cannot change size mid-sweep.
    const size_t      n    =size();
    const rvec_t&     reach=Reaches();
    const GaussianRF* last=nullptr;
    double            rad =0.0;
    for (size_t i=0;i<n;i++)
    {
        const GaussianRF& rf=*radials[i];
        rvec3_t dr=r-rf.GetCenter();
        if (dr.x*dr.x+dr.y*dr.y+dr.z*dr.z > reach[i]*reach[i]) continue;
        if (&rf!=last) { rad=rf(r); last=&rf; }
        ret[i]= ns[i]*pols[i](dr) * rad;
    }
    return ret;
}
rvec3vec_t IrrepBasisSet::Gradient   (const rvec3_t& r) const
{
    const size_t n=size();          // virtual + cross-.so: hoisted, same reason as operator() above
    rvec3vec_t ret(n);
    for (size_t i=0;i<n;i++)
    {
        const GaussianRF& rf=*radials[i];
        ret[i]= ns[i]*rf.Gradient(pols[i], r);   // Enzyme's (GaussianRF_AD.C)
    }
    return ret;
}

// THE JETS of this block's functions on a slab of points riding with atom `owner` (GeometryOptimizationPlan §13.2): φ from operator(), ∇φ
// from Gradient(r) (the user's ruling, 2026-11-08; Enzyme's since 2026-11-09, GaussianRF_AD.C), and each function's centre matched to an atom of `st`.  The jet itself
// (φ, (centre, -∇φ), (owner, +∇φ)) is PointJets' to assemble.
Outcome<size_t,AD::NoGradient> IrrepBasisSet::JetsOn(const rvec3vec_t& pts, size_t owner, const Structure* st, qchem::BasisSet::PointJets<double>& out) const
{
    if (!st) throw std::invalid_argument("PG_Cart::IrrepBasisSet::JetsOn: the jets need the structure (the atoms the functions sit on)");
    std::vector<rvec3_t> R; st->ForEachSite([&](int, const rvec3_t& r, bool){ R.push_back(r); });
    if (owner>=R.size()) throw std::out_of_range("PG_Cart::IrrepBasisSet::JetsOn: owner "+std::to_string(owner)+" of "+std::to_string(R.size())+" atoms");
    const size_t n=size();
    std::vector<size_t> centre(n);
    for (size_t i=0;i<n;i++)
    {
        const rvec3_t c=radials[i]->GetCenter();
        size_t best=R.size(); double d2=0.0;
        for (size_t A=0;A<R.size();A++) { const rvec3_t d=R[A]-c; const double e=d*d; if (best==R.size() || e<d2) {best=A; d2=e;} }
        if (best==R.size() || d2>1e-24) throw std::logic_error("PG_Cart::IrrepBasisSet::JetsOn: function "+std::to_string(i)+" sits on no atom of the structure");
        centre[i]=best;
    }
    auto t=Evaluatable_IBS<double>::FillJets(out, owner, pts.size(), centre);
    for (size_t q=0;q<pts.size();q++)
    {
        const rvec_t     v=(*this)(pts[q]);
        const rvec3vec_t g=Gradient(pts[q]);
        for (size_t i=0;i<n;i++) { t.phi(q,i)=v[i]; t.grad[0](q,i)=g[i].x; t.grad[1](q,i)=g[i].y; t.grad[2](q,i)=g[i].z; }
    }
    return Outcome<size_t,AD::NoGradient>::Ok(pts.size());
}

std::ostream& IrrepBasisSet::Write(std::ostream& os) const
{
    // No UT coverage
    return os << itsBlocks;
}

std::string IrrepBasisSet::Name() const 
{
    return "Pol. Gaussian ";
}

//----------------------------------------------------------------
//
// Orbital PG basis set.
//
Orbital_IBS::Orbital_IBS(Reader* bsr, const Structure* cl)
    : IrrepBasisSet(bsr,cl)
{};
Orbital_IBS::Orbital_IBS(   const rvec_t& exponents, size_t L, const Structure* cl)
    : IrrepBasisSet(exponents,L,cl)
{};
Orbital_IBS::Orbital_IBS(   const rvec_t& exponents, size_t L)
    : IrrepBasisSet(exponents,L)
{};
    
rFIT_CD_ABS* Orbital_IBS::CreateCDFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    // The A1 files support Z=1-54 (H-Te)  A2 version only go up to Zn
    Gaussian94Reader reader(BasisFile("A1_coul.bsd"));
    return new EFit_IBS(&reader,cl,mp,"A1_coul");
}
rFIT_SF_ABS* Orbital_IBS::CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    // The A1 files support Z=1-54 (H-Te)  A2 version only go up to Zn
    Gaussian94Reader reader(BasisFile("A1_exch.bsd"));
    return new EFit_IBS(&reader,cl,mp,"A1_exch");
}
// Orbital_1E_IBS::GetAoShells: this orbital IBS IS-A PGData, so it hands its own Cartesian data to the extractor.
std::vector<Symmetry::Molecule::AoShell> Orbital_IBS::GetAoShells() const {return ExtractAoShells(*this);}

// Gaussian::LatticeSum1E: the orbital IBS IS-A NR_Evaluator, which owns the radials/pols/ns and the shifted
// two-centre kernels; forward the lattice sums straight to it (the GPW periodic-1E seam; the offsets are
// enumerated internally per shell pair -- there is no cut in R).
// The Bloch orbitals on a point set: Phi(q,i) = chi^k_i(pts[q]) = Sum_R phase(R) chi_i(pts[q]-R).
// The image sum lives HERE, not in the periodic caller, for the two reasons LatticeSum1E documents --
// a TRANSFORMED view can then transform once per POINT instead of once per IMAGE, and the reach data
// that bounds the enumeration is owned on this side.
//
// The enumeration REPRODUCES the caller's former rule exactly (same eps-derived reach, same
// LatticePointsInSphere radius, same per-point bounding-sphere screen, same POINT-outer/offset-inner order), so
// the Cartesian result is BIT-IDENTICAL to the old GPW_Evaluator::Eval path.  That is deliberate: this
// increment moves WHERE the sum happens, and nothing else, so any number that moves afterwards belongs
// to the transform reordering in the spherical view and not to a changed image set.
void Orbital_IBS::BlochPointValues(const rvec3vec_t& pts, const cellphase_t& phase, const UnitCell& A,
                                   mat_t<dcmplx>& Phi) const
{
    const size_t n=GetNumFunctions(), npts=pts.size();
    Phi=mat_t<dcmplx>(npts,n,dcmplx(0.0));
    if (npts==0) return;

    // A single orbital reaches sqrt(-ln eps/alpha_min) and its centre sits within a cell span of any
    // evaluation point, so reach + span covers every image the magnitude screen keeps (the caller's own
    // note; the per-function screen inside operator() then prunes it sparse per point).
    const double maxReach=sqrt(-log(1e-10)/MinExponent());   // qchem.CMath, per the include convention
    const rvec3_t ctr=A.ToCartesian(rvec3_t(0.5,0.5,0.5));
    double cellRad=0.0;
    for (int cx=0;cx<=1;cx++) for (int cy=0;cy<=1;cy++) for (int cz=0;cz<=1;cz++)
        cellRad=std::max(cellRad, norm(A.ToCartesian(rvec3_t(double(cx),double(cy),double(cz)))-ctr));
    const double Rcut=std::max(2.0*maxReach+2.0*A.GetMaximumCellEdge(), maxReach+2.0*cellRad);
    const double rr  =cellRad+maxReach;

    std::vector<ivec3_t> offs=LatticePointsInSphere(A, Rcut);
    std::vector<rvec3_t> Rc(offs.size());
    cvec_t               ph(offs.size());
    for (size_t k=0;k<offs.size();k++)
    {
        Rc[k]=A.ToCartesian(rvec3_t(double(offs[k].x),double(offs[k].y),double(offs[k].z)));
        ph[k]=phase(offs[k]);
    }

    for (size_t q=0;q<npts;q++)
        for (size_t k=0;k<Rc.size();k++)
        {
            const rvec3_t d=pts[q]-Rc[k]-ctr;
            if (d.x*d.x+d.y*d.y+d.z*d.z > rr*rr) continue;   // image cannot reach this point
            const rvec_t chi=(*this)(pts[q]-Rc[k]);
            for (size_t i=0;i<n;i++) Phi(q,i)+=ph[k]*chi[i];
        }
}

chmat_t Orbital_IBS::MakeOverlap(const cellphase_t& phase, const UnitCell& A) const {return NR_Evaluator::MakeOverlap(phase,A);}
cvec_t  Orbital_IBS::MakeOverlap(const cellphase_t& phase, const UnitCell& A,
                                 const Gaussian::LatticeSum1E::GaussianFunction& g) const {return NR_Evaluator::MakeOverlap(phase,A,g);}
cvec_t  Orbital_IBS::MakeOverlap(const Gaussian::LatticeSum1E::GaussianFunction& g) const {return NR_Evaluator::MakeOverlap(g);}
dcmplx  Orbital_IBS::MakeOverlap(const cellphase_t& phase, const UnitCell& A, const Gaussian::LatticeSum1E::GaussianFunction& g,
                                 const Gaussian::LatticeSum1E::GaussianFunction& h) const {return NR_Evaluator::MakeOverlap(phase,A,g,h);}
double  Orbital_IBS::MakeOverlap(const Gaussian::LatticeSum1E::GaussianFunction& g,
                                 const Gaussian::LatticeSum1E::GaussianFunction& h) const {return NR_Evaluator::MakeOverlap(g,h);}
chmat_t Orbital_IBS::MakeKinetic(const cellphase_t& phase, const UnitCell& A) const {return NR_Evaluator::MakeKinetic(phase,A);}
chmat_t Orbital_IBS::MakeNuclear(const cellphase_t& phase, const UnitCell& A, const Structure* cl) const {return NR_Evaluator::MakeNuclear(phase,A,cl);}

// THE LATTICE PULLBACKS (GeometryOptimizationPlan §13.2; doc/Records/CellHalf_Design.md §5; lane B 2026-11-10).  ⚠ TEMPORARY: both go
// through the ONE MnD/libcint seam (LibcintPullbackWorkaround, the user's 2026-11-10 ruling) until Enzyme can differentiate the MnD kernel
// (EnzymeAD/Enzyme#3747); then these two bodies become Enzyme over MnD and the seam is deleted.  They hand it MnD's own screened image list.
Gaussian::LatticeSum1E::GeometryPullback Orbital_IBS::LatticePullback(int which, const chmat_t& Abar, const cellphase_t& phase,
                                                                      const UnitCell& A, const Structure* st) const
{
    const LibcintWorkaround::ImageEnumerator images=[&](size_t i, size_t j, const std::function<void(const ivec3_t&, const rvec3_t&)>& cb)
        { ForSeriesImages(which==0 ? SeriesKernel::Overlap : SeriesKernel::Kinetic, i, j, A, cb); };   // the energy series' own list
    return LibcintWorkaround::LatticeOneElectronPullback(which, static_cast<const PGData&>(*this), images, Abar, phase, st);
}
Gaussian::LatticeSum1E::GeometryPullback Orbital_IBS::MakeLatticeOverlapPullback(const chmat_t& Abar, const cellphase_t& phase,
                                                                                 const UnitCell& A, const Structure* st) const
{   return LatticePullback(0, Abar, phase, A, st); }
Gaussian::LatticeSum1E::GeometryPullback Orbital_IBS::MakeLatticeKineticPullback(const chmat_t& Abar, const cellphase_t& phase,
                                                                                 const UnitCell& A, const Structure* st) const
{   return LatticePullback(1, Abar, phase, A, st); }
chmat_t Orbital_IBS::MakeLocalGaussian(const cellphase_t& phase, const UnitCell& A, const Structure* cl,
                                       const std::function<Gaussian::LatticeSum1E::GaussianFunction(int)>& opForZ) const
{   return NR_Evaluator::MakeLocalGaussian(phase,A,cl,opForZ); }
chmat_t Orbital_IBS::MakeLocalGaussian(const Structure* cl,
                                       const std::function<Gaussian::LatticeSum1E::GaussianFunction(int)>& opForZ) const
{   return NR_Evaluator::MakeLocalGaussian(cl,opForZ); }
double  Orbital_IBS::MaxExponent() const {return NR_Evaluator::MaxExponent();}
double  Orbital_IBS::MinExponent() const {return NR_Evaluator::MinExponent();}
const GPWTolerances& Orbital_IBS::Tolerances() const {return NR_Evaluator::Tolerances();}
void    Orbital_IBS::BindTolerances(const GPWTolerances& t) {NR_Evaluator::BindTolerances(t);}
std::vector<rvec_t> Orbital_IBS::CollocateDensity(const chmat_t& D, const cellphase_t& phase, const UnitCell& A,
                                                  const std::vector<ivec3_t>& N_L,
                                                  const std::vector<double>& ecut_L,
                                                  const LatticeScreener& screener, double relFieldSharp) const
{   return NR_Evaluator::CollocateDensity(D,phase,A,N_L,ecut_L,screener,relFieldSharp); }
chmat_t Orbital_IBS::IntegratePotential(const std::vector<rvec_t>& V_L, const cellphase_t& phase, const UnitCell& A,
                                        const std::vector<ivec3_t>& N_L,
                                        const std::vector<double>& ecut_L,
                                        const LatticeScreener& screener, double relCutoffScale,
                                        const chmat_t* screenD, double fieldSharpness, double relFieldSharp,
                                        const std::vector<size_t>* pairLevels) const
{   return NR_Evaluator::IntegratePotential(V_L,phase,A,N_L,ecut_L,screener,relCutoffScale,screenD,fieldSharpness,relFieldSharp,pairLevels); }
rvec3vec_t Orbital_IBS::IntegratePotentialPullback(const std::vector<rvec_t>& V_L, const chmat_t& D, const cellphase_t& phase,
                                                   const UnitCell& A, const std::vector<ivec3_t>& N_L,
                                                   const std::vector<double>& ecut_L, const LatticeScreener& screener,
                                                   double relFieldSharp, const Structure* st) const
{   return NR_Evaluator::IntegratePotentialPullback(V_L,D,phase,A,N_L,ecut_L,screener,relFieldSharp,st); }
rmat3d_t Orbital_IBS::IntegratePotentialVirial(const std::vector<rvec_t>& V_L, const chmat_t& D, const cellphase_t& phase,
                                               const UnitCell& A, const std::vector<ivec3_t>& N_L,
                                               const std::vector<double>& ecut_L, const LatticeScreener& screener,
                                               double relFieldSharp) const
{   return NR_Evaluator::IntegratePotentialVirial(V_L,D,phase,A,N_L,ecut_L,screener,relFieldSharp); }
std::vector<cvec_t> Orbital_IBS::CollocateTransition(const mat_t<dcmplx>& dD, const cellphase_t& ketPhase, const rvec3_t& q,
                                                     const UnitCell& A, const std::vector<ivec3_t>& N_L,
                                                     const std::vector<double>& ecut_L, double relFieldSharp) const
{   return NR_Evaluator::CollocateTransition(dD,ketPhase,q,A,N_L,ecut_L,relFieldSharp); }
mat_t<dcmplx> Orbital_IBS::IntegrateTransition(const std::vector<cvec_t>& v_L, const cellphase_t& ketPhase, const rvec3_t& q,
                                               const UnitCell& A, const std::vector<ivec3_t>& N_L,
                                               const std::vector<double>& ecut_L, double relFieldSharp) const
{   return NR_Evaluator::IntegrateTransition(v_L,ketPhase,q,A,N_L,ecut_L,relFieldSharp); }
std::vector<size_t> Orbital_IBS::StaticFieldPairLevels(const std::vector<double>& ecut_L,
                                                       double beta, double lnEps) const
{   return NR_Evaluator::StaticFieldPairLevels(ecut_L,beta,lnEps); }
void Orbital_IBS::EmitLatticeSumReport(const UnitCell& A) const
{   NR_Evaluator::EmitLatticeSumReport(A); }
size_t Orbital_IBS::SetStreamSymmetryOps(const std::vector<Symmetry::Lattice_3D::DirectOp>& ops,
                                         const UnitCell& A, const rvec3_t& kFrac) const
{   return NR_Evaluator::SetStreamSymmetryOps(ops,A,kFrac); }
size_t Orbital_IBS::StreamFoldOrder() const
{   return NR_Evaluator::StreamFoldOrder(); }

//----------------------------------------------------------------
//
//  Fit PG basis set.
//
EFit_IBS::EFit_IBS(Reader* bsr, const Structure* cl, const qcMesh::MeshParams& mp, const std::string& name)
: Fit_IBS(*cl,mp), IrrepBasisSet(bsr,cl), itsName(name)
{};

 
} //namespace qchem::BasisSet::Gaussian::PG_Cart
