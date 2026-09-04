// File: BasisSet.C  Polarized Gaussian basis set, for MO calculations.
module;
#include <cassert>
#include <complex>   // dcmplx operators (phase*chi in BlochPointValues); same as GPW's Evaluator.C
#include <algorithm> //Need std::max
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <vector>
#include <functional>

module qchem.BasisSet.Molecule.PG_Cart;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.BasisSet.Molecule.Readers.Gaussian94;
import qchem.BasisSet.Molecule.BasisFiles;   // the auto fit-basis files (path owned by BasisFiles)
import qchem.BasisSet;
import qchem.Structure;
import qchem.UnitCell;   // UnitCell (CollocateDensity grid<->cell map)
import qchem.Symmetry.Unit;
import qchem.BasisSet.Molecule.PG_Cart.Symmetry;   // ExtractAoShells(const PGData&) -- for GetAoShells()
import qchem.stl_io;
import qchem.Streamable;
import qchem.Math;
import qchem.Blaze;

namespace qchem::BasisSet::Molecule::PG_Cart
{
using namespace ::qchem::BasisSet::Molecule::Evaluators::PG_Cart_MnD;  // Cartesian glue moved out to PG_Cart_MnD

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
    //  Read in all the radial functions.  These are usually contracted Gaussians, but could also
    //  be single Gaussians.
    //
    std::vector<std::unique_ptr<GaussianRF>> radials;   // owns the radials: no manual delete, leak-proof
    std::vector<std::vector<int> >    Ls;
    for (auto atom:*cl) //Loop over atoms.
    {
        bsr->FindAtom(*atom);
        while (std::unique_ptr<GaussianRF> rf{bsr->ReadNext(*atom)}) //Read in the radial function.
        {
            bool duplicate=false;
            for (size_t i=0; i<radials.size(); i++)
                if (*radials[i]==*rf) //Check for a duplicate, ignoring Lmax.
                {
                    duplicate=true;
                    std::vector<int> newLs=bsr->GetLs();
                    bool UseNewRF=Max(newLs) > Max(Ls[i]);
                    for (auto l:newLs)
                        if (std::find(Ls[i].begin(),Ls[i].end(),l)!=Ls[i].end()) Ls[i].push_back(l); //Add elements not in common.
                        // KNOWN BUG (flagged, NOT fixed -- needs a careful pass): this condition is INVERTED (should be
                        // `==`).  A genuinely-new L is never added, so a same-exponent multi-L shell (s+p sharing one
                        // even-tempered window, or a shared-exponent density-fit shell) silently drops its higher L.
                        // Naively flipping it to `==` is correct for the ORBITAL basis (F loads 8s+8p=32 not 8) but
                        // it changes the DENSITY-FIT bases and DEGRADES Si2 pseudo-atom size-consistency 0.0026->0.024
                        // Ha (the isolated-atom and molecular fits shift by DIFFERENT amounts) -> the merged-shell fit
                        // path is non-additive and must be understood before flipping.  The valence generator
                        // sidesteps this by emitting DISJOINT per-l exponent windows (qchem.ValenceBasisGen).
                    if (UseNewRF) radials[i]=std::move(rf);  // replace lower-Lmax dup: old freed, new owned (no erase/insert/delete)
                    // else: rf is the lower-Lmax duplicate -> freed automatically at scope exit
                    break;                                   // a radial appears at most once in radials (deduped)
                }
            if(!duplicate)
            {
                radials.push_back(std::move(rf));
                Ls     .push_back(bsr->GetLs());
            }
        }
    }

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
    for (auto atom:*cl)
    {
        for (size_t L=0;L<=LMax;L++)
        {
            std::vector<Polarization> Ps=MakePolarizations(L);
            for (auto e:es)
            {
                GaussianRF* r=new GaussianRF(e,atom->itsR,L);
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
        rvec3_t dr=r-rf.GetCenter();
        ret[i]= ns[i]*(pols[i].Gradient(dr) * rf(r) + pols[i](dr) * rf.Gradient(r));
    }
    return ret;
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
    return new EFit_IBS(&reader,cl,mp);
}
rFIT_SF_ABS* Orbital_IBS::CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    // The A1 files support Z=1-54 (H-Te)  A2 version only go up to Zn
    Gaussian94Reader reader(BasisFile("A1_exch.bsd"));
    return new EFit_IBS(&reader,cl,mp);
}
// Orbital_1E_IBS::GetAoShells: this orbital IBS IS-A PGData, so it hands its own Cartesian data to the extractor.
std::vector<Symmetry::Molecule::AoShell> Orbital_IBS::GetAoShells() const {return ExtractAoShells(*this);}

// Molecule::LatticeSum1E: the orbital IBS IS-A NR_Evaluator, which owns the radials/pols/ns and the shifted
// two-centre kernels; forward the lattice sums straight to it (the GPW periodic-1E seam; the offsets are
// enumerated internally per shell pair -- there is no cut in R).
// The Bloch orbitals on a point set: Phi(q,i) = chi^k_i(pts[q]) = Sum_R phase(R) chi_i(pts[q]-R).
// The image sum lives HERE, not in the periodic caller, for the two reasons LatticeSum1E documents --
// a TRANSFORMED view can then transform once per POINT instead of once per IMAGE, and the reach data
// that bounds the enumeration is owned on this side.
//
// The enumeration REPRODUCES the caller's former rule exactly (same eps-derived reach, same
// CellsInSphere radius, same per-point bounding-sphere screen, same POINT-outer/offset-inner order), so
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

    std::vector<ivec3_t> offs=A.CellsInSphere(Rcut);
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
                                 const Molecule::LatticeSum1E::GaussianFunction& g) const {return NR_Evaluator::MakeOverlap(phase,A,g);}
cvec_t  Orbital_IBS::MakeOverlap(const Molecule::LatticeSum1E::GaussianFunction& g) const {return NR_Evaluator::MakeOverlap(g);}
chmat_t Orbital_IBS::MakeKinetic(const cellphase_t& phase, const UnitCell& A) const {return NR_Evaluator::MakeKinetic(phase,A);}
chmat_t Orbital_IBS::MakeNuclear(const cellphase_t& phase, const UnitCell& A, const Structure* cl) const {return NR_Evaluator::MakeNuclear(phase,A,cl);}
chmat_t Orbital_IBS::MakeLocalGaussian(const cellphase_t& phase, const UnitCell& A, const Structure* cl,
                                       const std::function<Molecule::LatticeSum1E::GaussianFunction(int)>& opForZ) const
{   return NR_Evaluator::MakeLocalGaussian(phase,A,cl,opForZ); }
chmat_t Orbital_IBS::MakeLocalGaussian(const Structure* cl,
                                       const std::function<Molecule::LatticeSum1E::GaussianFunction(int)>& opForZ) const
{   return NR_Evaluator::MakeLocalGaussian(cl,opForZ); }
double  Orbital_IBS::MaxExponent() const {return NR_Evaluator::MaxExponent();}
double  Orbital_IBS::MinExponent() const {return NR_Evaluator::MinExponent();}
double  Orbital_IBS::RelCutoffSafety() const {return NR_Evaluator::RelCutoffSafety();}
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
EFit_IBS::EFit_IBS(Reader* bsr, const Structure* cl, const qcMesh::MeshParams& mp)
: Fit_IBS(*cl,mp), IrrepBasisSet(bsr,cl)
{};

 
} //namespace qchem::BasisSet::Molecule::PG_Cart
