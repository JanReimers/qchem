// File: BasisSet/Molecule/PG_LibCint/Imp/IrrepBasisSet.C  Cartesian PG basis, integrated by libcint.
//
// The radial-read / de-dup loop and the Cartesian angular expansion (MakePolarizations) are copied verbatim
// from PG_Cart so the (radial, polarization) component set -- and therefore the component ORDER -- is
// identical; only the integral ENGINE differs.  After PGData::Init flattens the components, NR_Evaluator::
// Init(cl) packs them into libcint's atm/bas/env.
module;
#include <cassert>
#include <algorithm> //std::max, std::find
#include <string>
#include <memory>
#include <vector>
#include <stdexcept>

module qchem.BasisSet.Molecule.PG_LibCint;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Polarization;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.Internal.Block;
import qchem.BasisSet.Molecule.Reader;
import qchem.BasisSet.Molecule.PG_Cart.Symmetry;   // ExtractAoShells(const PGData&) -- libcint-Cartesian reuses it
import qchem.Structure;
import qchem.Symmetry.Unit;
import qchem.stl_io;
import qchem.Math;
import qchem.Blaze;

namespace qchem::BasisSet::Molecule::PG_LibCint
{
using Cart::GaussianRF;
using Cart::Polarization;
using Cart::Block;

template <class T> T Max(const std::vector<T>& v) {return *std::max_element(v.begin(), v.end());}

// Cartesian angular expansion -- verbatim from PG_Cart (m outer, l middle; n = L-m-l).
static std::vector<Polarization> MakePolarizations(const std::vector<int>& Ls)
{
    std::vector<Polarization> ret;
    for (int bl:Ls)
        for(int m=0; m<=bl; m++)
            for(int l=0; l<=bl-m; l++)
                ret.push_back(Polarization(bl-m-l,l,m));
    return ret;
}
static std::vector<Polarization> MakePolarizations(int L)
{
    std::vector<Polarization> ret;
    for(int m=0; m<=L; m++)
        for(int l=0; l<=L-m; l++)
            ret.push_back(Polarization(L-m-l,l,m));
    return ret;
}

//----------------------------------------------------------------
//  Orbital Cartesian-Gaussian basis set, read identically to PG_Cart.
//
IrrepBasisSet::IrrepBasisSet(Reader* bsr, const Structure* cl, bool spherical)
    : IrrepBasisSetImp<double>(sym_t(new UnitQN)), itsSpherical(spherical)
{
    std::vector<GaussianRF*> radials;
    std::vector<std::vector<int> >    Ls;
    for (auto atom:*cl)
    {
        bsr->FindAtom(*atom);
        GaussianRF* rf=0;
        while ((rf=bsr->ReadNext(*atom)))
        {
            bool duplicate=false;
            std::vector<GaussianRF*>::iterator b(radials.begin());
            for (size_t i=0; b!=radials.end(); i++,b++)
                if (**b==*rf) //Check for a duplicate, ignoring Lmax.
                {
                    duplicate=true;
                    std::vector<int> newLs=bsr->GetLs();
                    bool UseNewRF=Max(newLs) > Max(Ls[i]);
                    for (auto l:newLs)
                        if (std::find(Ls[i].begin(),Ls[i].end(),l)!=Ls[i].end()) Ls[i].push_back(l);
                    if (UseNewRF) { radials.erase(b); radials.insert(b,rf); }
                    else          { delete rf; }
                }
            if(!duplicate) { radials.push_back(rf); Ls.push_back(bsr->GetLs()); }
        }
    }

    int i=0;
    for (auto r:radials)
    {
        Block* bfb=new Block(r);
        for (auto& p:MakePolarizations(Ls[i])) bfb->Add(p);
        itsBlocks.push_back(std::unique_ptr<Block>(bfb));
        i++;
    }
    std::vector<const Block*> bls;
    for (auto& bl:itsBlocks) bls.push_back(bl.get());
    PGData::Init(bls);
    NR_Evaluator::Init(cl, spherical);   // pack the flattened components into libcint (cart or sph)
};

IrrepBasisSet::IrrepBasisSet(const rvec_t& es, size_t LMax, const Structure* cl, bool spherical)
    : IrrepBasisSetImp<double>(sym_t(new UnitQN)), itsSpherical(spherical)
{
    for (auto atom:*cl)
        for (size_t L=0;L<=LMax;L++)
        {
            std::vector<Polarization> Ps=MakePolarizations((int)L);
            for (auto e:es)
            {
                Block* bfb=new Block(new GaussianRF(e,atom->itsR,(int)L));
                for (auto& p:Ps) bfb->Add(p);
                itsBlocks.push_back(std::unique_ptr<Block>(bfb));
            }
        }
    std::vector<const Block*> bls;
    for (auto& bl:itsBlocks) bls.push_back(bl.get());
    PGData::Init(bls);
    NR_Evaluator::Init(cl, spherical);
}

IrrepBasisSet::~IrrepBasisSet() {};

// Cartesian grid-eval (used only on the DFT mesh, not HF).  In spherical mode the delivered functions are
// libcint's real solid harmonics in libcint's own order/convention, which we deliberately do NOT reproduce
// here (the spherical libcint path is an HF-only oracle, validated by the basis-ordering-invariant energy);
// so grid-eval is unavailable there.
rvec_t IrrepBasisSet::operator() (const rvec3_t& r) const
{
    assert(!itsSpherical && "grid-eval unavailable for the spherical libcint oracle (HF-only)");
    rvec_t ret(size());
    for (size_t i=0;i<size();i++)
    {
        const GaussianRF& rf=*radials[i];
        rvec3_t dr=r-rf.GetCenter();
        ret[i]= ns[i]*pols[i](dr) * rf(r);
    }
    return ret;
}
rvec3vec_t IrrepBasisSet::Gradient (const rvec3_t& r) const
{
    assert(!itsSpherical && "grid-eval unavailable for the spherical libcint oracle (HF-only)");
    rvec3vec_t ret(size());
    for (size_t i=0;i<size();i++)
    {
        const GaussianRF& rf=*radials[i];
        rvec3_t dr=r-rf.GetCenter();
        ret[i]= ns[i]*(pols[i].Gradient(dr) * rf(r) + pols[i](dr) * rf.Gradient(r));
    }
    return ret;
}

std::ostream& IrrepBasisSet::Write(std::ostream& os) const {return os << BasisSetID();}

//----------------------------------------------------------------
Orbital_IBS::Orbital_IBS(Reader* bsr, const Structure* cl, bool sph)                : IrrepBasisSet(bsr,cl,sph) {};
Orbital_IBS::Orbital_IBS(const rvec_t& es, size_t L, const Structure* cl, bool sph) : IrrepBasisSet(es,L,cl,sph) {};

// Orbital_1E_IBS::GetAoShells: libcint-Cartesian shares PG_Cart's Cartesian PGData layout, so it reuses that extractor.
// Spherical libcint carries libcint's own real-harmonic order/norm (S3b, not convention-matched) yet its
// PGData base still holds the Cartesian layout -- reading it as Cartesian is the silent trap, so throw.
std::vector<Symmetry::Molecule::AoShell> Orbital_IBS::GetAoShells() const
{
    if (IsSpherical())
        throw std::runtime_error("PG_LibCint::GetAoShells: libcint-spherical SALC is not wired (S3b) -- its "
                                 "real-harmonic convention is not matched; use engine=MnD for spherical symmetry.");
    return PG_Cart::ExtractAoShells(*this);   // *this IS-A PGData with Cartesian components
}

} //namespace qchem::BasisSet::Molecule::PG_LibCint
