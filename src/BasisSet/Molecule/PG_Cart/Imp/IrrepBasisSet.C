// File: BasisSet.C  Polarized Gaussian basis set, for MO calculations.
module;
#include <cassert>
#include <algorithm> //Need std::max
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <vector>

module qchem.BasisSet.Molecule.PG_Cart;
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.GaussianRF;
import qchem.BasisSet.Molecule.Readers.Gaussian94;
import qchem.BasisSet.Molecule.BasisFiles;   // the auto fit-basis files (path owned by BasisFiles)
import qchem.BasisSet;
import qchem.Structure;
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
//  This contructor is used by Clone(RVec); only.
//
IrrepBasisSet::IrrepBasisSet(const IrrepBasisSet* bs, const bv_t& theBlocks)
    : IrrepBasisSetImp<double>(*bs)
{
    assert(false);
    // No UT coverage
}

IrrepBasisSet::~IrrepBasisSet() {};


rvec_t IrrepBasisSet::operator() (const rvec3_t& r) const
{
    rvec_t ret(size());
    for (size_t i=0;i<size();i++)
    {
        const GaussianRF& rf=*radials[i];
        rvec3_t dr=r-rf.GetCenter();
        ret[i]= ns[i]*pols[i](dr) * rf(r);
    }
    return ret;
}
rvec3vec_t IrrepBasisSet::Gradient   (const rvec3_t& r) const
{
    rvec3vec_t ret(size());
    for (size_t i=0;i<size();i++)
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
    
FIT_CD_ABS* Orbital_IBS::CreateCDFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    // The A1 files support Z=1-54 (H-Te)  A2 version only go up to Zn
    Gaussian94Reader reader(BasisFile("A1_coul.bsd"));
    auto* f = new EFit_IBS(&reader,cl);
    f->SetMesh(*cl, mp);
    return f;
}
FIT_SF_ABS* Orbital_IBS::CreateVxcFitBasisSet(const Structure* cl, const qcMesh::MeshParams& mp) const
{
    // The A1 files support Z=1-54 (H-Te)  A2 version only go up to Zn
    Gaussian94Reader reader(BasisFile("A1_exch.bsd"));
    auto* f = new EFit_IBS(&reader,cl);
    f->SetMesh(*cl, mp);
    return f;
}
// Orbital_1E_IBS_ABS::GetAoShells: this orbital IBS IS-A PGData, so it hands its own Cartesian data to the extractor.
std::vector<Symmetry::AoShell> Orbital_IBS::GetAoShells() const {return ExtractAoShells(*this);}

//----------------------------------------------------------------
//
//  Fit PG basis set.
//
EFit_IBS::EFit_IBS(Reader* bsr, const Structure* cl)
: IrrepBasisSet(bsr,cl)
{};

 
} //namespace qchem::BasisSet::Molecule::PG_Cart
