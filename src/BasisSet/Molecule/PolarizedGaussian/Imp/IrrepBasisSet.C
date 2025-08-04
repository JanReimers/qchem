// File: BasisSet.C  Polarized Gaussian basis set, for MO calculations.
module;
#include <cassert>
#include <algorithm> //Need std::max
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <memory>
#include <vector>


module qchem.BasisSet.Molecule.PolarizedGaussian;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.GaussianRF;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.IntegralEngine;
import qchem.BasisSet.Molecule.PolarizedGaussian.Internal.Readers.Gaussian94;
import qchem.BasisSet;
import qchem.Cluster;
import qchem.Atom;
import qchem.Symmetry.Unit;
import qchem.stl_io;

namespace PolarizedGaussian
{

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
IrrepBasisSet::IrrepBasisSet(Reader* bsr, const Cluster* cl)
    : IrrepBasisSet_Common<double>(new UnitQN)
{
    //
    //  Read in all the radial functions.  These are usually contracted Gaussians, but could also
    //  be single Gaussians.
    //
    std::vector<RadialFunction*> radials;
    std::vector<std::vector<int> >    Ls;
    for (auto& atom:*cl) //Loop over atoms.
    {
        bsr->FindAtom(*atom);
        RadialFunction* rf=0;
        while ((rf=bsr->ReadNext(*atom))) //Read in the radial function/
        {
            bool duplicate=false;
            std::vector<RadialFunction*>::iterator b(radials.begin());
            for (size_t i=0; b!=radials.end(); i++,b++)
                if (**b==*rf) //Check for a duplicate, ingnoring Lmax.
                {
                    duplicate=true;
                    std::vector<int> newLs=bsr->GetLs();
                    //std::vector<int>& Li=Ls[i];
                    bool UseNewRF=Max(newLs) > Max(Ls[i]);
                    for (auto l:newLs)
                        if (std::find(Ls[i].begin(),Ls[i].end(),l)!=Ls[i].end()) Ls[i].push_back(l); //Add elements not in common.                        

                    if (UseNewRF)
                    {
                        radials.erase(b);
                        radials.insert(b,rf);
                    }
                    else
                    {
                        delete rf;
                    }
                }
            if(!duplicate)
            {
                radials.push_back(rf);
                Ls     .push_back(bsr->GetLs());
            }
        }
    }

//
//  Automatically build the basis set from a list of atoms and a basis function reader.
//
    int nbasis=1,i=0;
    for (auto r:radials)
    {
        Block* bfb=new Block(r,nbasis);
        for (auto& p:MakePolarizations(Ls[i]))
        {
            bfb->Add(p);
            nbasis++;
        }
        itsBlocks.push_back(std::unique_ptr<Block>(bfb));
        i++;
    }
    
    std::vector<const Block*> bls;
    for (auto& bl:itsBlocks) bls.push_back(bl.get());
    IrrepIEClient::Init(bls);
};
IrrepBasisSet::IrrepBasisSet(const Vector<double>& es, size_t LMax, const Cluster* cl)
    : IrrepBasisSet_Common<double>(new UnitQN)
{
    int nbasis=1;
    for (auto& atom:*cl)
    {
        for (size_t L=0;L<=LMax;L++)
        {
            std::vector<Polarization> Ps=MakePolarizations(L);
            for (auto e:es)
            {
                RadialFunction* r=new GaussianRF(e,atom->itsR,L);
                Block* bfb=new Block(r,nbasis);
                for (auto& p:Ps)
                {
                    bfb->Add(p);
                    nbasis++;
                }
                itsBlocks.push_back(std::unique_ptr<Block>(bfb));
            }
            
        }
    }
    std::vector<const Block*> bls;
    for (auto& bl:itsBlocks) bls.push_back(bl.get());
    IrrepIEClient::Init(bls);
    
}
// Single atom version
IrrepBasisSet::IrrepBasisSet(const Vector<double>& es, size_t L)
    : IrrepBasisSet_Common<double>(new UnitQN())
{
    int nbasis=1;
    std::vector<Polarization> Ps=MakePolarizations(L);
    for (auto e:es)
    {
        RadialFunction* r=new GaussianRF(e,RVec3(0,0,0),L);
        Block* bfb=new Block(r,nbasis);
        for (auto& p:Ps)
        {
            bfb->Add(p);
            nbasis++;
        }
        itsBlocks.push_back(std::unique_ptr<Block>(bfb));
    }
    std::vector<const Block*> bls;
    for (auto& bl:itsBlocks) bls.push_back(bl.get());
    IrrepIEClient::Init(bls);
}
//----------------------------------------------------------------
//
//  This contructor is used by Clone(RVec); only.
//
IrrepBasisSet::IrrepBasisSet(const IrrepBasisSet* bs, const bv_t& theBlocks)
    : IrrepBasisSet_Common<double>(*bs)
    // , Orbital_IBS_Common<double>(bs->itsLAParams,theDB)
    // , IE_Common(db)
    //, itsBlocks(theBlocks) //Need to clone all the blocks.
{
    assert(false);
    //Need to clone all the blocks.
    // No UT coverage
    //MakeBasisFunctions(); //Compiler says these calls are ambiguous.  BUG
//    TBasisSetImplementation<double>::Insert(bs->GetIntegralEngine()->Clone());
}

IrrepBasisSet::Vec     IrrepBasisSet::operator() (const RVec3& r) const
{
    Vec ret(size());
    for (size_t i=0;i<size();i++)
    {
        const RadialFunction& rf=*radials[i];
        RVec3 dr=r-rf.GetCenter();
        ret(i+1)= ns[i]*pols[i](dr) * rf(r);
    }
    return ret;
}
IrrepBasisSet::Vec3Vec IrrepBasisSet::Gradient   (const RVec3& r) const
{
    Vec3Vec ret(size());
    for (size_t i=0;i<size();i++)
    {
        const RadialFunction& rf=*radials[i];
        RVec3 dr=r-rf.GetCenter();
        ret(i+1)= ns[i]*(pols[i].Gradient(dr) * rf(r) + pols[i](dr) * rf.Gradient(r));
    }
    return ret;
}

std::ostream& IrrepBasisSet::Write(std::ostream& os) const
{
    // No UT coverage
    return os << itsBlocks;
}

//----------------------------------------------------------------
//
// Orbital PG basis set.
//
Orbital_IBS::Orbital_IBS(const db_t* db, Reader* bsr, const Cluster* cl)
    : IrrepBasisSet(bsr,cl)
    , Orbital_IBS_Common<double>()
    , Orbital_HF_IBS_Common<double>(db)
    , Orbital_IE(db)
{};
Orbital_IBS::Orbital_IBS(const db_t* db, const Vector<double>& exponents, size_t L, const Cluster* cl)
    : IrrepBasisSet(exponents,L,cl)
    , Orbital_IBS_Common<double>()
    , Orbital_HF_IBS_Common<double>(db)
    , Orbital_IE(db)
{};
Orbital_IBS::Orbital_IBS(const db_t* db, const Vector<double>& exponents, size_t L)
    : IrrepBasisSet(exponents,L)
    , Orbital_IBS_Common<double>()
    , Orbital_HF_IBS_Common<double>(db)
    , Orbital_IE(db)
{};
    
::Fit_IBS* Orbital_IBS::CreateCDFitBasisSet(const ::BasisSet* bs,const Cluster* cl) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    // The A1 files support Z=1-54 (H-Te)  A2 version only go up to Zn
    PolarizedGaussian::Gaussian94Reader reader("../../../BasisSetData/A1_coul.bsd");
    return new Fit_IBS(db,&reader,cl);
}
::Fit_IBS* Orbital_IBS::CreateVxcFitBasisSet(const ::BasisSet* bs,const Cluster* cl) const
{
    auto db=dynamic_cast<const DB_cache<double>*>(bs);
    // The A1 files support Z=1-54 (H-Te)  A2 version only go up to Zn
    PolarizedGaussian::Gaussian94Reader reader("../../../BasisSetData/A1_exch.bsd");
    return new Fit_IBS(db,&reader,cl);
}

//----------------------------------------------------------------
//
//  Fit PG basis set.
//
Fit_IBS::Fit_IBS(const DB_cache<double>* db , Reader* bsr, const Cluster* cl)
: IrrepBasisSet(bsr,cl)
, Fit_IE(db)
{};

 
} //namespace PolarizedGaussian
