// File: BasisSet.C  Polarized Gaussian basis set, for MO calculations.

#include "Common/stl_io.h"
#include <cassert>
#include <algorithm> //Need std::max
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "Symmetry/Unit.H"
#include <BasisSet/BasisSet.H>


#include "PolarizedGaussian/BasisFunction.H"
#include "PolarizedGaussian/IntegralEngine.H"
#include "PolarizedGaussian/Readers/Gaussian94.H"
#include "PolarizedGaussian/Radial/GaussianRF.H"

#include "PolarizedGaussian/IrrepBasisSet.H"

import qchem.Cluster;
import qchem.Atom;

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
    : IBS_Common(new UnitQN)
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
            for (index_t i=0; b!=radials.end(); i++,b++)
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
//
//  Now insert the basis functions.
//
    MakeBasisFunctions(ns); //ns from PolarizedGaussianIEClient
};
IrrepBasisSet::IrrepBasisSet(const Vector<double>& es, size_t LMax, const Cluster* cl)
    : IBS_Common(new UnitQN)
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
//
//  Now insert the basis functions.
//
    MakeBasisFunctions(ns); //ns from PolarizedGaussianIEClient
    
}
// Single atom version
IrrepBasisSet::IrrepBasisSet(const Vector<double>& es, size_t L)
    : IBS_Common(new UnitQN())
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
//
//  Now insert the basis functions.
//
    MakeBasisFunctions(ns); //ns from PolarizedGaussianIEClient
}
//----------------------------------------------------------------
//
//  This contructor is used by Clone(RVec); only.
//
IrrepBasisSet::IrrepBasisSet(const IrrepBasisSet* bs, const bv_t& theBlocks)
    : IBS_Common(*bs)
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


std::ostream& IrrepBasisSet::Write(std::ostream& os) const
{
    // No UT coverage
    IBS_Common::Write(os);
    return os << itsBlocks;
}

void IrrepBasisSet::MakeBasisFunctions(const RVec& norms)
{
    EmptyBasisFunctions();
    size_t i=1;
    for (auto& bl:itsBlocks)
        for (std::vector<Polarization>::const_iterator p(bl->itsPols.begin()); p!=bl->itsPols.end(); p++)
            IBS_Common::Insert(new BasisFunction(bl->itsRadial,*p,norms(i++)));
}

//----------------------------------------------------------------
//
// Orbital PG basis set.
//
Orbital_IBS::Orbital_IBS(const db_t* db, Reader* bsr, const Cluster* cl)
    : IrrepBasisSet(bsr,cl)
    , Orbital_IBS_Common<double>()
    , Orbital_IE(db)
{};
Orbital_IBS::Orbital_IBS(const db_t* db, const Vector<double>& exponents, size_t L, const Cluster* cl)
    : IrrepBasisSet(exponents,L,cl)
    , Orbital_IBS_Common<double>()
    , Orbital_IE(db)
{};
Orbital_IBS::Orbital_IBS(const db_t* db, const Vector<double>& exponents, size_t L)
    : IrrepBasisSet(exponents,L)
    , Orbital_IBS_Common<double>()
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
IrrepBasisSet* Orbital_IBS::Clone(const RVec3& newCenter) const
{
    // No UT coverage
    assert(false);
    return 0;
}

//----------------------------------------------------------------
//
//  Fit PG basis set.
//
Fit_IBS::Fit_IBS(const DB_cache<double>* db , Reader* bsr, const Cluster* cl)
: IrrepBasisSet(bsr,cl)
, TIBS_Common<double>()
, Fit_IE(db)
{};

::Fit_IBS*Fit_IBS::Clone(const RVec3&) const
{
// No UT coverage
    assert(false);
    return 0;
}
 
} //namespace PolarizedGaussian
