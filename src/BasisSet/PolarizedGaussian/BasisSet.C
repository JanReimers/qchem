// File: BasisSet.C  Polarized Gaussian basis set, for MO calculations.



#include "Imp/BasisSet/PolarizedGaussian/BasisFunction.H"
#include "Imp/BasisSet/PolarizedGaussian/BasisSet.H"
#include "Imp/BasisSet/PolarizedGaussian/IntegralEngine.H"
#include "Imp/BasisSet/PolarizedGaussian/Readers/RadialFunction.H"
#include "BasisSetImplementation/UnitSymmetryQN.H"
#include <Cluster.H>
#include "Imp/Containers/ptr_vector_io.h"
#include <cassert>
#include <algorithm> //Need std::max

namespace PolarizedGaussian
{

std::vector<Polarization> MakePolarizations(const std::vector<int>& Ls);
template <class T> T Max(const std::vector<T>& v)
{
    return *std::max_element(v.begin(), v.end());
}



//#######################################################################
//
//  Concrete  gaussian basis set.
//
BasisSet::BasisSet()
    : BasisSetImplementation()
    , TBasisSetImplementation<double>()
{};

BasisSet::
BasisSet(const LinearAlgebraParams& lap,IntegralDataBase<double>* theDB, Reader* bsr, const Cluster* cl)
    : BasisSetImplementation(new UnitSymmetryQN)
    , TBasisSetImplementation<double>(lap,theDB)
{
//
//  Read in all the radial functions.  These are usually contracted Gaussians, but could also
//  be single Gaussians.
//
    std::vector<RadialFunction*> radials;
    std::vector<std::vector<int> >    Ls;
    for (auto atom:*cl) //Loop over atoms.
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
        itsBlocks.push_back(bfb);
        i++;
    }
    
    IntegralEngine::blocks_t bls;
    for (auto bl:itsBlocks) bls.push_back(bl);
    IrrepIEClient::Init(bls);
    TBasisSetImplementation<double>::Insert(new IntegralEngine(this));    
//
//  Now insert the basis functions.
//
    MakeBasisFunctions(ns); //ns from PolarizedGaussianIEClient
};


//----------------------------------------------------------------
//
//  This contructor is used by Clone(RVec); only.
//
BasisSet::BasisSet(const BasisSet* bs,
        IntegralDataBase<double>* theDB,
        const optr_vector1<Block*>& theBlocks)
    : BasisSetImplementation(*bs)
    , TBasisSetImplementation<double>(bs->itsLAParams,theDB)
    , itsBlocks(theBlocks)
{
    // No UT coverage
    //MakeBasisFunctions(); //Compiler says these calls are ambiguous.  BUG
//    TBasisSetImplementation<double>::Insert(bs->GetIntegralEngine()->Clone());
}

void BasisSet::MakeBasisFunctions(const RVec& norms)
{
    EmptyBasisFunctions();
    size_t i=1;
    for (optr_vector1<Block*>::const_iterator bl(itsBlocks.begin()); bl!=itsBlocks.end(); bl++)
        for (std::vector<Polarization>::const_iterator p((*bl)->itsPols.begin()); p!=(*bl)->itsPols.end(); p++)
            BasisSetImplementation::Insert(new BasisFunction((*bl)->itsRadial,*p,norms(i++)));
}//Compiler says these calls are ambiguous.  BUG

std::ostream& BasisSet::Write(std::ostream& os) const
{
    // No UT coverage
    if (!Pretty())
    {
        os << itsBlocks;
        BasisSetImplementation::Write(os);
        TBasisSetImplementation<double>::Write(os);
    }
    else
    {
        for (optr_vector1<Block*>::const_iterator bl(itsBlocks.begin()); bl!=itsBlocks.end(); bl++)
            os << **bl;
    }
    return os;
}

std::istream& BasisSet::Read (std::istream& is)
{
    // No UT coverage
    is >> itsBlocks;
//    MakeBasisFunctions();
    BasisSetImplementation::Read(is);
    TBasisSetImplementation<double>::Read(is);
    return is;
}

IrrepBasisSet* BasisSet::Clone() const
{
    return new BasisSet(*this);
}

IrrepBasisSet* BasisSet::Clone(const RVec3& newCenter) const
{
    // No UT coverage
//    optr_vector1<Block*> newBlocks;
//    for (optr_vector1<Block*>::const_iterator b(itsBlocks.begin()); b!=itsBlocks.end(); b++)
//        newBlocks.push_back((*b)->Clone(newCenter));
//    return new BasisSet(this,GetDataBase()->Clone(),newBlocks);
    assert(false);
    return 0;
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

} //namespace PolarizedGaussian
