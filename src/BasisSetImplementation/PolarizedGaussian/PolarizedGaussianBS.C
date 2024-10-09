// File: PolarizedGaussianBS.C  Polarized Gaussian basis set, for MO calculations.



#include "BasisSetImplementation/PolarizedGaussian/PolarizedGaussianBF.H"
#include "BasisSetImplementation/PolarizedGaussian/PolarizedGaussianBS.H"
#include "BasisSetImplementation/PolarizedGaussian/PolarizedGaussianIE1.H"
#include "BasisSetImplementation/PolarizedGaussian/RadialFunctionReader.H"
#include "BasisSetImplementation/UnitSymmetryQN.H"
#include "BasisSetImplementation/NumericalIEImp.H"
#include "Cluster.H"
#include "Misc/ptr_vector1_io.h"
#include <cassert>
#include <algorithm> //Need std::max

std::vector<Polarization> MakePolarizations(const std::vector<int>& Ls);
template <class T> T Max(const std::vector<T>& v)
{
    return *std::max_element(v.begin(), v.end());
}



//#######################################################################
//
//  Concrete  gaussian basis set.
//
PolarizedGaussianBS::PolarizedGaussianBS()
    : BasisSetImplementation()
    , TBasisSetImplementation<double>()
{};

PolarizedGaussianBS::
PolarizedGaussianBS(const LinearAlgebraParams& lap,IntegralDataBase<double>* theDB, RadialFunctionReader* bsr, const Cluster* cl, Mesh * theMesh)
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
        BasisFunctionBlock* bfb=new BasisFunctionBlock(r,nbasis);
        for (auto& p:MakePolarizations(Ls[i]))
        {
            bfb->Add(p);
            nbasis++;
        }
        itsBlocks.push_back(bfb);
        i++;
    }
//
//  Now insert the basis functions.
//
    MakeBasisFunctions();
//
//  Make the integral engine.  Can't do this until all the basis functions and
//  blocks are in place.
//
    {
        PolarizedGaussianIE1::blocks_t bls;
        for (auto bl:itsBlocks) bls.push_back(bl);
        PolarizedGaussianIEClient::Init(bls);   
        TBasisSetImplementation<double>::Insert(new PolarizedGaussianIE1(this));    
    }

    if (theMesh)  
        TBasisSetImplementation<double>::Insert(new NumericalIEImp<double>(theMesh) );
        
};


//----------------------------------------------------------------
//
//  This contructor is used by Clone(RVec); only.
//
PolarizedGaussianBS::PolarizedGaussianBS(const PolarizedGaussianBS* bs,
        IntegralDataBase<double>* theDB,
        const optr_vector1<BasisFunctionBlock*>& theBlocks)
    : BasisSetImplementation(*bs)
    , TBasisSetImplementation<double>(bs->itsLAParams,theDB)
    , itsBlocks(theBlocks)
{
    // No UT coverage
    MakeBasisFunctions(); //Compiler says these calls are ambiguous.  BUG
//    TBasisSetImplementation<double>::Insert(bs->GetIntegralEngine()->Clone());
}

void PolarizedGaussianBS::MakeBasisFunctions()
{
    EmptyBasisFunctions();
    for (optr_vector1<BasisFunctionBlock*>::const_iterator bl(itsBlocks.begin()); bl!=itsBlocks.end(); bl++)
        for (std::vector<Polarization>::const_iterator p((*bl)->itsPols.begin()); p!=(*bl)->itsPols.end(); p++)
            BasisSetImplementation::Insert(new PolarizedGaussianBF((*bl)->itsRadial,*p));
}//Compiler says these calls are ambiguous.  BUG

std::ostream& PolarizedGaussianBS::Write(std::ostream& os) const
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
        for (optr_vector1<BasisFunctionBlock*>::const_iterator bl(itsBlocks.begin()); bl!=itsBlocks.end(); bl++)
            os << **bl;
    }
    return os;
}

std::istream& PolarizedGaussianBS::Read (std::istream& is)
{
    // No UT coverage
    is >> itsBlocks;
    MakeBasisFunctions();
    BasisSetImplementation::Read(is);
    TBasisSetImplementation<double>::Read(is);
    return is;
}

IrrepBasisSet* PolarizedGaussianBS::Clone() const
{
    return new PolarizedGaussianBS(*this);
}

IrrepBasisSet* PolarizedGaussianBS::Clone(const RVec3& newCenter) const
{
    // No UT coverage
//    optr_vector1<BasisFunctionBlock*> newBlocks;
//    for (optr_vector1<BasisFunctionBlock*>::const_iterator b(itsBlocks.begin()); b!=itsBlocks.end(); b++)
//        newBlocks.push_back((*b)->Clone(newCenter));
//    return new PolarizedGaussianBS(this,GetDataBase()->Clone(),newBlocks);
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

