// File: BasisSet.C  Polarized Gaussian basis set, for MO calculations.



#include "Imp/BasisSet/PolarizedGaussian/BasisFunction.H"
#include "Imp/BasisSet/PolarizedGaussian/IrrepBasisSet.H"
#include "Imp/BasisSet/PolarizedGaussian/IntegralEngine.H"
//#include "Imp/BasisSet/PolarizedGaussian/Readers/Reader.H"
#include "Imp/BasisSet/PolarizedGaussian/Readers/Gaussian94.H"
#include "Imp/BasisSet/PolarizedGaussian/Radial/GaussianRF.H"
#include "Imp/BasisSet/SphericalGaussian/QuantumNumber.H"
#include "Imp/BasisSet/GaussianScaler.H"
#include "Imp/Cluster/Atom.H"
#include <UnitSymmetryQN.H>
#include <Cluster.H>
#include "Imp/Containers/ptr_vector_io.h"
#include <cassert>
#include <algorithm> //Need std::max

template <class T> void FillPower(Vector<T>& arr,T start, T stop);

namespace PolarizedGaussian
{

std::vector<Polarization> MakePolarizations(const std::vector<int>& Ls);
std::vector<Polarization> MakePolarizations(int L);
template <class T> T Max(const std::vector<T>& v)
{
    return *std::max_element(v.begin(), v.end());
}



//#######################################################################
//
//  Concrete  gaussian basis set.
//
IrrepBasisSet::IrrepBasisSet()
    : IrrepBasisSetCommon()
    , TIrrepBasisSetCommon<double>()
{};

IrrepBasisSet::
IrrepBasisSet(const LAParams& lap,IntegralDataBase<double>* theDB, Reader* bsr, const Cluster* cl)
    : IrrepBasisSetCommon(new UnitSymmetryQN)
    , TIrrepBasisSetCommon<double>(lap,theDB)
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
    
    std::vector<const Block*> bls;
    for (auto bl:itsBlocks) bls.push_back(bl);
    IrrepIEClient::Init(bls);
    TIrrepBasisSetCommon<double>::Insert(new IntegralEngine());    
//
//  Now insert the basis functions.
//
    MakeBasisFunctions(ns); //ns from PolarizedGaussianIEClient
};



IrrepBasisSet::
IrrepBasisSet(const LAParams& lap,IntegralDataBase<double>* theDB,   size_t N, double emin, double emax, size_t LMax, const Cluster* cl)
    : IrrepBasisSetCommon(new UnitSymmetryQN)
    , TIrrepBasisSetCommon<double>(lap,theDB)
{
    GaussianScaler gs(N,emin,emax,LMax);
    int nbasis=1;
    for (auto atom:*cl)
    {
        for (size_t L=0;L<=LMax;L++)
        {
            Vector<double> es(gs.N(L));
            FillPower(es,gs.emin(L),gs.emax(L));
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
                itsBlocks.push_back(bfb);
            }
            
        }
    }
    std::vector<const Block*> bls;
    for (auto bl:itsBlocks) bls.push_back(bl);
    IrrepIEClient::Init(bls);
    TIrrepBasisSetCommon<double>::Insert(new IntegralEngine());    
//
//  Now insert the basis functions.
//
    MakeBasisFunctions(ns); //ns from PolarizedGaussianIEClient
    
}

// Single atom version
IrrepBasisSet::
IrrepBasisSet(const LAParams& lap,IntegralDataBase<double>* theDB,   size_t N, double emin, double emax, size_t L)
    : IrrepBasisSetCommon(new UnitSymmetryQN())
    , TIrrepBasisSetCommon<double>(lap,theDB)
{
    Vector<double> es(N);
    FillPower(es,emin,emax);
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
        itsBlocks.push_back(bfb);
    }
    std::vector<const Block*> bls;
    for (auto bl:itsBlocks) bls.push_back(bl);
    IrrepIEClient::Init(bls);
    TIrrepBasisSetCommon<double>::Insert(new IntegralEngine());    
//
//  Now insert the basis functions.
//
    MakeBasisFunctions(ns); //ns from PolarizedGaussianIEClient
}
//----------------------------------------------------------------
//
//  This contructor is used by Clone(RVec); only.
//
IrrepBasisSet::IrrepBasisSet(const IrrepBasisSet* bs,
        IntegralDataBase<double>* theDB,
        const optr_vector1<Block*>& theBlocks)
    : IrrepBasisSetCommon(*bs)
    , TIrrepBasisSetCommon<double>(bs->itsLAParams,theDB)
    , itsBlocks(theBlocks)
{
    // No UT coverage
    //MakeBasisFunctions(); //Compiler says these calls are ambiguous.  BUG
//    TBasisSetImplementation<double>::Insert(bs->GetIntegralEngine()->Clone());
}

void IrrepBasisSet::MakeBasisFunctions(const RVec& norms)
{
    EmptyBasisFunctions();
    size_t i=1;
    for (optr_vector1<Block*>::const_iterator bl(itsBlocks.begin()); bl!=itsBlocks.end(); bl++)
        for (std::vector<Polarization>::const_iterator p((*bl)->itsPols.begin()); p!=(*bl)->itsPols.end(); p++)
            IrrepBasisSetCommon::Insert(new BasisFunction((*bl)->itsRadial,*p,norms(i++)));
}//Compiler says these calls are ambiguous.  BUG


IrrepBasisSet* IrrepBasisSet::CreateCDFitBasisSet(const Cluster* cl) const
{
    // The A1 files support Z=1-54 (H-Te)  A2 version only go up to Zn
    PolarizedGaussian::Gaussian94Reader reader("../BasisSetData/A1_coul.bsd");
    return new IrrepBasisSet(itsLAParams,GetDataBase(),&reader,cl);
}

IrrepBasisSet* IrrepBasisSet::CreateVxcFitBasisSet(const Cluster* cl) const
{
    // The A1 files support Z=1-54 (H-Te)  A2 version only go up to Zn
    PolarizedGaussian::Gaussian94Reader reader("../BasisSetData/A1_exch.bsd");
    return new IrrepBasisSet(itsLAParams,GetDataBase(),&reader,cl);    
}


std::ostream& IrrepBasisSet::Write(std::ostream& os) const
{
    // No UT coverage
    IrrepBasisSetCommon::Write(os);
    TIrrepBasisSetCommon<double>::Write(os);
    if (!Pretty())
    {
        os << itsBlocks;
       
    }
    else
    {
//        for (optr_vector1<Block*>::const_iterator bl(itsBlocks.begin()); bl!=itsBlocks.end(); bl++)
//            os << **bl;
    }
    return os;
}

std::istream& IrrepBasisSet::Read (std::istream& is)
{
    // No UT coverage
    is >> itsBlocks;
//    MakeBasisFunctions();
    IrrepBasisSetCommon::Read(is);
    TIrrepBasisSetCommon<double>::Read(is);
    return is;
}

IrrepBasisSet* IrrepBasisSet::Clone() const
{
    return new IrrepBasisSet(*this);
}

IrrepBasisSet* IrrepBasisSet::Clone(const RVec3& newCenter) const
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
std::vector<Polarization> MakePolarizations(int L)
{
    std::vector<Polarization> ret;
        for(int m=0; m<=L; m++)
            for(int l=0; l<=L-m; l++)
                ret.push_back(Polarization(L-m-l,l,m));

    return ret;
}

} //namespace PolarizedGaussian
