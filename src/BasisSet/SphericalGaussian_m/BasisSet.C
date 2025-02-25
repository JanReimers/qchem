// File SphericalGaussian_m/BasisSet.H

#include "Imp/BasisSet/SphericalGaussian_m/BasisSet.H"
#include "Imp/BasisSet/SphericalGaussian_m/IrrepBasisSet.H"
#include "Imp/BasisSet/SphericalGaussian_m/IntegralEngine.H"
#include "Imp/BasisSet/PolarizedGaussian/Readers/Reader.H"
#include "Imp/BasisSet/PolarizedGaussian/RadialFunction.H"
#include "Imp/BasisSet/GaussianScaler.H"
#include <algorithm>

namespace SphericalGaussian_m
{


BasisSet::BasisSet(const LAParams& lap,size_t N, double emin, double emax, size_t LMax)
: BasisSetImp(new IntegralEngine) // this makes a integral DB
{
    GaussianScaler gs(N,emin,emax,LMax);
    for (size_t L=0;L<=LMax;L++)
        for (int m=-L;m<=(int)L;m++)
            Insert(new IrrepBasisSet(lap,GetDataBase(),gs.Get_es(L),L,m));            
        
}

using PolarizedGaussian::Reader;
using PolarizedGaussian::RadialFunction;

BasisSet::BasisSet(const LAParams& lap, Reader* reader, const Atom* atom)
: BasisSetImp(new IntegralEngine) // this makes a integral DB
{
    std::map<int,std::set<double> > Lexponents;
    reader->FindAtom(*atom);
    RadialFunction* rf=0;
    while ((rf=reader->ReadNext(*atom))) //Read in the radial function.
    {
        std::set<double> es=rf->GetExponents();
        std::vector<int> Ls=reader->GetLs();
        for (auto l:Ls)
        {
            if (const auto& i=Lexponents.find(l);i==Lexponents.end())
            {
                Lexponents[l]=es;
            }
            else
            {
                std::set<double>& esL=i->second;
                std::set<double> unions;
                set_union(esL.begin(), esL.end(), es.begin(),es.end(), inserter(unions, unions.begin()));
                i->second = unions;
            }
        }
    }
    for (auto& le:Lexponents)
    {
        int L=le.first;
        for (int m=-L;m<=L;m++)
            Insert(new IrrepBasisSet(lap,GetDataBase(),le.second,L,m)); //Common with optr_vector     
    }
}


} //namespace
