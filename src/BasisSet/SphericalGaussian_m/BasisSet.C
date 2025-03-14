// File SphericalGaussian_m/BasisSet.H

#include "Imp/BasisSet/SphericalGaussian_m/BasisSet.H"
#include "Imp/BasisSet/SphericalGaussian_m/IrrepBasisSet.H"
#include "Imp/BasisSet/PolarizedGaussian/Readers/Reader.H"
#include "Imp/BasisSet/PolarizedGaussian/RadialFunction.H"
#include "Imp/BasisSet/GaussianScaler.H"
#include "Imp/Integrals/AngularIntegrals.H"
#include <algorithm>

namespace SphericalGaussian_m
{


BasisSet::BasisSet(const LAParams& lap,size_t N, double emin, double emax, size_t LMax)
{
    GaussianScaler gs(N,emin,emax,LMax);
    const DB_BS_2E<double>* db=this;
    for (size_t L=0;L<=LMax;L++)
        for (int m=-L;m<=(int)L;m++)
            Insert(new Orbital_IBS(lap,db,gs.Get_es(L),L,m));            
        
}

using PolarizedGaussian::Reader;
using PolarizedGaussian::RadialFunction;

BasisSet::BasisSet(const LAParams& lap, Reader* reader, const Atom* atom)
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
        size_t L=le.first;
        for (int m=-L;m<=(int)L;m++)
            Insert(new Orbital_IBS(lap,this,le.second,L,m));   
    }
}



} //namespace
