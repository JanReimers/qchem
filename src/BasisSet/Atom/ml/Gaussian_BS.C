// File: Atom/ml/Gaussian_BS.H  r^l exp(-ar^2)*Y_lm type basis set.

#include "Imp/BasisSet/Atom/ml/Gaussian_BS.H"
#include "Imp/BasisSet/Atom/ml/Gaussian_IBS.H"
#include "Imp/BasisSet/PolarizedGaussian/Readers/Reader.H"
#include "Imp/BasisSet/PolarizedGaussian/RadialFunction.H"
#include "Imp/BasisSet/Atom/radial/Gaussian/ExponentScaler.H"
#include "Imp/Integrals/AngularIntegrals.H"
#include <algorithm>

namespace Atom_ml
{
namespace Gaussian
{


BasisSet::BasisSet(size_t N, double emin, double emax, size_t LMax)
{
    ::Gaussian::ExponentScaler gs(N,emin,emax,LMax);
    for (size_t L=0;L<=LMax;L++)
        for (int m=-L;m<=(int)L;m++)
            Insert(new Orbital_IBS(this,gs.Get_es(L),L,m));
        
}

using PolarizedGaussian::Reader;
using PolarizedGaussian::RadialFunction;

BasisSet::BasisSet(Reader* reader, const Atom* atom)
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
        Vector<double> es(le.second.size());
        int i=0;
        for (auto e:le.second) es(++i)=e; //Convert to vector.

        for (int m=-L;m<=(int)L;m++)
            Insert(new Orbital_IBS(this,es,L,m));   
    }
}



} //namespace
} //namespace
