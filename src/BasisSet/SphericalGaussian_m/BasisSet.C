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
    GaussianScaler ss(N,emin,emax,LMax);
    for (int L=0;L<=LMax;L++)
    {
        size_t  NL=ss.N(L);
        for (int m=-L;m<=L;m++)
        {
            IrrepBasisSet* ibs=0;
            if (L==0)
                ibs=new IrrepBasisSet(lap,GetDataBase(),N,ss.emin(L),ss.emax(L),L,m);
            else
            {
                const AtomIrrepIEClient* ibs0=(*this)[1];
//                ibs=new IrrepBasisSet(lap,GetDataBase(),N,ibs0->es(1),ibs0->es(N),L,m);             
                ibs=new IrrepBasisSet(lap,GetDataBase(),N-2*L,ibs0->es(L+1),ibs0->es(N-L),L,m);             
            }
            Append(ibs);
            Insert(ibs);            
        }
    }
        
    
        
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
//                if (l==2)
//                {
//                    std::set<double>& es0=Lexponents[0];
//                    for (auto& e:esL) es0.insert(e); //Stick all the d exponents into the s list
//                }
            }
        }
    }
    for (auto& le:Lexponents)
    {
        int L=le.first;
        for (int m=-L;m<=L;m++)
        {
            IrrepBasisSet* ibs=new IrrepBasisSet(lap,GetDataBase(),le.second,L,m);
            Append(ibs); //IECleint
            Insert(ibs); //Common with optr_vector     
        }
    }
}


} //namespace
