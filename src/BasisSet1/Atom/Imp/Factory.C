// File::BasisSet/Atom/Factory.C  Factory function for atom basis sets.
module; 
#include <cassert>
#include <nlohmann/json.hpp>

module qchem.BasisSet.Atom.Factory;
import qchem.Symmetry.AtomEC;
import qchem.BasisSet.Atom.BasisSet;
import qchem.BasisSet.Atom.Evaluators.BSpline.IBS;
import qchem.BasisSet.Atom.Evaluators.Gaussian.IBS; 
import qchem.BasisSet.Atom.Evaluators.Slater.IBS;

import qchem.BasisSet.DB_Cache;

using json = nlohmann::json;

namespace BasisSet::Atom
{

using namespace Evaluators;

Real_BS* Factory(const nlohmann::json& js,size_t Z)
{
    return Factory(js,Atom_EC(Z));
}

Real_BS* Factory(const nlohmann::json& js,const Atom_EC& aec)
{
    if (::BasisSet::theGlobalCache==0)
        ::BasisSet::theGlobalCache=new ::BasisSet::IntegralsCache_RAM<double>(true);     
    Type type=js["type"].template get<Type>();
    Real_BS* bs=0;

    switch (type)
    {
    case Type::Slater:
    {
        if (js.contains("exponents"))
        {
            auto es1=js["exponents"].template get<std::vector<double>>();
            rvec_t es(es1.size(),&es1[0]);
            bs=new BasisSet_HF<Slater::Slater_IBS_Evaluator>(es,aec);

        }
        else
        {
            size_t N=js["N"];
            double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
            bs=new BasisSet_HF<Slater::Slater_IBS_Evaluator>(N,emin,emax,aec);
        }
        break;
    }
    case Type::Gaussian:
    {
        if (js.contains("exponents"))
        {
            auto es1=js["exponents"].template get<std::vector<double>>();
            rvec_t es(es1.size(),&es1[0]);
            bs=new BasisSet_HF<Gaussian::Gaussian_IBS_Evaluator>(es,aec);
        }
        else
        {
            size_t N=js["N"];
            double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
            bs=new BasisSet_HF<Gaussian::Gaussian_IBS_Evaluator>(N,emin,emax,aec);
        }
        break;
    }
    case Type::BSpline6:
    {
        size_t N=js["N"];
        double rmin=js["rmin"].template get<double>(),rmax=js["rmax"].template get<double>();
        bs=new BasisSet_1E_HF<BSpline::BSpline_IBS_Evaluator<6>>(N,rmin,rmax,aec);
        break;   
    }

    case Type::Slater_RKB:
    {
        size_t N=js["N"];
        double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
        bs=new BasisSet_RKB<Slater::Slater_IBS_Evaluator,Slater::Slater_RKBS_IBS_Evaluator>(N,emin,emax,aec);
        break;
    }
    case Type::Gaussian_RKB:
    {
        size_t N=js["N"];
        double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
        bs=new BasisSet_RKB<Gaussian::Gaussian_IBS_Evaluator,Gaussian::Gaussian_RKBS_IBS_Evaluator>(N,emin,emax,aec);
        break;
    }
//     case Type::BSpline_RKB:
//     {
//         assert(false);
//         //double rmin=js["rmin"].template get<double>(),rmax=js["rmax"].template get<double>();
//         // bs=new Atom_kappa::BSpline::BasisSet<6>(N,rmin,rmax,LMax);
//         break;   
//     } 




    }
    
    assert(bs);
    return bs;
}

} //namespace 
