// File::BasisSet/Atom/Factory.C  Factory function for atom basis sets.
module; 
#include <cassert>
#include <nlohmann/json.hpp>

module qchem.BasisSet.Atom.Factory;
import qchem.Symmetry.AtomEC;
import qchem.BasisSet.Atom.BasisSet;
import qchem.BasisSet.Atom.Evaluators.BSpline.IBS;
import qchem.BasisSet.Atom.Evaluators.BSpline.IBS_r;
import qchem.BasisSet.Atom.Evaluators.Gaussian.IBS; 
import qchem.BasisSet.Atom.Evaluators.Slater.IBS;

using json = nlohmann::json;

namespace BasisSet::Atom
{

using namespace Evaluators;

Real_BS* Factory(const nlohmann::json& js,size_t Z)
{
    return Factory(js,Atom_EC(Z));
}

Real_BS* Factory(const nlohmann::json& js,const ElectronConfiguration& aec)
{
    Type type=js["type"].template get<Type>();
    Real_BS* bs=0;

    switch (type)
    {
    case Type::Slater:
    {
        if (js.contains("exponents"))
        {
            auto es1=js["exponents"].template get<std::vector<double>>();
            size_t ltrim=js["ltrim"].template get<size_t>();
            rvec_t es(es1.size(),&es1[0]);
            bs=new BasisSet_HF<Slater::Evaluator>(es,aec,ltrim);

        }
        else
        {
            size_t N=js["N"];
            double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
            bs=new BasisSet_HF<Slater::Evaluator>(N,emin,emax,aec);
        }
        break;
    }
    case Type::Gaussian:
    {
        if (js.contains("exponents"))
        {
            auto es1=js["exponents"].template get<std::vector<double>>();
            size_t ltrim=js["ltrim"].template get<size_t>();
            rvec_t es(es1.size(),&es1[0]);
            bs=new BasisSet_HF<Gaussian::Evaluator>(es,aec,ltrim);
        }
        else
        {
            size_t N=js["N"];
            double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
            bs=new BasisSet_HF<Gaussian::Evaluator>(N,emin,emax,aec);
        }
        break;
    }
    case Type::BSpline6:
    {
        size_t N=js["N"];
        double rmin=js["rmin"].template get<double>(),rmax=js["rmax"].template get<double>();
        bs=new BasisSet_1E_HF<BSpline::Evaluator<6>>(N,rmin,rmax,aec);
        break;   
    }
    case Type::BSpliner6:
    {
        size_t N=js["N"];
        double rmin=js["rmin"].template get<double>(),rmax=js["rmax"].template get<double>();
        bs=new BasisSet_1E_HF<BSpline::Evaluator_r<6>>(N,rmin,rmax,aec);
        break;   
    }

    case Type::Slater_RKB:
    {
        size_t N=js["N"];
        double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
        bs=new BasisSet_RKB<Slater::Evaluator,Slater::RKBS_Evaluator>(N,emin,emax,aec);
        break;
    }
    case Type::Gaussian_RKB:
    {
        size_t N=js["N"];
        double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
        bs=new BasisSet_RKB<Gaussian::Evaluator,Gaussian::RKBS_Evaluator>(N,emin,emax,aec);
        break;
    }

    }
    
    assert(bs);
    return bs;
}

} //namespace 
