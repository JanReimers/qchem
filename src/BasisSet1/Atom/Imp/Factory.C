// File::BasisSet/Atom/Factory.C  Factory function for atom basis sets.
module; 
#include <cassert>
#include <nlohmann/json.hpp>

module qchem.BasisSet1.Atom.Factory;
import qchem.Symmetry.AtomEC;
// import qchem.BasisSet1.Atom.Slater.RKB.BS;
// import qchem.BasisSet1.Atom.Gaussian.NR.BS;
// import qchem.BasisSet1.Atom.Slater.NR.BS;
// import qchem.BasisSet1.Atom.BSpline.NR.BS;

import qchem.BasisSet1.Atom.BasisSet;
import BasisSet1.Atom.Evaluators.BSpline.BS;
import BasisSet.Atom.Gaussian_BS_Evaluator;
import BasisSet.Atom.Slater.NR.BS_Evaluator;
import BasisSet.Atom.Slater.RKB.IBS_Evaluator;
import BasisSet.Atom.Slater.NR.IBS_Evaluator;

import qchem.BasisSet1.DB_Cache;

using json = nlohmann::json;

namespace BasisSet1::Atom
{

BasisSet1::Real_BS* Factory(const nlohmann::json& js,size_t Z)
{
    return Factory(js,Atom_EC(Z));
}

BasisSet1::Real_BS* Factory(const nlohmann::json& js,const ElectronConfiguration& ec)
{
    if (BasisSet1::theGlobalCache==0)
        BasisSet1::theGlobalCache=new BasisSet1::IntegralsCache_RAM<double>(true);     
    Type type=js["type"].template get<Type>();
    const Atom_EC& aec=dynamic_cast<const Atom_EC&>(ec);
    size_t LMax=aec.GetLMax();
    BasisSet1::Real_BS* bs=0;
    switch (type)
    {
    case Type::Slater:
    {
        if (js.contains("exponents"))
        {
            auto es1=js["exponents"].template get<std::vector<double>>();
            rvec_t es(es1.size(),&es1[0]);
            bs=new BasisSet<Slater_BS_Evaluator>(es,ec);

        }
        else
        {
            size_t N=js["N"];
            double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
            bs=new BasisSet<Slater_BS_Evaluator>(N,emin,emax,ec);
        }
        break;
    }
    case Type::Gaussian:
    {
        if (js.contains("exponents"))
        {
            auto es1=js["exponents"].template get<std::vector<double>>();
            rvec_t es(es1.size(),&es1[0]);
            bs=new BasisSet<Gaussian_BS_Evaluator>(es,ec);
        }
        else
        {
            size_t N=js["N"];
            double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
            bs=new BasisSet<Gaussian_BS_Evaluator>(N,emin,emax,ec);
        }
        break;
    }
    case Type::BSpline6:
    {
        size_t N=js["N"];
        double rmin=js["rmin"].template get<double>(),rmax=js["rmax"].template get<double>();
        bs=new BasisSet<BSpline_BS_Evaluator<6>>(N,rmin,rmax,ec);
        break;   
    }
//     case Type::BSpline16:
//     {
//         double rmin=js["rmin"].template get<double>(),rmax=js["rmax"].template get<double>();
//         bs=new AtomBS::BSpline1::BasisSet<6,BSpline_r_BS_Evaluator>(N,rmin,rmax,ec);
//         break;   
//     } 
//     case Type::BSpline9:
//     {
//         double rmin=js["rmin"].template get<double>(),rmax=js["rmax"].template get<double>();
//         bs=new AtomBS::BSpline::BasisSet<9>(N,rmin,rmax,ec);
//         break;   
//     } 
// case Type::BSpliner6:
//     {
//         double rmin=js["rmin"].template get<double>(),rmax=js["rmax"].template get<double>();
//         bs=new AtomBS::BSpline::BasisSet_r<6>(N,rmin,rmax,ec);
//         break;   
//     } 
//     case Type::BSpliner9:
//     {
//         double rmin=js["rmin"].template get<double>(),rmax=js["rmax"].template get<double>();
//         bs=new AtomBS::BSpline::BasisSet_r<9>(N,rmin,rmax,ec);
//         break;   
//     } 

    case Type::Slater_RKB:
    {
        size_t N=js["N"];
        double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
        bs=new BasisSet_RKB<Slater_IBS_Evaluator,Slater_RKBS_IBS_Evaluator>(N,emin,emax,ec);
        break;
    }
    case Type::Gaussian_RKB:
    {
        size_t N=js["N"];
        double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
        bs=new BasisSet_RKB<Gaussian_IBS_Evaluator,Gaussian_RKBS_IBS_Evaluator>(N,emin,emax,ec);
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
