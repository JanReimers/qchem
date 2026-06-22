// File::BasisSet/Atom/Factory.C  Factory function for atom basis sets.
module; 
#include <cassert>
#include <nlohmann/json.hpp>
#include <blaze/math/dense/DenseIterator.h> //In order for std::sort to work.

module qchem.BasisSet.Atom.Factory;
import qchem.ElectronConfiguration.AtomNR;
import qchem.ElectronConfiguration.AtomDirac;
import qchem.BasisSet.Atom.BasisSet;
import qchem.BasisSet.Atom.Evaluators.BSpline.IBS;
import qchem.BasisSet.Atom.Evaluators.BSpline.IBS_r;
import qchem.BasisSet.Atom.Evaluators.Gaussian.IBS; 
import qchem.BasisSet.Atom.Evaluators.Slater.IBS;
import qchem.Math;
import qchem.Blaze;

using json = nlohmann::json;

namespace BasisSet::Atom
{

using namespace Evaluators;

Real_BS* Factory(const nlohmann::json& js,size_t Z)
{
    Type type=js["type"].template get<Type>();
    if (type==Type::Slater_RKB || type==Type::Gaussian_RKB)
        return Factory(js,AtomDirac_EC(Z));
    else
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
            bs=new BasisSet_HF<Slater::NR_Evaluator>(es,aec,ltrim);

        }
        else
        {
            size_t N=js["N"];
            double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
            bs=new BasisSet_HF<Slater::NR_Evaluator>(N,emin,emax,aec);
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
            bs=new BasisSet_HF<Gaussian::NR_Evaluator>(es,aec,ltrim);
        }
        else
        {
            size_t N=js["N"];
            double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
            bs=new BasisSet_HF<Gaussian::NR_Evaluator>(N,emin,emax,aec);
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
        if (js.contains("exponents"))
        {
            auto es1=js["exponents"].template get<std::vector<double>>();
            size_t ltrim=js["ltrim"].template get<size_t>();
            rvec_t es(es1.size(),&es1[0]);
            bs=new BasisSet_RKB<Slater::RKBL_Evaluator,Slater::RKBS_Evaluator>(es,aec,ltrim);
        }
        else
        {
            size_t N=js["N"];
            double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
            bs=new BasisSet_RKB<Slater::RKBL_Evaluator,Slater::RKBS_Evaluator>(N,emin,emax,aec);
        }
        break;
    }
    case Type::Gaussian_RKB:
    {
        if (js.contains("exponents"))
        {
            auto es1=js["exponents"].template get<std::vector<double>>();
            size_t ltrim=js["ltrim"].template get<size_t>();
            rvec_t es(es1.size(),&es1[0]);
            bs=new BasisSet_RKB<Gaussian::RKBL_Evaluator,Gaussian::RKBS_Evaluator>(es,aec,ltrim);
        }
        else
        {
            size_t N=js["N"];
            double emin=js["emin"].template get<double>(),emax=js["emax"].template get<double>();
            bs=new BasisSet_RKB<Gaussian::RKBL_Evaluator,Gaussian::RKBS_Evaluator>(N,emin,emax,aec);
        }
        break;
    }

    }
    
    assert(bs);
    return bs;
}

rvec_t stride(const rvec_t& v, size_t s)
{
    rvec_t ret((v.size()+1)/s);
    auto i=v.begin();
    for (auto& iret:ret)
    {
        iret=*i;
        i+=s;
    }
    return ret;
}

rvec_t SlaterExponents(BasisSetAccuracy acc,size_t Z)
{
    using enum BasisSetAccuracy;
    // beta               needs to be empircally adjusted to get non singular overlap matrices (smin>=tol) for acc=High.
    //    Using SVD     decomposition of S, we need beta >=1.27 for Slater functions (tol=1e-13).
    //    Using Eigen   decomposition of S, we need beta >=1.27 for Slater functions (tol=1e-13)..
    //    Using Cholesky decomposition of S, we need beta >=1.23 for Slater functions (tol=0.00043)..
    // These betas are too low fo actuall calculations.  They get lost iterating in noise.
    // tol, emin and emax needs to be empircally adjusted to get good HF ground state energies for *all* atoms!
    // 
    // Good for He emin=0.1, beta=1.30, emax=500. NZ=floor(N-14+Z*14/100.);
    // Good for Rn emin=0.1, beta=1.40, emax=1500. NZ=floor(N-14+Z*14/100.);
    // with these setting Cas Z=20 does not work well at all though!  Too much exponent trimming.  Switch to trimming
    // based on (1.-sqrt(Z/100.) which give Ca a few more exponents to play with.
    double emin=0.1, beta=1.44, emax=1500.; //Assume Z=100 for emax this gets a good cusp for 1s Z=100 orbital.
    size_t N=floor(log(emax/emin)/log(beta));
    rvec_t exponents(N);
    double e=emin;
    for (auto i:iv_t(0,N))
    {
        exponents[i]=e;
        e*=beta;
    }
    //  prune for Z
    
    size_t NZ=floor(N-14*(1.-sqrt(Z/100.)));
    // std::cout << "Z,N,NZ=" << Z << " " << N << " " << NZ << std::endl;
    exponents=blazem::subvector(exponents,0,NZ);  //Trim off all the large exponents for smaller Z
    if (acc==BasisSetAccuracy::Medium)
    {
        size_t N1= NZ/6.;
        exponents=blazem::subvector(exponents,N1,NZ-N1-1); //Take awau N1 small exponents, and one large exponent.
    }
    else if(acc==BasisSetAccuracy::Low)
    {
        exponents=stride(exponents,2); //pick every other 
    }
    return exponents;
}

rvec_t GaussianExponents(BasisSetAccuracy acc,size_t Z)
{
    using enum BasisSetAccuracy;
    // beta               needs to be empircally adjusted to get non singular overlap matrices (smin>=tol) for acc=Extreme.
    //    Using SVD     decomposition of S, we need beta >=1.27 for Slater functions (tol=1e-13).
    //    Using Eigen   decomposition of S, we need beta >=1.27 for Slater functions (tol=1e-13)..
    //    Using Cholesky decomposition of S, we need beta >=1.23 for Slater functions (tol=0.00043)..
    // tol, emin and emax needs to be empircally adjusted to get good HF ground state energies for *all* atoms!
    // 
    // Good for He emin=0.1, beta=1.30, emax=500. NZ=floor(N-14+Z*14/100.);
    // Good for Rn emin=0.1, beta=1.40, emax=1500. NZ=floor(N-14+Z*14/100.);
    // with these setting Ca Z=20 does not work well at all though!
    double emin=0.01, beta=1.6, emax=4e7; //Assume Z=100 for emax this gets a good cusp for 1s Z=100 orbital.
    size_t N=floor(log(emax/emin)/log(beta));
    rvec_t exponents(N);
    double e=emin;
    for (auto i:iv_t(0,N))
    {
        exponents[i]=e;
        e*=beta;
    }
    //  prune for Z
    
    size_t NZ=floor(N-18*(1.-sqrt(Z/100.)));
    // std::cout << "Z,N,NZ=" << Z << " " << N << " " << NZ << std::endl;
    exponents=blazem::subvector(exponents,0,NZ);  //Trim off all the large exponents for smaller Z
    if (acc==BasisSetAccuracy::Medium)
    {
        size_t N1=NZ/12.-2;
        exponents=blazem::subvector(exponents,N1,NZ-2*N1); //Take away N1 small exponents, and one large exponent.
    }
    else if(acc==BasisSetAccuracy::Low)
    {
        exponents=stride(exponents,2); //pick every other 
    }
    return exponents;
}

Real_BS* Factory(BasisSetAccuracy acc, Type type,size_t Z)
{
    using enum BasisSetAccuracy;
    nlohmann::json js;
    switch (type)
    {
        case Type::Slater_RKB: 
        case Type::Slater: 
        {
            switch (acc)
            {
                case N3:
                    js={{"type",type},{"ltrim",0},{"exponents",{0.5,1,2.0}}};break;
                case N5:
                    js={{"type",type},{"ltrim",0},{"exponents",{0.25,0.5,1,2.0,4.0}}};break;
                case Low:
                case Medium:
                case High:
                    js={{"type",type},{"ltrim",1},{"exponents",SlaterExponents(acc,Z)}};break;
            }
            break;
        }
        case Type::Gaussian_RKB:
        case Type::Gaussian:
        {
            switch (acc)
            {
                case N3:
                    js={{"type",type},{"ltrim",0},{"exponents",{0.5,1,2.0}}};break;
                case N5:
                    js={{"type",type},{"ltrim",0},{"exponents",{0.25,0.5,1,2.0,4.0}}};break;
                case Low:
                case Medium:
                case High:
                    js={{"type",type},{"ltrim",1},{"exponents",GaussianExponents(acc,Z)}};break;
            }
            break;
            break;
        }
        case Type::BSpline6:
        case Type::BSpliner6:
        {
            switch (acc)
            {
                case N3:
                    js={{"type",type},{"rmin",0.5},{"rmax",2},{"N",3}};break;
                case N5:
                    js={{"type",type},{"rmin",0.25},{"rmax",4},{"N",5}};break;
                case Low:
                    js={{"type",type},{"rmin",0.03 },{"rmax",20},{"N",15}};break;
                case Medium:
                    js={{"type",type},{"rmin",0.003},{"rmax",40},{"N",30}};break;
                case High:
                    js={{"type",type},{"rmin",0.003},{"rmax",40},{"N",50}};break;
            }
            break;
        }
    }
    return Factory(js,Z);
}

} //namespace 
