// File: UnitTests/BasisSetPool.H  Define a restricted set of basis sets for all unittests.  
// The object here is to define restrict collection of atom basis sets so that the global DB_Cache to not get too large.
#include <ranges>
#include <nlohmann/json.hpp>
#include <blaze/Math.h>
#include <cmath>
#include <iostream>
import qchem.BasisSet.Atom.Factory;

// _GLIBCXX_RANGES
// Extreme: As big and dense as you can get without becoming numercially unstable.
// High:    Extrame with 1/4 of the largest exponents removed.  Giving up core cusp but not exponent density
// Medium   Extreme with every other expoenent removed.  Giving up density but not core cusp.
// Low      Extremen with 1/4 of the largest exponents removed *and* every other expoenent removed.
//
enum class BasisSetAccuracy {N3,N5,Low,Medium,High,Extreme};
std::string BasisSetAccuracyStrs[]={"N3     ","N5     ","Low    ", "Medium ","High   ","Extreme"};

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


//
//  THe proceedure here is to make a big list for the extreme case, for Z=100, emax=500.
//  which is then trimmed based on Z.  And then further trimmed for lower accuracy cases.
//  Results shoud alway be even temperated and as much exponent sharing as possible,
//
rvec_t SlaterExponents(BasisSetAccuracy acc,size_t Z)
{
    using enum BasisSetAccuracy;
    // beta               needs to be empircally adjusted to get non singular overlap matrices (smin>=tol) for acc=Extreme.
    //    Using SVD     decomposition of S, we need beta >=1.27 for Slater functions (tol=1e-13).
    //    Using Eigen   decomposition of S, we need beta >=1.27 for Slater functions (tol=1e-13)..
    //    Using Cholsky decomposition of S, we need beta >=1.23 for Slater functions (tol=0.00043)..
    // tol, emin and emax needs to be empircally adjusted to get good HF ground state energies for *all* atoms!
    // 
    double emin=0.25, beta=1.23, emax=500.; //Assume Z=100 for emax this gets a good cusp for 1s Z=100 orbital.
    size_t N=floor(log(emax/emin)/log(beta));
    rvec_t exponents(N);
    double e=emin;
    for (auto i:iv_t(0,N))
    {
        exponents[i]=e;
        e*=beta;
    }
    //  prune for Z
    size_t NZ=floor(N-25+Z/4.);
    // std::cout << "Z,N,NZ=" << Z << " " << N << " " << NZ << std::endl;
    exponents=blaze::subvector(exponents,0,NZ);  //Trim off all the large exponents for smaller Z
    if (acc==BasisSetAccuracy::High)
    {
        size_t N1= 3*NZ/4;
        exponents=blaze::subvector(exponents,0,N1); //Trim 8 biggrst exponents.
    }
    else if(acc==BasisSetAccuracy::Medium)
    {
        exponents=stride(exponents,2); //pick every other 
    }
    else if(acc==BasisSetAccuracy::Low)
    {
        size_t N1= NZ>=16 ? NZ-8 : NZ/2;
        exponents=blaze::subvector(exponents,0,N1); //Trim 8 biggest exponents and then pick every other 
    }
    return exponents;
}
BasisSet::Real_BS* Factory(BasisSetAccuracy acc, BasisSet::Atom::Type type,size_t Z)
{
    using enum BasisSet::Atom::Type;
    using enum BasisSetAccuracy;
    nlohmann::json js;
    switch (type)
    {
        case Slater_RKB: 
        case Slater: 
        {
            switch (acc)
            {
                case N3:
                    js={{"type",type},{"exponents",{0.5,1,2.0}}};break;
                case N5:
                    js={{"type",type},{"exponents",{0.25,0.5,1,2.0,4.0}}};break;
                case Low:
                case Medium:
                case High:
                case Extreme:
                    js={{"type",type},{"exponents",SlaterExponents(acc,Z)}};break;
            }
            break;
        }
        case Gaussian_RKB:
        case Gaussian:
        {
            break;
        }
        case BSpline6:
        {
            break;
        }
    }
    return BasisSet::Atom::Factory(js,Z);
}

