// File: UnitTests/BasisSetPool.C  Define a restricted set of basis sets for all unit tests.
module;
#include <ranges>
#include <nlohmann/json.hpp>
#include <blaze/Math.h>
#include <cmath>
#include <iostream>
module qchem.Unittests.BasisSetPool;

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
    //    Using Cholsky decomposition of S, we need beta >=1.23 for Slater functions (tol=0.00043)..
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
    exponents=blaze::subvector(exponents,0,NZ);  //Trim off all the large exponents for smaller Z
    if (acc==BasisSetAccuracy::Medium)
    {
        size_t N1= NZ/6.;
        exponents=blaze::subvector(exponents,N1,NZ-N1-1); //Take awau N1 small exponents, and one large exponent.
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
    //    Using Cholsky decomposition of S, we need beta >=1.23 for Slater functions (tol=0.00043)..
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
    std::cout << "Z,N,NZ=" << Z << " " << N << " " << NZ << std::endl;
    exponents=blaze::subvector(exponents,0,NZ);  //Trim off all the large exponents for smaller Z
    if (acc==BasisSetAccuracy::Medium)
    {
        size_t N1=NZ/12.-2;
        exponents=blaze::subvector(exponents,N1,NZ-2*N1); //Take away N1 small exponents, and one large exponent.
    }
    else if(acc==BasisSetAccuracy::Low)
    {
        exponents=stride(exponents,2); //pick every other 
    }
    return exponents;
}

BasisSet::Real_BS* PoolFactory(BasisSetAccuracy acc, BasisSet::Atom::Type type,size_t Z)
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
        case Gaussian_RKB:
        case Gaussian:
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
        case BSpline6:
        case BSpliner6:
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
    return BasisSet::Atom::Factory(js,Z);
}

