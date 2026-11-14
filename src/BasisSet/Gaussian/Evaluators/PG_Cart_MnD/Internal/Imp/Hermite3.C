// File: Hermite3.C   Class for managing 3 function Hermite coefficients.
module;

#include <iostream>
#include <cassert>
#include <vector>

module qchem.BasisSet.Gaussian.Evaluators.PG_Cart_MnD.Internal.Hermite;
namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD
{
//----------------------------------------------------------------------------------------
//
//  None of this will make any sense unless you read McMurchie and Davidson,
//  J. Comp. Phys. 26, 218 (1978).  Even then it may not make sense!!!
//


//----------------------------------------------------------------------------------------
//
//  Construction zone.
//

Hermite3::Hermite3()
    : itsLA   (0)
    , itsLB   (0)
    , itsLC   (0)
    , itsScale(-1.0)
{};

// The constructor, GetAny and operator() are INLINE in the module interface (the 3-centre twin Enzyme differentiates builds them).

std::ostream& operator<<(std::ostream& os,const Hermite3& h3)
{
    os << "scale " << h3.itsScale << ", LA " << h3.itsLA << ", LB " << h3.itsLB << ", LC " << h3.itsLC << std::endl;

    os << "N=" << 0 << std::endl;
    for (int na=0; na<=h3.itsLA; na++)
    {
        for (int nb=0; nb<=h3.itsLB; nb++)
        {
            for (int nc=0; nc<=h3.itsLC; nc++)
            {
                Polarization p(na,nb,nc),z(0,0,0);
                os << h3(p,p,z) << " ";
            }
            os << std::endl;
        }
        os << std::endl;
    }
    os << std::endl;

    os << "N=" << 0 << std::endl;
    for (int na=0; na<=h3.itsLA; na++)
    {
        for (int nb=0; nb<=h3.itsLB; nb++)
        {
            for (int nc=0; nc<=h3.itsLC; nc++) os << h3.d[h3.Index(0,na,nb,nc)] << " ";
            os << std::endl;
        }
        os << std::endl;
    }
    os << std::endl;
    os << "L=" << 0 << std::endl;
    for (int na=0; na<=h3.itsLA; na++)
    {
        for (int nb=0; nb<=h3.itsLB; nb++)
        {
            for (int nc=0; nc<=h3.itsLC; nc++) os << h3.e[h3.Index(0,na,nb,nc)] << " ";
            os << std::endl;
        }
        os << std::endl;
    }
    os << std::endl;
    os << "M=" << 0 << std::endl;
    for (int na=0; na<=h3.itsLA; na++)
    {
        for (int nb=0; nb<=h3.itsLB; nb++)
        {
            for (int nc=0; nc<=h3.itsLC; nc++) os << h3.f[h3.Index(0,na,nb,nc)] << " ";
            os << std::endl;
        }
        os << std::endl;
    }
    os << std::endl;

    /*for (int N=0; N<=h3.itsLA+h3.itsLB+h3.itsLC;N++)
    {
    os << "N=" << N << std::endl;
      for (int na=0;na<=h3.itsLA;na++)
      {
        for (int nb=0;nb<=h3.itsLB;nb++)
        {
    for (int nc=0;nc<=h3.itsLC;nc++) os << h3.d[Index(N,na,nb,nc)] << " ";
    os << std::endl;
        }
        os << std::endl;
      }
      os << std::endl;
    }
    for (int N=0; N<=h3.itsLA+h3.itsLB+h3.itsLC;N++)
    {
    os << "L=" << N << std::endl;
      for (int na=0;na<=h3.itsLA;na++)
      {
        for (int nb=0;nb<=h3.itsLB;nb++)
        {
    for (int nc=0;nc<=h3.itsLC;nc++) os << h3.e[Index(N,na,nb,nc)] << " ";
    os << std::endl;
        }
        os << std::endl;
      }
      os << std::endl;
    }
    for (int N=0; N<=h3.itsLA+h3.itsLB+h3.itsLC;N++)
    {
    os << "M=" << N << std::endl;
      for (int na=0;na<=h3.itsLA;na++)
      {
        for (int nb=0;nb<=h3.itsLB;nb++)
        {
    for (int nc=0;nc<=h3.itsLC;nc++) os << h3.f[Index(N,na,nb,nc)]  << " ";
    os << std::endl;
        }
        os << std::endl;
      }
      os << std::endl;
    }
    */
    return os;
}

} //namespace qchem::BasisSet::Gaussian::Evaluators::PG_Cart_MnD
