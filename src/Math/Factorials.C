// File: Factorials.C  Compile-time numeric lookup tables (factorials, double factorials, powers of two).
module;
#include <array>

export module qchem.Factorials;

// These are pure integer/double recurrences of fixed size, so they are built at COMPILE time by constexpr
// functions and exposed as immutable constants -- no runtime "Initializing..." pass, no forced static-init-
// before-main hack, no mutable globals.  Index range is [0,NMax]; callers index them like the old C arrays
// (Fact[k]), which std::array::operator[] supports unchanged.
export namespace qchem
{
    inline constexpr int NMax=19; //Should enable up to f orbitals.

    namespace detail
    {
        constexpr std::array<double,NMax+1> MakeFact()    //!< n!  : 1,1,2,6,24,...
        {
            std::array<double,NMax+1> a{};
            a[0]=1.0;
            for (int n=1; n<=NMax; n++) a[n]=a[n-1]*n;
            return a;
        }
        constexpr std::array<double,NMax+1> MakeDFact()   //!< n!! : 1,1,2,3,8,15,...
        {
            std::array<double,NMax+1> a{};
            a[0]=1.0; a[1]=1.0;
            for (int n=2; n<=NMax; n++) a[n]=a[n-2]*n;
            return a;
        }
        constexpr std::array<double,NMax+1> MakeTwon()     //!< 2^n : 1,2,4,8,...
        {
            std::array<double,NMax+1> a{};
            a[0]=1.0;
            for (int n=1; n<=NMax; n++) a[n]=a[n-1]*2;
            return a;
        }
    }

    inline constexpr std::array<double,NMax+1> Fact =detail::MakeFact();  //!< Factorials n! lookup table.
    inline constexpr std::array<double,NMax+1> DFact=detail::MakeDFact(); //!< Double factorials n!! lookup table.
    inline constexpr std::array<double,NMax+1> Twon =detail::MakeTwon();  //!< 2^n lookup table.
}
