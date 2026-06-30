// File: BasisSet/Atom/Evaluators/Internal/PascalTriangle.C
module;
#include <array>
#include <cstddef>
export module qchem.BasisSet.Atom.Internal.PascalTriangle;

namespace qchem {

//
//  Pascal's triangle of binomial coefficients C(row,j), 0 <= j <= row <= N, stored row-major in a flat
//  array (row r starts at r*(r+1)/2).  A pure integer recurrence of fixed size, so the whole table is built
//  at COMPILE time by the constexpr constructor and exposed as an immutable constant -- no runtime
//  "Initializing..." pass, no mutable singleton.
//
export class PascalTriangle
{
public:
    static constexpr int         N    = 9;                 //!< rows 0..N  (C(0,0) .. C(9,j))
    static constexpr std::size_t Size = (N+1)*(N+2)/2;

    constexpr PascalTriangle()
    {
        itsData[0]=1.0;
        for (int row=1; row<=N; row++)
        {
            (*this)(row,0)=(*this)(row,row)=1.0;
            for (int j=1; j<row; j++)
                (*this)(row,j)=(*this)(row-1,j-1)+(*this)(row-1,j);
        }
    }

    constexpr double  operator()(int row,int j) const {return itsData[Index(row,j)];}
    constexpr double& operator()(int row,int j)       {return itsData[Index(row,j)];}

    static const PascalTriangle thePascalTriangle;

private:
    static constexpr std::size_t Index(int row,int j) {return row*(row+1)/2+j;}

    std::array<double,Size> itsData{};
};

inline constexpr PascalTriangle PascalTriangle::thePascalTriangle{};

} // namespace qchem