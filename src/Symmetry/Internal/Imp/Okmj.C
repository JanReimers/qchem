module;
#include <iostream>
#include <iomanip>
#include <cassert>
#include <blaze/Math.h>
module qchem.Symmetry.Internal.Spherical;
import qchem.Common.Strings; //To get SPDFG string table.
import qchem.stl_io;

using std::cout;
using std::endl;
namespace Symmetry::Internal::Spherical
{

Ωκ::Ωκ(int _κ) : κ(_κ) 
{
    assert(abs(κ)<10);
};

size_t Ωκ::SequenceIndex() const //Used for op<
 {
    assert(abs(κ)<=LMax+1);
    return κ+LMax;
 }


size_t Ωκ::GetDegeneracy() const
{
    return Getj()+0.5; //(2j+1)/2 degeneracy for one spin state.
}


std::ostream& Ωκ::Write(std::ostream& os) const
{
    int jindex=Getj()-0.5;
    os << SPDFG[Getl()] << j2s[jindex] << " κ=" << std::setw(2) << κ << " ";
        
    return os;
}


Ωκmj::Ωκmj(int κ, const rvec_t& _mjs) : Ωκ(κ), mjs(_mjs) {};

size_t Ωκmj::SequenceIndex() const //Used for op<
 {
    assert(abs(κ)<=LMax+1);
    // double mjmax=*std::max_element(mjs.begin(),mjs.end());
    double mjmax=blaze::max(mjs);
    size_t offset=2*(LMax+1); //End of Ωκ sequence indexes.
    for (int k1=-(int)LMax-1;k1<κ;k1++) offset+=2*j(k1)+1; //add up all the degneracies below κ.
    return mjmax+Getj()+offset;
 }


size_t Ωκmj::GetDegeneracy() const
{
    return mjs.size();
}

std::ostream& Ωκmj::Write(std::ostream& os) const
{
    int jindex=Getj()-0.5;
    os << SPDFG[Getl()] << j2s[jindex] << " κ=" << std::setw(2) << κ << " mj=" << std::setw(4) << std::setprecision(1) << mjs << " ";
    return os;
}

} // namespace


