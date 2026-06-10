module;
#include <iostream>
#include <iomanip>
#include <cassert>
#include <algorithm>
#include <cmath>
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


Ωκmj::Ωκmj(int κ, const rvec_t& _mjs) : κ(κ), mjs(_mjs) {};

size_t Ωκmj::SequenceIndex() const //Used for op<
 {
    assert(abs(κ)<=LMax+1);

    // Keep Ωκmj indices separate from pure Ωκ indices.
    constexpr size_t pure_omega_offset = 1024;

    // Canonical κ code (non-negative)
    size_t kappa_code = static_cast<size_t>(κ + static_cast<int>(LMax) + 1);

    auto sorted_mjs = mjs;
    std::sort(sorted_mjs.begin(), sorted_mjs.end());

    double j_val = Getj();
    int mj_steps = static_cast<int>(std::lround(2.0 * j_val));
    size_t value_range = static_cast<size_t>(2 * mj_steps + 1);
    size_t size_range = 16;

    size_t code = kappa_code;
    code = code * size_range + sorted_mjs.size();
    for (double mj : sorted_mjs) {
        int mj_code = static_cast<int>(std::lround((mj + j_val) * 2.0));
        assert(mj_code >= 0 && static_cast<size_t>(mj_code) < value_range);
        code = code * value_range + static_cast<size_t>(mj_code);
    }

    return pure_omega_offset + code;
 }


size_t Ωκmj::GetDegeneracy() const
{
    return mjs.size();
}

std::ostream& Ωκmj::Write(std::ostream& os) const
{
    int jindex=Getj()-0.5;
    os << SPDFG[Getl()] << j2s[jindex] << " κ=" << std::setw(2) << κ << " mj={";
    for (double mj:mjs) os << std::fixed << std::setw(4) << std::setprecision(1) << mj << " ";
    os << "}";
    return os;
}

} // namespace


