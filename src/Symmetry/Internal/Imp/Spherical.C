// File: Symmetry/Internal/Imp/Spherical.C  The concrete atomic (spherical / Dirac-spinor) symmetry containers.
//
// Implementations of the four m-resolved / m-degenerate angular-symmetry quantum-number carriers that the
// `qchem.Symmetry.Internal.Spherical` module exports (all constructed via the Factory):
//   Yl    -- non-magnetic (m-degenerate) spherical harmonic Y_l                (l)
//   Ylm   -- m-resolved spherical harmonics                                    (l + a set of m_l)
//   Omega-kappa      -- non-magnetic Dirac spinor                              (kappa)
//   Omega-kappa-mj   -- m_j-resolved Dirac spinor                              (kappa + a set of m_j)
// Each is pure quantum-number storage: SequenceIndex() is just an ordering/caching key, Write() the label.
// (Merged from the former Yl.C / Ylm.C / Okmj.C -- three impl units of one module.)
module;
#include <iostream>
#include <iomanip>
#include <cassert>
#include <algorithm>
#include <blaze/math/dense/DenseIterator.h> // so std::sort can see the Blaze iterator op==/op!=
module qchem.Symmetry.Internal.Spherical;
import qchem.Strings;   // SPDFG / j2s label tables
import qchem.stl_io;
import qchem.Blaze;
import qchem.Math;

using std::cout;
using std::endl;

namespace qchem::Symmetry::Internal::Spherical
{

//================================================================================================
// Yl -- m-degenerate spherical harmonic symmetry (l only).
//================================================================================================
std::ostream& Yl::Write(std::ostream& os) const
{
    return os << SPDFG[itsL];
}

//================================================================================================
// Ylm -- m-resolved spherical harmonics (l + a set of m_l).
//================================================================================================
Ylm::Ylm(size_t l, const ivec_t& _mls)
: itsL(l), mls(_mls)
{
    assert(mls.size()>0);
};

size_t Ylm::SequenceIndex() const //Used for op<
{
    static size_t start=LMax+1;  //Start after all the Yl Sequence Indexes

    // Sort mls to get canonical ordering (lexicographic)
    auto sorted_mls = mls;
    std::sort(sorted_mls.begin(), sorted_mls.end());

    // Large offset per l value to ensure no collisions between different l values
    size_t l_offset = start + (size_t)itsL * 100000000UL;

    // Encode sorted_mls as a positional number with base (2*l+2)
    // This ensures each unique sorted_mls maps to a unique index
    size_t mls_code = 0;
    for (int ml : sorted_mls) {
        int ml_encoded = ml + (int)itsL + 1;  // Map [-l, +l] to [0, 2*l+1]
        mls_code = mls_code * (2*(int)itsL + 2) + ml_encoded;
    }

    return l_offset + mls_code;
}

size_t Ylm::GetDegeneracy() const
{
    return mls.size();
}

static inline size_t width(int ml) {return ml<0 ? 2 : 1;}
std::ostream& Ylm::Write(std::ostream& os) const
{
    os << SPDFG[itsL] << " ";
    if (mls.size()<2*(size_t)itsL+1)
    {
        os << "{";
        for (auto ml:mls)
        {
            size_t w=width(ml);
            os << std::setw(w) << ml << " ";
        }
        os << "}";
    }

    return os;
}

//================================================================================================
// Omega-kappa -- m-degenerate Dirac spinor (kappa only).
//================================================================================================
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
    os << SPDFG[Getl()] << j2s[jindex] << " κ=" << std::setw(2) << κ;

    return os;
}

//================================================================================================
// Omega-kappa-mj -- m_j-resolved Dirac spinor (kappa + a set of m_j).
//================================================================================================
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
    int mj_steps = static_cast<int>(lround(2.0 * j_val));
    size_t value_range = static_cast<size_t>(2 * mj_steps + 1);
    size_t size_range = 16;

    size_t code = kappa_code;
    code = code * size_range + sorted_mjs.size();
    for (double mj : sorted_mjs) {
        int mj_code = static_cast<int>(lround((mj + j_val) * 2.0));
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
