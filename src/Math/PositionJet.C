// File: Math/PositionJet.C  PositionJet<T>: a value and its first derivatives in the atom positions, sparse in atoms (lane B, 2026-11-08;
// doc/GeometryOptimizationPlan.md §13.1, approved by the user).
//
// THE OBJECT.  The first-order JET of a quantity in the atom positions R_1..R_N (Ehresmann's term in differential geometry): its value and
// its differential, here {v, ∂v/∂R_A}.  A Rall number, or a dual number, is the same object: first-order forward differentiation by
// arithmetic.  Every derivative-aware face returns it, and chaining across a face is ordinary arithmetic on it, so a library boundary or a
// virtual call costs the derivative nothing (the user's Rall-algebra input, 2026-11-08).  Enzyme fills a jet INSIDE an Imp; the chain rule
// between faces is written ONCE, here.
//
// SPARSE IN ATOMS.  A basis function depends on its own centre (and on a mesh point's owner), a density at a point on the atoms whose
// functions reach it.  So a jet holds only the atoms it depends on, sorted by atom index, with no duplicates.  An absent atom has a zero
// derivative, exactly.  A dense 3N jet would be plain forward mode, ∝ 3N on every value (§13.1's cost section).
//
// THE TWO DIRECTIONS.  Forward: build jets of a face's outputs and combine them.  Pullback: a scalar E = Σ_g a_g x_g, with the adjoint
// a_g = ∂E/∂x_g already known (a potential, fit coefficients), is AddWeightedTo(F, a_g) over the x_g's jets.
//
// Second order (a Hessian-vector product) is deferred (§13.1).
module;
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>
export module qchem.PositionJet;
export import qchem.Vector3D;

export namespace qchem
{

template <class T> class PositionJet
{
public:
    //! One atom's part of the differential: ∂value/∂R_atom.
    struct Entry { size_t atom; Vector3D<T> d; };

    PositionJet() = default;                                           //!< zero, independent of every position
    explicit PositionJet(T v) : itsValue(v) {}                         //!< a constant: no derivatives
    PositionJet(T v, size_t atom, const Vector3D<T>& d) : itsValue(v) {if (d.x!=T(0) || d.y!=T(0) || d.z!=T(0)) itsD.push_back({atom,d});}
    //! The coordinate \a axis (0, 1, 2 = x, y, z) of atom \a atom, at value \a v: ∂/∂R_atom is the unit vector on that axis.
    static PositionJet Coordinate(size_t atom, int axis, T v)
    {
        if (axis<0 || axis>2) throw std::invalid_argument("PositionJet::Coordinate: axis is 0, 1 or 2");
        return PositionJet(v, atom, Vector3D<T>(axis==0, axis==1, axis==2));
    }

    T Value() const {return itsValue;}
    //! The non-zero part of the differential, sorted by atom, no duplicates.
    const std::vector<Entry>& Derivatives() const {return itsD;}
    //! ∂value/∂R_atom (zero for an atom the jet does not depend on).
    Vector3D<T> D(size_t atom) const
    {
        auto it=std::lower_bound(itsD.begin(), itsD.end(), atom, [](const Entry& e, size_t a){return e.atom<a;});
        return (it!=itsD.end() && it->atom==atom) ? it->d : Vector3D<T>(0,0,0);
    }

    //! THE PULLBACK STEP: F[A] += a ∂value/∂R_A, for every atom the jet depends on (F sized to the structure's atom count).
    void AddWeightedTo(std::vector<Vector3D<T>>& F, T a) const
    {
        for (const Entry& e : itsD)
        {
            if (e.atom>=F.size()) throw std::out_of_range("PositionJet::AddWeightedTo: atom "+std::to_string(e.atom)+" is past the force array");
            F[e.atom]+=a*e.d;
        }
    }

    // --- arithmetic: the chain rule, once ---
    PositionJet& operator+=(const PositionJet& b) {Merge(b, T(1), T(1)); itsValue+=b.itsValue; return *this;}
    PositionJet& operator-=(const PositionJet& b) {Merge(b, T(1), T(-1)); itsValue-=b.itsValue; return *this;}
    PositionJet& operator*=(const PositionJet& b) {Merge(b, b.itsValue, itsValue); itsValue*=b.itsValue; return *this;}   // (uv)' = u'v + uv'
    PositionJet& operator*=(T s) {itsValue*=s; for (Entry& e : itsD) e.d*=s; return *this;}
    PositionJet& operator/=(T s) {return *this*=(T(1)/s);}
    PositionJet& operator+=(T s) {itsValue+=s; return *this;}
    PositionJet& operator-=(T s) {itsValue-=s; return *this;}

    friend PositionJet operator+(PositionJet a, const PositionJet& b) {return a+=b;}
    friend PositionJet operator-(PositionJet a, const PositionJet& b) {return a-=b;}
    friend PositionJet operator*(PositionJet a, const PositionJet& b) {return a*=b;}
    friend PositionJet operator/(const PositionJet& a, const PositionJet& b) {return a*Reciprocal(b);}
    friend PositionJet operator*(PositionJet a, T s) {return a*=s;}
    friend PositionJet operator*(T s, PositionJet a) {return a*=s;}
    friend PositionJet operator/(PositionJet a, T s) {return a/=s;}
    friend PositionJet operator+(PositionJet a, T s) {return a+=s;}
    friend PositionJet operator-(PositionJet a, T s) {return a-=s;}
    friend PositionJet operator-(PositionJet a) {return a*=T(-1);}
    friend PositionJet operator+(T s, PositionJet a) {return a+=s;}
    friend PositionJet operator-(T s, PositionJet a) {a*=T(-1); return a+=s;}

    //! THE CHAIN RULE through a scalar function: the jet of f(x), given f(x) = \a fx and f′(x) = \a dfx.
    friend PositionJet Chain(PositionJet x, T fx, T dfx)
    {
        for (Entry& e : x.itsD) e.d*=dfx;
        x.itsValue=fx;
        return x;
    }
    friend PositionJet Reciprocal(const PositionJet& x) {return Chain(x, T(1)/x.itsValue, -T(1)/(x.itsValue*x.itsValue));}
    friend PositionJet exp (const PositionJet& x) {const T e=std::exp(x.itsValue); return Chain(x, e, e);}
    friend PositionJet sqrt(const PositionJet& x) {const T s=std::sqrt(x.itsValue); return Chain(x, s, T(0.5)/s);}

private:
    //! this ← sa this + sb b on the differential: a merge of the two sorted atom lists.
    void Merge(const PositionJet& b, T sa, T sb)
    {
        std::vector<Entry> out;
        out.reserve(itsD.size()+b.itsD.size());
        size_t i=0, j=0;
        while (i<itsD.size() || j<b.itsD.size())
        {
            if (j==b.itsD.size() || (i<itsD.size() && itsD[i].atom<b.itsD[j].atom)) { out.push_back({itsD[i].atom, sa*itsD[i].d}); ++i; }
            else if (i==itsD.size() || b.itsD[j].atom<itsD[i].atom)            { out.push_back({b.itsD[j].atom, sb*b.itsD[j].d}); ++j; }
            else { out.push_back({itsD[i].atom, sa*itsD[i].d+sb*b.itsD[j].d}); ++i; ++j; }
        }
        itsD=std::move(out);
    }

    T                  itsValue{};
    std::vector<Entry> itsD;   //!< sorted by atom, no duplicates
};

} // namespace qchem
