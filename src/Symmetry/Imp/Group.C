// File: Symmetry/Imp/Group.C  The Group face's realizations, factories and intersection table (migration step 1).
module;
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdlib>
#include <memory>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
module qchem.Symmetry.Group;

namespace qchem::Symmetry
{
using Mat=Matrix3D<double>;

// ---------------------------------------------------------------------------------------------------------------------------
//  Generator
Generator Generator::Op(const SymOp& g)                     {Generator r(Kind::Op); r.itsOp=g; return r;}
Generator Generator::Translations(const Matrix3D<double>& A)
{
    if (std::fabs(Determinant(A))==0.0) throw std::invalid_argument("Generator::Translations: the lattice is not full rank");
    Generator r(Kind::Translations); r.itsLattice=A; return r;
}
Generator Generator::FullRotations(const rvec3_t& c)        {Generator r(Kind::FullRotations); r.itsCentre=c; return r;}
Generator Generator::AxialRotations(const rvec3_t& c, const rvec3_t& n)
{
    const double len=std::sqrt(n*n);
    if (len==0.0) throw std::invalid_argument("Generator::AxialRotations: a zero axis");
    Generator r(Kind::AxialRotations); r.itsCentre=c; r.itsAxis=n/len; return r;
}
Generator Generator::SpinRotations(const SpinRotationGroup& s) {Generator r(Kind::SpinRotations); r.itsSpin=s; return r;}

namespace
{
template <class T> const T& Typed(const std::optional<T>& slot, const char* what)
{
    if (!slot) throw std::logic_error(std::string("Generator: asked for its ")+what+" but that realization is not the populated one");
    return *slot;
}
} //anon
const SymOp&             Generator::GetOp       () const {return Typed(itsOp,      "op");}
const Matrix3D<double>&  Generator::GetLattice  () const {return Typed(itsLattice, "lattice");}
const rvec3_t&           Generator::GetCentre   () const {return Typed(itsCentre,  "centre");}
const rvec3_t&           Generator::GetAxis     () const {return Typed(itsAxis,    "axis");}
const SpinRotationGroup& Generator::GetSpinGroup() const {return Typed(itsSpin,    "spin group");}

// ---------------------------------------------------------------------------------------------------------------------------
//  Small geometry
namespace
{
double  Norm  (const rvec3_t& v) {return std::sqrt(v*v);}
rvec3_t Column(const Mat& A, size_t j) {return rvec3_t(A(1,j+1), A(2,j+1), A(3,j+1));}
Mat     FromColumns(const rvec3_t& a, const rvec3_t& b, const rvec3_t& c)
{
    return Mat(a.x, b.x, c.x,
               a.y, b.y, c.y,
               a.z, b.z, c.z);
}
Mat Identity() {return Mat(1,0,0, 0,1,0, 0,0,1);}
double MaxAbsDiff(const Mat& a, const Mat& b)
{
    double m=0.0;
    for (size_t i=1;i<=3;i++)
        for (size_t j=1;j<=3;j++) m=std::max(m, std::fabs(a(i,j)-b(i,j)));
    return m;
}
double MaxAbs(const rvec3_t& v) {return std::max({std::fabs(v.x), std::fabs(v.y), std::fabs(v.z)});}

GroupTolerance Looser(const GroupTolerance& a, const GroupTolerance& b) {return {std::max(a.length, b.length), std::max(a.angle, b.angle)};}

SymOp E() {return {Identity(), rvec3_t(0,0,0), SpinAction::None};}
SpinAction Times(SpinAction a, SpinAction b) {return a==b ? SpinAction::None : SpinAction::Flip;}
//! (a∘b)(r) = R_a (R_b r + t_b) + t_a.
SymOp Compose(const SymOp& a, const SymOp& b) {return {a.W*b.W, a.W*b.tau + a.tau, Times(a.sigma, b.sigma)};}
bool  SameOp (const SymOp& a, const SymOp& b, const GroupTolerance& tol)
{
    return a.sigma==b.sigma && MaxAbsDiff(a.W, b.W)<=tol.angle && Norm(a.tau-b.tau)<=tol.length;
}
//! The mirror through the plane containing \a c with unit normal \a m.
SymOp Mirror(const rvec3_t& c, const rvec3_t& m)
{
    Mat R=Identity();
    R-=2.0*Mat(m.x*m.x, m.x*m.y, m.x*m.z,
               m.y*m.x, m.y*m.y, m.y*m.z,
               m.z*m.x, m.z*m.y, m.z*m.z);
    return {R, c-R*c, SpinAction::None};
}
//! Any unit vector perpendicular to the unit \a d.
rvec3_t Perpendicular(const rvec3_t& d)
{
    const rvec3_t trial=std::fabs(d.x)<0.9 ? rvec3_t(1,0,0) : rvec3_t(0,1,0);
    const rvec3_t p=Cross(d, trial);
    return p/Norm(p);
}
std::string Str(const rvec3_t& v) {std::ostringstream s; s<<"("<<v.x<<", "<<v.y<<", "<<v.z<<")"; return s.str();}

// ---------------------------------------------------------------------------------------------------------------------------
//  Exact lattice algebra.  A lattice is handed over as Cartesian basis columns; two lattices are compared in one's FRACTIONAL
//  frame, where a commensurate partner has RATIONAL coordinates N/q.  Rationalising is the only tolerance-dependent step (the
//  Cartesian residual of N/q against the true coordinates must be within the length tolerance); everything after it is integer.

//! The largest common denominator searched when rationalising a lattice's coordinates in another's frame.  Commensurate cells in
//! use are small supercells (MnO AFM-II doubles along [111], q 2); 1000 leaves room for any supercell a run can afford.
constexpr long long kMaxDenominator=1000;
//! The most candidate elements \c Intersect enumerates for two crystals (cosets x the index of the common sublattice).
constexpr size_t    kMaxIntersectCandidates=size_t(1)<<20;

using IVec=std::array<long long,3>;
struct Rational { std::vector<IVec> N; long long q; };   // fractional coordinates N/q, one IVec per vector

//! Rationalise fractional vectors \a f (in \a A's frame) with ONE common denominator; nullopt if none up to kMaxDenominator fits.
std::optional<Rational> Rationalise(const std::vector<rvec3_t>& f, const Mat& A, double tolLength)
{
    for (long long q=1;q<=kMaxDenominator;q++)
    {
        Rational r{{}, q};
        bool ok=true;
        for (const rvec3_t& v:f)
        {
            const IVec n{std::llround(q*v.x), std::llround(q*v.y), std::llround(q*v.z)};
            const rvec3_t back(double(n[0])/q, double(n[1])/q, double(n[2])/q);
            if (Norm(A*(v-back))>tolLength) {ok=false; break;}
            r.N.push_back(n);
        }
        if (ok) return r;
    }
    return std::nullopt;
}

//! Integer column echelon of the 3 x m matrix whose COLUMNS are \a cols, by unimodular column operations, tracking the m x m
//! transform in \a U (its columns, each of length m).  On return column r (r < rank) has its pivot in row r and zeros above it
//! (lower-triangular), every later column is zero in rows <= r, and the pivots are positive.
void ColumnEchelon(std::vector<IVec>& cols, std::vector<std::vector<long long>>& U)
{
    const size_t m=cols.size();
    U.assign(m, std::vector<long long>(m, 0));
    for (size_t j=0;j<m;j++) U[j][j]=1;
    auto axpy=[&](size_t dst, size_t src, long long k)   // col[dst] -= k col[src]
    {
        for (int i=0;i<3;i++) cols[dst][i]-=k*cols[src][i];
        for (size_t i=0;i<m;i++) U[dst][i]-=k*U[src][i];
        for (int i=0;i<3;i++)
            if (std::llabs(cols[dst][i])>(1LL<<40)) throw std::overflow_error("Group lattice algebra: integer entries grew past 2^40");
    };
    size_t r=0;
    for (int row=0;row<3 && r<m;row++)
    {
        for (;;)
        {
            size_t best=m;
            for (size_t j=r;j<m;j++)
                if (cols[j][row]!=0 && (best==m || std::llabs(cols[j][row])<std::llabs(cols[best][row]))) best=j;
            if (best==m) break;                                           // the row is empty from r on
            std::swap(cols[r], cols[best]); std::swap(U[r], U[best]);
            bool done=true;
            for (size_t j=r+1;j<m;j++)
                if (cols[j][row]!=0) {axpy(j, r, cols[j][row]/cols[r][row]); if (cols[j][row]!=0) done=false;}
            if (done) break;
        }
        if (r<m && cols[r][row]!=0)
        {
            if (cols[r][row]<0) {for (auto& x:cols[r]) x=-x; for (auto& x:U[r]) x=-x;}
            r++;
        }
    }
}

//! T(L1) ∩ T(L2) = T(L1 ∩ L2), exactly.  In A1's frame scaled by q, L1 = qZ^3 and L2 = N Z^3; the common vectors x = q u = N v are
//! the kernel of [qI | -N], and the intersection's basis is A1 u over that kernel's basis.  nullopt when the lattices are
//! incommensurate at this tolerance.
std::optional<Mat> LatticeIntersection(const Mat& A1, const Mat& A2, double tolLength)
{
    const Mat M=Invert(A1)*A2;
    const auto rat=Rationalise({Column(M,0), Column(M,1), Column(M,2)}, A1, tolLength);
    if (!rat) return std::nullopt;
    std::vector<IVec> cols;
    for (int j=0;j<3;j++) {IVec c{0,0,0}; c[j]=rat->q; cols.push_back(c);}
    for (const IVec& n:rat->N) cols.push_back(IVec{-n[0], -n[1], -n[2]});
    std::vector<std::vector<long long>> U;
    ColumnEchelon(cols, U);
    for (size_t j=3;j<6;j++)
        if (cols[j]!=IVec{0,0,0}) throw std::logic_error("LatticeIntersection: [qI|-N] is not rank 3");
    rvec3_t b[3];
    for (int k=0;k<3;k++) b[k]=A1*rvec3_t(double(U[3+k][0]), double(U[3+k][1]), double(U[3+k][2]));
    return FromColumns(b[0], b[1], b[2]);
}

//! The lattice generated by \a A's columns and the extra Cartesian translations \a extra (a centring); throws if they are not
//! rational in A's frame (a translation of a crystal's own group always is).
Mat LatticeSum(const Mat& A, const std::vector<rvec3_t>& extra, double tolLength)
{
    if (extra.empty()) return A;
    const Mat Ai=Invert(A);
    std::vector<rvec3_t> f;
    for (const rvec3_t& t:extra) f.push_back(Ai*t);
    const auto rat=Rationalise(f, A, tolLength);
    if (!rat) throw std::logic_error("Group: a crystal's pure translation is irrational in its own lattice frame");
    std::vector<IVec> cols;
    for (int j=0;j<3;j++) {IVec c{0,0,0}; c[j]=rat->q; cols.push_back(c);}
    for (const IVec& n:rat->N) cols.push_back(n);
    std::vector<std::vector<long long>> U;
    ColumnEchelon(cols, U);
    rvec3_t b[3];
    for (int k=0;k<3;k++) b[k]=A*rvec3_t(double(cols[k][0])/rat->q, double(cols[k][1])/rat->q, double(cols[k][2])/rat->q);
    return FromColumns(b[0], b[1], b[2]);
}

//! Representatives of L2 / L for a sublattice L ⊆ L2, as Cartesian vectors of L2: with L's basis in L2's frame brought to the
//! lower-triangular echelon H, the integer n with 0 <= n_i < H_ii are a complete residue system (|det H| of them).
std::vector<rvec3_t> CosetRepresentatives(const Mat& A2, const Mat& L, double tolLength)
{
    const Mat X=Invert(A2)*L;
    std::vector<IVec> cols;
    for (int j=0;j<3;j++)
    {
        const rvec3_t x=Column(X,j);
        const IVec n{std::llround(x.x), std::llround(x.y), std::llround(x.z)};
        if (Norm(A2*(x-rvec3_t(double(n[0]), double(n[1]), double(n[2]))))>tolLength)
            throw std::logic_error("CosetRepresentatives: L is not a sublattice of L2");
        cols.push_back(n);
    }
    std::vector<std::vector<long long>> U;
    ColumnEchelon(cols, U);
    const long long h0=cols[0][0], h1=cols[1][1], h2=cols[2][2];
    std::vector<rvec3_t> reps;
    for (long long a=0;a<h0;a++)
        for (long long b=0;b<h1;b++)
            for (long long c=0;c<h2;c++) reps.push_back(A2*rvec3_t(double(a), double(b), double(c)));
    return reps;
}

//! Is the Cartesian \a v a vector of the lattice \a A?
bool InLattice(const rvec3_t& v, const Mat& A, const Mat& Ainv, double tolLength)
{
    const rvec3_t f=Ainv*v;
    const rvec3_t n(std::round(f.x), std::round(f.y), std::round(f.z));
    return Norm(A*(f-n))<=tolLength;
}
//! \a t reduced into A's unit cell [0,1)^3 (fractionally), a canonical coset representative.
rvec3_t ReduceIntoCell(const rvec3_t& t, const Mat& A, const Mat& Ainv, double tolLength)
{
    rvec3_t f=Ainv*t;
    f=rvec3_t(f.x-std::floor(f.x), f.y-std::floor(f.y), f.z-std::floor(f.z));
    rvec3_t r=A*f;
    // a coordinate within tolerance of 1 is the same coset as 0: snap it, so equal cosets get equal reps.
    const rvec3_t a[3]={Column(A,0), Column(A,1), Column(A,2)};
    const double fc[3]={f.x, f.y, f.z};
    for (int i=0;i<3;i++)
        if ((1.0-fc[i])*Norm(a[i])<=tolLength) r=r-a[i];
    return r;
}

// ---------------------------------------------------------------------------------------------------------------------------
//  The three realizations.  Not exported: every client holds a GroupPtr from a factory.

class FiniteRealization : public Group
{
public:
    FiniteRealization(std::vector<SymOp> generators, const GroupTolerance& tol, std::string name={}) : itsName(std::move(name)), itsTol(tol)
    {
        itsOps.push_back(E());
        auto add=[&](const SymOp& g) {for (const SymOp& h:itsOps) if (SameOp(g, h, itsTol)) return false; itsOps.push_back(g); return true;};
        for (const SymOp& g:generators) add(g);
        for (bool grew=true; grew;)
        {
            grew=false;
            const size_t n=itsOps.size();
            for (size_t i=0;i<n;i++)
                for (size_t j=0;j<n;j++)
                    if (add(Compose(itsOps[i], itsOps[j])))
                    {
                        grew=true;
                        if (itsOps.size()>kMaxFiniteOrder)
                            throw std::invalid_argument("FiniteGroup: the closure of the ops passed "+std::to_string(kMaxFiniteOrder)
                                                        +" elements: they are not a finite group at this tolerance");
                    }
        }
    }
    using Group::Contains;
    bool Contains(const Generator& g) const override
    {
        if (g.GetKind()!=Generator::Kind::Op) return false;   // a finite group contains no continuous family, no lattice
        for (const SymOp& h:itsOps) if (SameOp(g.GetOp(), h, itsTol)) return true;
        return false;
    }
    Space ActsOn() const override
    {
        for (const SymOp& h:itsOps) if (h.sigma==SpinAction::Flip) return Space::RSpin;
        return Space::R;
    }
    std::string Describe() const override
    {
        const std::string order="order "+std::to_string(itsOps.size());
        return itsName.empty() ? "a finite group of "+order : itsName+", "+order;
    }
protected:
    Form GetForm() const override {return {Form::Kind::Finite, itsOps, Mat(), rvec3_t(0,0,0), {}, itsTol};}
private:
    std::string        itsName;   //!< what the issuer says the group IS (e.g. "D2h (abelian subgroup of D6h)"); empty = unnamed
    std::vector<SymOp> itsOps;
    GroupTolerance     itsTol;
};

class CrystalRealization : public Group
{
public:
    //! \a ops Cartesian; their translations may include a centring (R = I, no flip), which joins the lattice.
    CrystalRealization(const Mat& A, const std::vector<SymOp>& ops, const GroupTolerance& tol) : itsA(A), itsTol(tol)
    {
        std::vector<rvec3_t> centring;
        Mat Ai=Invert(itsA);
        for (const SymOp& g:ops)
            if (g.sigma==SpinAction::None && MaxAbsDiff(g.W, Identity())<=tol.angle && !InLattice(g.tau, itsA, Ai, tol.length))
                centring.push_back(g.tau);
        itsA=LatticeSum(itsA, centring, tol.length);
        itsAinv=Invert(itsA);
        itsCosets.push_back(E());
        for (const SymOp& g:ops)
        {
            const SymOp r{g.W, ReduceIntoCell(g.tau, itsA, itsAinv, tol.length), g.sigma};
            if (!ContainsOp(r)) itsCosets.push_back(r);
        }
    }
    using Group::Contains;
    bool Contains(const Generator& g) const override
    {
        switch (g.GetKind())
        {
        case Generator::Kind::Op:           return ContainsOp(g.GetOp());
        case Generator::Kind::Translations:
        {
            const Mat& L=g.GetLattice();
            for (size_t j=0;j<3;j++) if (!InLattice(Column(L,j), itsA, itsAinv, itsTol.length)) return false;
            return true;
        }
        default: return false;              // a crystal holds no continuous rotation family, and no spin group
        }
    }
    Space ActsOn() const override
    {
        for (const SymOp& h:itsCosets) if (h.sigma==SpinAction::Flip) return Space::RSpin;
        return Space::R;
    }
    std::string Describe() const override
    {
        return "a crystal group with "+std::to_string(itsCosets.size())+" cosets of its lattice (det "
               +std::to_string(Determinant(itsA))+" bohr^3)"+(ActsOn()==Space::RSpin ? ", magnetic" : "");
    }
protected:
    Form GetForm() const override {return {Form::Kind::Crystal, itsCosets, itsA, rvec3_t(0,0,0), {}, itsTol};}
private:
    bool ContainsOp(const SymOp& g) const
    {
        for (const SymOp& c:itsCosets)
            if (c.sigma==g.sigma && MaxAbsDiff(c.W, g.W)<=itsTol.angle && InLattice(g.tau-c.tau, itsA, itsAinv, itsTol.length))
                return true;
        return false;
    }
    Mat                itsA, itsAinv;
    std::vector<SymOp> itsCosets;
    GroupTolerance     itsTol;
};

//! The POINTWISE stabiliser of an affine subspace: a point (O(3)_c) or a line (C∞v).  Every orthogonal op fixing it pointwise.
class StabiliserRealization : public Group
{
public:
    StabiliserRealization(const rvec3_t& c, std::vector<rvec3_t> dirs, const GroupTolerance& tol) : itsC(c), itsDirs(std::move(dirs)), itsTol(tol)
    {
        if (itsDirs.size()>1) throw std::logic_error("StabiliserRealization: a plane or the whole space is a finite stabiliser");
        for (rvec3_t& d:itsDirs) d=d/Norm(d);
    }
    using Group::Contains;
    bool Contains(const Generator& g) const override
    {
        switch (g.GetKind())
        {
        case Generator::Kind::Op:
        {
            const SymOp& op=g.GetOp();
            if (op.sigma!=SpinAction::None) return false;
            if (Norm(op.W*itsC+op.tau-itsC)>itsTol.length) return false;
            for (const rvec3_t& d:itsDirs) if (MaxAbs(op.W*d-d)>itsTol.angle) return false;
            return true;
        }
        case Generator::Kind::FullRotations:  return itsDirs.empty() && Norm(g.GetCentre()-itsC)<=itsTol.length;
        case Generator::Kind::AxialRotations:
        {
            // SO(2) about the line l' fixes our subspace pointwise iff the subspace lies on l'.
            const rvec3_t n=g.GetAxis(), dc=itsC-g.GetCentre();
            if (Norm(dc-(dc*n)*n)>itsTol.length) return false;
            for (const rvec3_t& d:itsDirs) if (Norm(Cross(d, n))>itsTol.angle) return false;
            return true;
        }
        default: return false;
        }
    }
    Space ActsOn() const override {return Space::R;}
    std::string Describe() const override
    {
        return itsDirs.empty() ? "O(3) about "+Str(itsC) : "C∞v about the line through "+Str(itsC)+" along "+Str(itsDirs[0]);
    }
protected:
    Form GetForm() const override {return {Form::Kind::Stabiliser, {}, Mat(), itsC, itsDirs, itsTol};}
private:
    rvec3_t              itsC;
    std::vector<rvec3_t> itsDirs;
    GroupTolerance       itsTol;
};

// ---------------------------------------------------------------------------------------------------------------------------
//  The intersection table.

//! Two pointwise stabilisers meet in the stabiliser of the affine hull of both subspaces: a point -> O(3), a line -> C∞v, a
//! plane -> {E, its mirror}, the whole space -> {E}.
GroupPtr IntersectStabilisers(const Group::Form& a, const Group::Form& b, const GroupTolerance& tol)
{
    std::vector<rvec3_t> basis;
    auto addDirection=[&](rvec3_t v, double cut)
    {
        for (const rvec3_t& e:basis) v=v-(v*e)*e;
        const double len=Norm(v);
        if (len>cut) basis.push_back(v/len);
    };
    for (const rvec3_t& d:a.directions) addDirection(d, tol.angle);
    for (const rvec3_t& d:b.directions) addDirection(d, tol.angle);
    {
        // the offset between the two subspaces, less its components along their directions, measured as a length
        rvec3_t v=b.centre-a.centre;
        for (const rvec3_t& e:basis) v=v-(v*e)*e;
        const double len=Norm(v);
        if (len>tol.length) basis.push_back(v/len);
    }
    switch (basis.size())
    {
    case 0:  return std::make_shared<StabiliserRealization>(a.centre, std::vector<rvec3_t>{}, tol);
    case 1:  return std::make_shared<StabiliserRealization>(a.centre, std::vector<rvec3_t>{basis[0]}, tol);
    case 2:
    {
        const rvec3_t m=Cross(basis[0], basis[1]);
        return std::make_shared<FiniteRealization>(std::vector<SymOp>{Mirror(a.centre, m/Norm(m))}, tol);
    }
    default: return std::make_shared<FiniteRealization>(std::vector<SymOp>{}, tol);
    }
}

//! The elements of the crystal \a cr that fix the point \a c (its site stabiliser, SpaceGroup::SiteStabilizer's answer): a coset
//! {R|t} holds one iff c - R c - t is a lattice vector, and the element is then {R | c - R c}.
std::vector<SymOp> SiteStabiliser(const Group::Form& cr, const rvec3_t& c, const GroupTolerance& tol)
{
    const Mat Ai=Invert(cr.lattice);
    std::vector<SymOp> out;
    for (const SymOp& g:cr.ops)
        if (InLattice(c-g.W*c-g.tau, cr.lattice, Ai, tol.length)) out.push_back({g.W, c-g.W*c, g.sigma});
    return out;
}

} //anon

// ---------------------------------------------------------------------------------------------------------------------------
//  Group

std::vector<Generator> Group::Generators() const
{
    const Form f=GetForm();
    std::vector<Generator> gens;
    switch (f.kind)
    {
    case Form::Kind::Finite:     for (const SymOp& g:f.ops) gens.push_back(Generator::Op(g)); break;
    case Form::Kind::Crystal:    gens.push_back(Generator::Translations(f.lattice));
                                 for (const SymOp& g:f.ops) gens.push_back(Generator::Op(g)); break;
    case Form::Kind::Stabiliser:
        if (f.directions.empty()) gens.push_back(Generator::FullRotations(f.centre));
        else
        {
            gens.push_back(Generator::AxialRotations(f.centre, f.directions[0]));
            gens.push_back(Generator::Op(Mirror(f.centre, Perpendicular(f.directions[0]))));   // a mirror containing the axis
        }
        break;
    }
    return gens;
}

bool Group::HasLattice() const {return GetForm().kind==Form::Kind::Crystal;}
Frame Group::MakeFrame() const
{
    const Form f=GetForm();
    Frame out;
    if (f.kind==Form::Kind::Crystal) out.lattice=f.lattice;
    if (f.kind==Form::Kind::Stabiliser && f.directions.empty()) out.centre=f.centre;   // a point (O(3)_c), not a line
    return out;
}

bool Group::Contains(const Group& h) const
{
    for (const Generator& g:h.Generators()) if (!Contains(g)) return false;
    return true;
}

Outcome<GroupPtr,std::string> Group::Intersect(const Group& o) const
{
    using R=Outcome<GroupPtr,std::string>;
    const Form a=GetForm(), b=o.GetForm();
    const GroupTolerance tol=Looser(a.tol, b.tol);
    auto filter=[&](const std::vector<SymOp>& candidates, const Group& by)
    {
        std::vector<SymOp> kept;
        for (const SymOp& g:candidates) if (by.Contains(Generator::Op(g))) kept.push_back(g);
        return kept;
    };
    // A finite side enumerates; the other filters.
    if (a.kind==Form::Kind::Finite) return R::Ok(std::make_shared<FiniteRealization>(filter(a.ops, o), tol));
    if (b.kind==Form::Kind::Finite) return R::Ok(std::make_shared<FiniteRealization>(filter(b.ops, *this), tol));
    // Two continuous stabilisers: the closed table.
    if (a.kind==Form::Kind::Stabiliser && b.kind==Form::Kind::Stabiliser) return R::Ok(IntersectStabilisers(a, b, tol));
    // A crystal and a stabiliser: the crystal's site stabiliser at the stabiliser's point (finite), filtered by the stabiliser.
    if (a.kind==Form::Kind::Crystal && b.kind==Form::Kind::Stabiliser)
        return R::Ok(std::make_shared<FiniteRealization>(filter(SiteStabiliser(a, b.centre, tol), o), tol));
    if (b.kind==Form::Kind::Crystal && a.kind==Form::Kind::Stabiliser)
        return R::Ok(std::make_shared<FiniteRealization>(filter(SiteStabiliser(b, a.centre, tol), *this), tol));
    // Two crystals: the common sublattice L = L_a ∩ L_b (exact), then b's elements modulo L, filtered by a.
    const auto L=LatticeIntersection(a.lattice, b.lattice, tol.length);
    if (!L) return R::Fail("the two crystals' lattices have no common sublattice with denominators up to "
                           +std::to_string(kMaxDenominator)+" (incommensurate at length tolerance "+std::to_string(tol.length)
                           +" bohr): "+Describe()+" ∩ "+o.Describe());
    const std::vector<rvec3_t> reps=CosetRepresentatives(b.lattice, *L, tol.length);
    if (reps.size()*b.ops.size()>kMaxIntersectCandidates)
        return R::Fail("the intersection would enumerate "+std::to_string(reps.size()*b.ops.size())+" candidates (more than "
                       +std::to_string(kMaxIntersectCandidates)+"): "+Describe()+" ∩ "+o.Describe());
    std::vector<SymOp> candidates;
    for (const SymOp& g:b.ops)
        for (const rvec3_t& n:reps) candidates.push_back({g.W, g.tau+n, g.sigma});
    return R::Ok(std::make_shared<CrystalRealization>(*L, filter(candidates, *this), tol));
}

// ---------------------------------------------------------------------------------------------------------------------------
//  Factories

GroupPtr PointRotationGroup(const rvec3_t& c, const GroupTolerance& tol)
{
    return std::make_shared<StabiliserRealization>(c, std::vector<rvec3_t>{}, tol);
}
GroupPtr LineRotationGroup(const rvec3_t& c, const rvec3_t& n, const GroupTolerance& tol)
{
    if (Norm(n)==0.0) throw std::invalid_argument("LineRotationGroup: a zero axis");
    return std::make_shared<StabiliserRealization>(c, std::vector<rvec3_t>{n}, tol);
}
GroupPtr TranslationGroup(const Matrix3D<double>& A, const GroupTolerance& tol)
{
    return std::make_shared<CrystalRealization>(A, std::vector<SymOp>{}, tol);
}
GroupPtr TrivialGroup()
{
    return std::make_shared<FiniteRealization>(std::vector<SymOp>{}, kExactTolerance, "{e}");
}
GroupPtr FiniteGroup(const std::vector<SymOp>& cartesianOps, const GroupTolerance& tol)
{
    return std::make_shared<FiniteRealization>(cartesianOps, tol);
}
GroupPtr FiniteGroup(const std::vector<Molecule::SymOp>& ops, const rvec3_t& origin, const GroupTolerance& tol, const std::string& name)
{
    std::vector<SymOp> cart;
    for (const Molecule::SymOp& op:ops) cart.push_back({op.Matrix(), origin-op.Matrix()*origin, SpinAction::None});
    return std::make_shared<FiniteRealization>(cart, tol, name);
}
namespace
{
std::vector<SymOp> ToCartesian(const std::vector<SymOp>& fractional, const Mat& A)
{
    const Mat Ai=Invert(A);
    std::vector<SymOp> cart;
    for (const SymOp& g:fractional) cart.push_back({A*g.W*Ai, A*g.tau, g.sigma});
    return cart;
}
} //anon
GroupPtr CrystalGroup(const Lattice_3D::SpaceGroup& sg, const GroupTolerance& tol)
{
    std::vector<SymOp> frac;
    for (const Lattice_3D::SpaceGroupOp& op:sg.Ops()) frac.push_back({op.W, op.tau, SpinAction::None});
    return std::make_shared<CrystalRealization>(sg.CellMatrix(), ToCartesian(frac, sg.CellMatrix()), tol);
}
GroupPtr CrystalGroup(const std::vector<SymOp>& fractionalOps, const Matrix3D<double>& A, const GroupTolerance& tol)
{
    return std::make_shared<CrystalRealization>(A, ToCartesian(fractionalOps, A), tol);
}

} //namespace
