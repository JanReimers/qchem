// File: Structure/Lattice_3D/Imp/Lattice_3D.C Define a 3D infinite lattice.
module;
#include <iostream>
#include <cassert>
#include <algorithm> //sort
#include <cmath>     //round, abs (SiteRotations), pow (SiteEnvironmentRotations)
#include <stdexcept>
#include <vector>
#include <memory>    //make_shared (GetStructure)

module qchem.Lattice_3D;

import qchem.Streamable;
import qchem.Math;
import qchem.Blaze; //for op!= on blaze iterators (range-for over rvec3vec_t)

namespace qchem {

//--------------------------------------------------------------------------
//
//  Construction zone.
//
// Lattice_3D::Lattice(            )
//     : itsUnitCell (            )
//     , itsLimits   (0,0,0       )
//     , itsAtoms    (new Molecule)
//     , itsTolerence(0.0001      )
// {}

Lattice_3D::Lattice_3D(const UnitCell& cell, const Vector3D<int>& Limits)
    : itsUnitCell (cell        )
    , itsLimits   (Limits      )
    , itsTolerence(0.0001      )
{}


//---------------------------------------------------------
//
//  Structure stuff.
//

ReciprocalLattice Lattice_3D::Reciprocal() const
{
    return ReciprocalLattice(itsUnitCell.MakeReciprocalCell());
}

std::shared_ptr<const Structure> Lattice_3D::GetStructure() const
{
    return std::make_shared<UnitCell>(itsUnitCell); // deep-copies the atom basis (Cartesian a.u.)
}

const Symmetry::Lattice_3D::SpaceGroup& Lattice_3D::GetSpaceGroup(double tol) const
{
    namespace SL = Symmetry::Lattice_3D;
    if (!itsSpaceGroup)
    {
        // The UnitCell -> {A, fractional sites} adapter (formerly private to the GPW factory).
        const Structure& st = itsUnitCell;   // iterate atoms via the PUBLIC Structure interface
        std::vector<SL::AtomSite> sites;
        for (Atom* a : st)
            sites.push_back({a->itsZ, itsUnitCell.ToFractional(a->itsR)});
        itsSpaceGroup = std::make_shared<const SL::SpaceGroup>(
                            SL::SpaceGroup::Detect(itsUnitCell.GetCellMatrix(), sites, tol));
    }
    return *itsSpaceGroup;
}

std::vector<Symmetry::Lattice_3D::SymOp> Lattice_3D::ShubnikovOps(const std::vector<int>& spins,
                                                                  double tol) const
{
    namespace SL = Symmetry::Lattice_3D;
    const Structure& st = itsUnitCell;
    assert(spins.size()==st.GetNumAtoms() && "ShubnikovOps: one collinear label per atom, in atom order");
    std::vector<SL::AtomSite> decorated;
    size_t i=0;
    for (Atom* a : st)
        decorated.push_back({a->itsZ, itsUnitCell.ToFractional(a->itsR), spins[i++]});
    return GetSpaceGroup(tol).ShubnikovOps(decorated, tol);
}

std::vector<rmat3d_t> Lattice_3D::SiteRotations(size_t atom, const std::vector<int>& spins, double tol) const
{
    namespace SL = Symmetry::Lattice_3D;
    const Structure& st = itsUnitCell;
    assert(atom<st.GetNumAtoms() && "SiteRotations: atom index outside the cell");
    std::vector<int> sp = spins.empty() ? std::vector<int>(st.GetNumAtoms(),0) : spins;
    // The site's fractional position, in cell atom order.
    rvec3_t f; { size_t i=0; for (Atom* a : st) { if (i==atom) f=itsUnitCell.ToFractional(a->itsR); i++; } }
    const Matrix3D<double>& A=itsUnitCell.GetCellMatrix();
    const Matrix3D<double>  Ainv=Invert(A);
    std::vector<rmat3d_t> R;
    for (const SL::SymOp& op : ShubnikovOps(sp, tol))
    {
        if (op.sigma!=Symmetry::SpinAction::None) continue;
        const rvec3_t d = op.W*f + op.tau - f;
        auto integral=[tol](double x){ return std::abs(x-std::round(x))<=tol; };
        if (!(integral(d.x) && integral(d.y) && integral(d.z))) continue;
        R.push_back(A*op.W*Ainv);
    }
    assert(!R.empty() && "SiteRotations: the identity fixes every site");
    return R;
}

std::vector<rmat3d_t> Lattice_3D::SiteEnvironmentRotations(size_t atom, size_t shells, double tol) const
{
    const Structure& st = itsUnitCell;
    assert(atom<st.GetNumAtoms() && "SiteEnvironmentRotations: atom index outside the cell");
    assert(shells>=1);
    // The site and the cell's atoms, Cartesian.
    struct Pt { int Z; rvec3_t r; };
    std::vector<Pt> cell; rvec3_t site;
    { size_t i=0; for (Atom* a : st) { cell.push_back({a->itsZ, a->itsR}); if (i==atom) site=a->itsR; i++; } }
    const Matrix3D<double>& A=itsUnitCell.GetCellMatrix();
    // Every image within a generous radius -- the images of the (shells) nearest distances are certainly
    // inside 3 cells in each direction for any sane cell; the shell cut below is what actually selects.
    std::vector<Pt> env;
    for (int i=-3;i<=3;i++) for (int j=-3;j<=3;j++) for (int k=-3;k<=3;k++)
    {
        const rvec3_t t=A*rvec3_t(i,j,k);
        for (const Pt& p : cell)
        {
            const rvec3_t d=p.r+t-site;
            if (norm(d)>tol) env.push_back({p.Z, d});
        }
    }
    // The first (shells) distinct distances.
    std::vector<double> dist; for (const Pt& p : env) dist.push_back(norm(p.r));
    std::sort(dist.begin(), dist.end());
    std::vector<double> distinct;
    for (double d : dist) if (distinct.empty() || d-distinct.back()>tol*std::max(1.0,d)) distinct.push_back(d);
    const double rcut = (shells<=distinct.size() ? distinct[shells-1] : distinct.back()) * (1.0+tol);
    std::vector<Pt> cl; for (const Pt& p : env) if (norm(p.r)<=rcut) cl.push_back(p);
    std::vector<Pt> shell1; for (const Pt& p : cl) if (norm(p.r)<=distinct[0]*(1.0+tol)) shell1.push_back(p);
    // A non-coplanar reference triple from the first shell (fall back to the whole cluster if the first
    // shell is planar or too small -- a linear/planar coordination).
    auto triple=[&](const std::vector<Pt>& pts, size_t& a, size_t& b, size_t& c)
    {
        for (a=0;a<pts.size();a++) for (b=a+1;b<pts.size();b++) for (c=b+1;c<pts.size();c++)
        {
            Matrix3D<double> M(pts[a].r.x,pts[b].r.x,pts[c].r.x, pts[a].r.y,pts[b].r.y,pts[c].r.y, pts[a].r.z,pts[b].r.z,pts[c].r.z);
            if (std::abs(Determinant(M))>1e-6*std::pow(norm(pts[a].r),3)) return true;
        }
        return false;
    };
    size_t ia,ib,ic;
    const std::vector<Pt>& ref = triple(shell1, ia, ib, ic) ? shell1 : cl;
    if (&ref==&cl && !triple(cl, ia, ib, ic))
        throw std::runtime_error("SiteEnvironmentRotations: the site's environment is coplanar -- no finite point group from it");
    const Matrix3D<double> Am(ref[ia].r.x,ref[ib].r.x,ref[ic].r.x, ref[ia].r.y,ref[ib].r.y,ref[ic].r.y, ref[ia].r.z,ref[ib].r.z,ref[ic].r.z);
    const Matrix3D<double> Ainv=Invert(Am);
    auto isSymmetry=[&](const rmat3d_t& R)
    {
        for (const Pt& p : cl)
        {
            const rvec3_t q=R*p.r; bool hit=false;
            for (const Pt& o : cl) if (o.Z==p.Z && norm(o.r-q)<=tol*std::max(1.0,norm(q))) { hit=true; break; }
            if (!hit) return false;
        }
        return true;
    };
    std::vector<rmat3d_t> R;
    auto seen=[&](const rmat3d_t& X)
    {
        for (const rmat3d_t& Y : R)
        {
            bool same=true;
            for (int i=1;i<=3 && same;i++) for (int j=1;j<=3;j++) if (std::abs(X(i,j)-Y(i,j))>1e-8) { same=false; break; }
            if (same) return true;
        }
        return false;
    };
    for (size_t a=0;a<ref.size();a++) if (ref[a].Z==ref[ia].Z)
    for (size_t b=0;b<ref.size();b++) if (b!=a && ref[b].Z==ref[ib].Z)
    for (size_t c=0;c<ref.size();c++) if (c!=a && c!=b && ref[c].Z==ref[ic].Z)
    {
        const Matrix3D<double> Bm(ref[a].r.x,ref[b].r.x,ref[c].r.x, ref[a].r.y,ref[b].r.y,ref[c].r.y, ref[a].r.z,ref[b].r.z,ref[c].r.z);
        const rmat3d_t X=Bm*Ainv;
        // orthogonal?
        const rmat3d_t G=Transpose(X)*X; bool orth=true;
        for (int i=1;i<=3 && orth;i++) for (int j=1;j<=3;j++) if (std::abs(G(i,j)-(i==j ? 1.0 : 0.0))>1e-6) { orth=false; break; }
        if (!orth || seen(X) || !isSymmetry(X)) continue;
        R.push_back(X);
    }
    assert(!R.empty() && "SiteEnvironmentRotations: the identity fixes every environment");
    return R;
}


//----------------------------------------------------------
//
//  Simple lattice questions.
//
size_t Lattice_3D::GetNumSites() const
{
    return GetNumBasisSites() * GetNumUnitCells();
}

size_t Lattice_3D::GetNumBasisSites() const
{
    return itsUnitCell.GetNumAtoms();
}

size_t Lattice_3D::GetNumUnitCells() const
{
    return itsLimits.x * itsLimits.y * itsLimits.z;
}

//------------------------------------------------------------------
//
//  Coordinate to site number translations.
//
size_t Lattice_3D::GetSiteNumber(const rvec3_t& r) const
{
    rvec3_t basis;
    Vector3D<int> cell;
    SplitCoordinate(r,basis,cell);

    size_t ib = Find(basis); //Find within tolerence.
    assert (ib<GetNumBasisSites());

    size_t sitenum=ib + GetNumBasisSites()*(cell.z + itsLimits.z*(cell.y + itsLimits.y*cell.x));
    assert(sitenum<GetNumSites());
    // std::cout << "GetCoordinate(sitenum)="  << GetCoordinate(sitenum) << std::endl;
    // std::cout << "rvec3_t(cell.x,cell.y,cell.z)="  << rvec3_t(cell.x,cell.y,cell.z) << std::endl;
    // std::cout << "basis="  << basis << std::endl;
    assert(itsUnitCell.GetDistance(GetCoordinate(sitenum)-rvec3_t(cell.x,cell.y,cell.z)-basis) < itsTolerence);
    return sitenum;
}

size_t Lattice_3D::GetBasisNumber(const rvec3_t& r) const
{
    rvec3_t basis;
    Vector3D<int> cell;
    SplitCoordinate(r,basis,cell);

    size_t ret=Find(basis);
    assert(ret<GetNumBasisSites());
    return ret;
}

size_t Lattice_3D::GetBasisNumber(size_t SiteNumber) const
{
    assert(SiteNumber<GetNumSites());
    return SiteNumber%GetNumBasisSites();
}

void Lattice_3D::SplitCoordinate(const rvec3_t& r, rvec3_t& basis, Vector3D<int>& cell) const
{
    cell.x=(int)floor(r.x);
    cell.y=(int)floor(r.y);
    cell.z=(int)floor(r.z);

    basis.x=r.x-cell.x;
    basis.y=r.y-cell.y;
    basis.z=r.z-cell.z;

    cell.x = cell.x%itsLimits.x;
    cell.y = cell.y%itsLimits.y;
    cell.z = cell.z%itsLimits.z;

    if (cell.x < 0) cell.x+=itsLimits.x;
    if (cell.y < 0) cell.y+=itsLimits.y;
    if (cell.z < 0) cell.z+=itsLimits.z;
}

Vector3D<int> Lattice_3D::GetCellCoord (const rvec3_t& r) const
{
    rvec3_t basis;
    Vector3D<int> ret;
    SplitCoordinate(r,basis,ret);
    return ret;
}


rvec3_t Lattice_3D::GetCoordinate(size_t SiteNumber) const
{
    assert(SiteNumber<GetNumSites());
    size_t ib=GetBasisNumber(SiteNumber);
    rvec3_t ret=GetBasisVector(ib);

    SiteNumber-=ib;
    SiteNumber/=GetNumBasisSites();
    int iz=SiteNumber%itsLimits.z;
    assert(iz>=0);
    assert(iz<itsLimits.z);

    SiteNumber-=iz;
    SiteNumber/=itsLimits.z;
    int iy=SiteNumber%itsLimits.y;
    assert(iy>=0);
    assert(iy<itsLimits.y);

    SiteNumber-=iy;
    SiteNumber/=itsLimits.y;
    int ix=SiteNumber;
    assert(ix>=0);
    assert(ix<itsLimits.x);

    ret+=rvec3_t(ix,iy,iz);
    return ret;
}

//----------------------------------------------------------
//
//  Advanced lattice questions.
//
rvec_t Lattice_3D::GetDistances(size_t NumShells) const
{
    double maxd=itsUnitCell.GetMinimumCellEdge()*NumShells; //Initial guess.
    std::vector<double> distances; //scratch buffer: built by push_back, then sorted

    rvec3vec_t super_cells=GetSuperCells(maxd);

    for (auto a1:itsUnitCell)
        for (auto a2:itsUnitCell)
            for (auto& c:super_cells)
            {
                double d=itsUnitCell.GetDistance(c + a2->itsR - a1->itsR);
                if(d>0 && d<=maxd && Find(d,distances)==distances.size()) distances.push_back(d);
            }

    std::sort(distances.begin(),distances.end());
    assert(distances.size()>=NumShells); //guess radius minEdge*NumShells should be ample
    return rvec_t(std::min(NumShells,distances.size()),distances.data());
}

rvec3vec_t Lattice_3D::GetBonds(size_t BasisNumber, double Distance) const
{
    assert(BasisNumber<GetNumBasisSites());
    assert(Distance>0);

    std::vector<rvec3_t> ret; //scratch buffer: size not known up front
    rvec3_t rb=GetBasisVector(BasisNumber);
    rvec3vec_t super_cells=GetSuperCells(Distance);

    for (auto a:itsUnitCell)
        for (auto& c:super_cells)
        {
            rvec3_t bond = a->itsR + c - rb;
            double mbond=itsUnitCell.GetDistance(bond);
            if (fabs(mbond-Distance) < itsTolerence) ret.push_back(bond);
        }
    return rvec3vec_t(ret.size(),ret.data());
}

rvec3vec_t Lattice_3D::GetBondsInSphere(size_t BasisNumber, double Distance) const
{
    assert(BasisNumber<GetNumBasisSites());
    assert(Distance>0);

    std::vector<rvec3_t> ret; //scratch buffer: size not known up front
    rvec3_t rb=GetBasisVector(BasisNumber);
    rvec3vec_t super_cells=GetSuperCells(Distance);

    for (auto a:itsUnitCell)
        for (auto& c:super_cells)
        {
            rvec3_t bond = a->itsR + c - rb;
            double mbond=itsUnitCell.GetDistance(bond);
            if (mbond<Distance+itsTolerence) ret.push_back(bond);
        }
    return rvec3vec_t(ret.size(),ret.data());
}

std::vector<ivec3_t>  Lattice_3D::GetCellsInSphere(double rmax) const
{
    return itsUnitCell.CellsInSphere(rmax); //direct lattice vectors R within rmax
}
//--------------------------------------------------------
//
//  Private unitilities.
//
size_t  Lattice_3D::Find(const rvec3_t& r) const //Search within the primary unit cell.
{
    size_t ret=GetNumBasisSites();
    size_t i=0;
    for (auto a:itsUnitCell)
    {
        if (itsUnitCell.GetDistance(r - a->itsR) < itsTolerence)
        {
            ret=i;
            break;
        }
        i++;
    } 
    return ret;
}

size_t  Lattice_3D::Find(double r,const std::vector<double>& lis) const
{
    size_t ret=lis.size();
    size_t i=0;
    for (std::vector<double>::const_iterator b(lis.begin()); b!=lis.end(); b++,i++) if (fabs(r-*b) < itsTolerence)
        {
            ret=i;
            break;
        }
    return ret;
}

rvec3vec_t Lattice_3D::GetSuperCells(double MaxDistance) const
{
    Vector3D<int> nc=itsUnitCell.GetNumCells(MaxDistance);
    rvec3vec_t ret((2*nc.x+1)*(2*nc.y+1)*(2*nc.z+1)); //full box, size known up front
    size_t k=0;
    for (int ix=-nc.x; ix<=nc.x; ix++)
        for (int iy=-nc.y; iy<=nc.y; iy++)
            for (int iz=-nc.z; iz<=nc.z; iz++)
                ret[k++]=ivec3_t(ix,iy,iz);
    return ret;
}

rvec3_t Lattice_3D::GetBasisVector(size_t BasisNumber) const
{
    assert(BasisNumber<GetNumBasisSites());
    rvec3_t ret;
    {
        for (auto b:itsUnitCell)
        {
            if (BasisNumber==0)
            {
                ret=b->itsR;
                break;
            }
            BasisNumber--;
        }
    }
    return ret;
}

using std::endl;
//------------------------------------------------------------
//
//  Streamable stuff.
//
std::ostream& Lattice_3D::Write(std::ostream& os) const
{
    os << itsUnitCell << endl << itsLimits;

    return os;
}





} // namespace qchem