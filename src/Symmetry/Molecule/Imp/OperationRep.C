// File: Symmetry/Imp/OperationRep.C  Whole-basis operation rep: center permutation + per-shell angular rep.
module;
#include <vector>
#include <algorithm>
#include <cassert>
module qchem.Symmetry.Molecule.OperationRep;

namespace qchem::Symmetry::Molecule
{

rmat_t BuildOperationRep(const std::vector<AoShell>& shells, const rmat3d_t& R,
                         const rvec3_t& origin, double tol)
{
    size_t nAO = 0;
    for (const auto& s : shells) nAO = std::max(nAO, s.offset + s.nComponents());
    rmat_t M(nAO, nAO, 0.0);

    double tol2 = tol*tol;
    for (const auto& B : shells)
    {
        rvec3_t img = origin + R*(B.center - origin);     // where this shell's center maps to
        const AoShell* Bp = nullptr;                      // the image shell (same type, at img)
        for (const auto& C : shells)
            if (C.shellType==B.shellType)
            {
                rvec3_t d = C.center - img;
                if (d*d <= tol2) { Bp = &C; break; }
            }
        assert(Bp && "BuildOperationRep: no image shell -- R is not a symmetry of the basis");

        // The shell's own angular rep D(b,a) -- Cartesian or spherical, the builder does not distinguish.
        rmat_t D = B.rep->Rep(R);
        size_t nc = B.nComponents();
        for (size_t a=0;a<nc;a++) for (size_t b=0;b<nc;b++)
            // normalized rep: (N^B_a / N^B'_b) D(b,a), placed at (image row, source col)
            M(Bp->offset+b, B.offset+a) = (B.norm[a]/Bp->norm[b]) * D(b,a);
    }
    return M;
}

} //namespace
