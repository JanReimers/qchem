// File: Mesh/Imp/Recipe.C  The bodies of qchem.Mesh.Recipe: the pin-30 descriptor of a stated point set and the uniform-mesh resolution bridge.
module;
#include <cassert>
#include <string>
module qchem.Mesh.Recipe;
import qchem.Symmetry.Group;

namespace qchem::qcMesh
{

Outcome<Symmetry::GroupPtr,std::string> Forces(const MeshParams& mp, const Symmetry::Frame& frame)
{
    using R=Outcome<Symmetry::GroupPtr,std::string>;
    if (mp.cellKind==UnitCellKind::Uniform)
    {
        if (!frame.lattice) return R::Fail("a uniform cell mesh imposes H_R = T, and there is no lattice here to build it on -- To fix this, please state the XC "
                                           "mesh cellKind Becke or Auto, the atom-centred mesh");
        return R::Ok(Symmetry::TranslationGroup(*frame.lattice, Symmetry::kExactTolerance));
    }
    if (mp.angular==AngularKind::None)
    {
        if (!frame.centre) return R::Fail("the radial-only point set (angular none: one direction at 4 pi, exact only for a spherical integrand) "
                                          "imposes H_R = O(3) about one centre, and there is no such centre here -- To fix this, please state an angular rule "
                                          "(Lebedev or GaussLegendre)");
        return R::Ok(Symmetry::PointRotationGroup(*frame.centre, Symmetry::kExactTolerance));
    }
    return R::Ok(Symmetry::TrivialGroup());
}

int UniformDivisions(double cellEdge, double eCut)
{
    assert(cellEdge>0.0 && "UniformDivisions: a cell edge is a positive length");
    assert(eCut>0.0 && "UniformDivisions: a cutoff of 0 does not select a resolution (use nUniform)");
    return max(1, int(ceil(2.0*cellEdge*sqrt(2.0*eCut)/Pi)));
}

double UniformCutoff(double cellEdge, int n)
{
    assert(cellEdge>0.0 && "UniformCutoff: a cell edge is a positive length");
    assert(n>0 && "UniformCutoff: a mesh has at least one division per axis");
    const double G=n*Pi/(2.0*cellEdge);   // invert n = 2 a sqrt(2E)/pi for sqrt(2E)
    return 0.5*G*G;
}

} //namespace qchem::qcMesh
