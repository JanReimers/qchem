// File: BasisSet/Molecule/tests/M_MEvaluator.C
//
// Guard for the isM_1E_Evaluator (matrix-delivery) evaluator category and the Orbital_1E_IBS dispatch.
// A Matrix1E_Adapter wraps any scalar is1E_Evaluator and pre-assembles its 1E matrices -- a stand-in for a
// block-oriented integral library / disk source (the eventual libcint wrapper).  We check (a) it satisfies
// isM_1E_Evaluator, (b) the matrices it delivers equal the ones the scalar IBS builds element-by-element,
// and (c) the 1E mixin's matrix-forwarding branch compiles (explicit instantiation below).
#include "gtest/gtest.h"
#include <ostream>
#include <string>
#include <vector>

import qchem.BasisSet.Molecule.Evaluators;                    // Evaluator, is1E_/isM_1E_Evaluator
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD;        // NR_Evaluator (a scalar 1E evaluator)
import qchem.BasisSet.Molecule.PG_Cart;                       // Orbital_IBS
import qchem.BasisSet.Molecule.IBS;                           // Orbital_1E_IBS<E> mixin (the dispatch)
import qchem.BasisSet.Orbital_1E_IBS;                         // cached Overlap()/Kinetic()/Nuclear()
import qchem.Cluster;
import qchem.Types;
import qchem.Blaze;

using namespace BasisSet::Molecule::Evaluators;
using BasisSet::Molecule::PG_Cart::Orbital_IBS;

// A matrix-delivery 1E evaluator: wraps a scalar is1E_Evaluator and assembles its 1E matrices on demand.
// (A real isM_1E source -- libcint, a disk cache -- would produce these without an underlying kernel
// evaluator; here we build on the scalar one so the test has a known-correct reference to compare to.)
template <is1E_Evaluator E>
class Matrix1E_Adapter : public virtual Evaluator
{
    const E& its;
    template <class K> rsmat_t Build(K k) const
    {
        rsmat_t S(its.size());
        for (auto i:its.indices()) for (auto j:its.indices(i)) S(i,j)=k(i,j);
        return S;
    }
public:
    Matrix1E_Adapter(const E& e) : its(e) {}
    virtual size_t        size() const {return its.size();}
    virtual rvec_t        Norm() const {return its.Norm();}
    virtual std::string   Name() const {return "Matrix1E_Adapter";}
    virtual std::ostream& Write(std::ostream& os) const {return os << "Matrix1E_Adapter[" << size() << "]";}

    rsmat_t OverlapMatrix()                  const {return Build([&](size_t i,size_t j){return its.Overlap(i,j);});}
    rsmat_t KineticMatrix()                  const {return Build([&](size_t i,size_t j){return its.Grad2  (i,j);});}
    rsmat_t NuclearMatrix(const Cluster* cl) const {return Build([&](size_t i,size_t j){return its.Nuclear(i,j,cl);});}
};

static_assert( isM_1E_Evaluator<Matrix1E_Adapter<PG_Cart_MnD::NR_Evaluator>>,
              "Matrix1E_Adapter must satisfy isM_1E_Evaluator");
static_assert(!is1E_Evaluator <Matrix1E_Adapter<PG_Cart_MnD::NR_Evaluator>>,
              "a pure matrix evaluator does NOT provide the per-element scalar kernels");

// Compile-check the mixin's matrix-forwarding branch (the scalar branch is already instantiated by PG_Cart).
template class BasisSet::Molecule::Orbital_1E_IBS<Matrix1E_Adapter<PG_Cart_MnD::NR_Evaluator>>;

TEST(M_MEvaluator, matrix_matches_scalar)
{
    Molecule h2o;
    h2o.Insert(new Atom(8, 0, rvec3_t(0, 0,      0.117)));   // same geometry as M_Evaluator
    h2o.Insert(new Atom(1, 0, rvec3_t(0, 0.757, -0.467)));
    h2o.Insert(new Atom(1, 0, rvec3_t(0,-0.757, -0.467)));

    rvec_t exps{1.0, 0.25};
    Orbital_IBS ibs(exps, 1, &h2o);                          // scalar PG evaluator (s + p, 2 exponents)
    const PG_Cart_MnD::NR_Evaluator& sev = ibs;

    Matrix1E_Adapter<PG_Cart_MnD::NR_Evaluator> mev(sev);    // the matrix-delivery view of it
    ASSERT_EQ(mev.size(), ibs.GetNumFunctions());

    const BasisSet::Orbital_1E_IBS<double>& bs1e = ibs;      // the scalar-built reference matrices
    const rsmat_t& S = bs1e.Overlap();
    const rsmat_t& K = bs1e.Kinetic();
    const rsmat_t& V = bs1e.Nuclear(&h2o);

    rsmat_t So=mev.OverlapMatrix(), Ko=mev.KineticMatrix(), Vo=mev.NuclearMatrix(&h2o);
    for (size_t i=0;i<mev.size();++i)
        for (size_t j=i;j<mev.size();++j)
        {
            EXPECT_NEAR(So(i,j), S(i,j), 1e-12) << "overlap ("<<i<<","<<j<<")";
            EXPECT_NEAR(Ko(i,j), K(i,j), 1e-12) << "kinetic ("<<i<<","<<j<<")";
            EXPECT_NEAR(Vo(i,j), V(i,j), 1e-12) << "nuclear ("<<i<<","<<j<<")";
        }
}
