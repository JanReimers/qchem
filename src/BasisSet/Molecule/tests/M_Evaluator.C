// File: BasisSet/Molecule/tests/M_Evaluator.C
//
// Increment-2 guard for the molecular Evaluator (Goal B): the concrete PG_Evaluator's 1E kernels must
// reproduce, element for element, the integrals the PG IBS already produces via MakeIntegrals.  An
// Orbital_IBS IS-A PGData, so we wrap it directly and compare PG_Evaluator.{Overlap,Grad2,Nuclear}(i,j)
// against ibs.{Overlap,Kinetic,Nuclear}().  (Kinetic() is the <p^2>=<-nabla^2> block, no 1/2 -- so it
// must equal Grad2.)  Absolute correctness of the underlying primitives is covered by M_PG_Oracle; this
// just proves the evaluator is a faithful, normalized view of them.
#include "gtest/gtest.h"
#include <vector>

import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD;                 // NR_Evaluator
import qchem.BasisSet.Molecule.PG_Cart;                      // Orbital_IBS
import qchem.BasisSet.Orbital_1E_IBS;                                  // cached Overlap()/Kinetic()/Nuclear() accessors
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD.PGData;      // PGData (the base of Orbital_IBS)
import qchem.Structure;                                                  // Molecule, Atom
import qchem.Types;
import qchem.Blaze;
using namespace qchem;

using BasisSet::Molecule::Evaluators::PG_Cart_MnD::NR_Evaluator;
using BasisSet::Molecule::PG_Cart::Orbital_IBS;

TEST(M_Evaluator, kernels_match_IBS_integrals)
{
    Molecule h2o;
    h2o.Insert(new Atom(8, 0, rvec3_t(0, 0,      0.117)));   // O on the C2 (z) axis
    h2o.Insert(new Atom(1, 0, rvec3_t(0, 0.757, -0.467)));
    h2o.Insert(new Atom(1, 0, rvec3_t(0,-0.757, -0.467)));

    rvec_t exps{1.0, 0.25};
    Orbital_IBS ibs(exps, 1, &h2o);                          // s + p shells, 2 exponents each

    const NR_Evaluator& ev = ibs;                            // the IBS IS-A evaluator; structure is per-call
    ASSERT_EQ(ev.size(), ibs.GetNumFunctions());

    // The cached 1E matrix accessors live on the IBS interface.  Reach them through the interface
    // reference: on the concrete Orbital_IBS the bare names now collide with the evaluator's element
    // kernels (Overlap(i,j), Nuclear(i,j,cl)), since the IBS IS-A PG_Evaluator.
    const BasisSet::Orbital_1E_IBS<double>& bs1e = ibs;
    const rsmat_t& S = bs1e.Overlap();       // normalized 1E matrices the IBS already builds
    const rsmat_t& K = bs1e.Kinetic();       // the <p^2> block (no 1/2)
    const rsmat_t& V = bs1e.Nuclear(&h2o);

    for (size_t i=0; i<ev.size(); ++i)
        for (size_t j=i; j<ev.size(); ++j)
        {
            EXPECT_NEAR(ev.Overlap(i,j), S(i,j), 1e-12) << "overlap ("<<i<<","<<j<<")";
            EXPECT_NEAR(ev.Grad2  (i,j), K(i,j), 1e-12) << "grad2/kinetic ("<<i<<","<<j<<")";
            EXPECT_NEAR(ev.Nuclear(i,j,&h2o), V(i,j), 1e-12) << "nuclear ("<<i<<","<<j<<")";
        }
}
