// File: BasisSet/Molecule/tests/M_MEvaluator.C
//
// Guard for the matrix-delivery evaluator category (isM_1E / isM_DFT / isM_HF) and the IBS mixins'
// dispatch.  A Matrix_Adapter wraps a scalar is1E_DFT_HF_Evaluator and pre-assembles its 1E matrices,
// 3-centre ERI3s and 4-centre ERI4s -- a stand-in for a block-oriented integral library / disk source
// (the eventual libcint wrapper).  We check it satisfies the isM_ concepts (and NOT the scalar ones), that
// every matrix/ERI it delivers equals what the scalar IBS builds via the public accessors, and that all
// three mixins' matrix branches compile (explicit instantiations).
#include "gtest/gtest.h"
#include <ostream>
#include <string>
#include <vector>

import qchem.BasisSet.Molecule.Evaluators;                    // Evaluator + the is*_/isM_* concepts
import qchem.BasisSet.Molecule.Evaluators.PG_Cart_MnD;        // NR_Evaluator (a scalar 1E+DFT+HF evaluator)
import qchem.BasisSet.Molecule.PG_Cart;                       // Orbital_IBS
import qchem.BasisSet.Molecule.IBS;                           // Orbital_{1E,DFT,HF}_IBS<E> mixins (dispatch)
import qchem.BasisSet.Orbital_1E_IBS;                         // public Overlap()/Kinetic()/Nuclear()
import qchem.BasisSet.Orbital_DFT_IBS;                        // public Overlap3C()/Repulsion3C() + Fit_IBS
import qchem.BasisSet.Orbital_HF_IBS;                         // public Direct()/Exchange()
import qchem.BasisSet.Fit_IBS;
import qchem.BasisSet.Internal.ERI3;                          // ERI3<double>, fnorm
import qchem.BasisSet.Internal.ERI4;                          // ERI4, fnorm
import qchem.Structure;
import qchem.Types;
import qchem.Blaze;

using namespace BasisSet::Molecule::Evaluators;
using BasisSet::Molecule::PG_Cart::Orbital_IBS;
using NRE = PG_Cart_MnD::NR_Evaluator;

// A full matrix-delivery evaluator: wraps a scalar is1E_DFT_HF_Evaluator and assembles its matrices / ERIs
// on demand.  (A real isM_ source -- libcint, a disk cache -- would produce these directly; here we build on
// the scalar evaluator so the test has a known-correct reference, the scalar IBS, to compare against.)
template <is1E_DFT_HF_Evaluator E>
class Matrix_Adapter : public virtual Evaluator
{
    const E& its;
    template <class K> rsmat_t Build1(K k) const
    {
        rsmat_t S(its.size());
        for (auto i:its.indices()) for (auto j:its.indices(i)) S(i,j)=k(i,j);
        return S;
    }
    template <class K> ERI3<double> Build3(const Matrix_Adapter& fit, K k) const
    {
        size_t Na=its.size(), Nc=fit.its.size();
        ERI3<double> s3;
        for (size_t ic=0; ic<Nc; ic++)
        {
            rsmat_t s(Na);
            for (size_t ia=0; ia<Na; ia++) for (size_t ib=ia; ib<Na; ib++) s(ia,ib)=k(its, ia, ib, fit.its, ic);
            s3.push_back(s);
        }
        return s3;
    }
public:
    Matrix_Adapter(const E& e) : its(e) {}
    virtual size_t        size() const {return its.size();}
    virtual rvec_t        Norm() const {return its.Norm();}
    virtual std::string   Name() const {return "Matrix_Adapter";}
    virtual std::ostream& Write(std::ostream& os) const {return os << "Matrix_Adapter[" << size() << "]";}

    // --- isM_1E ---
    rsmat_t OverlapMatrix()                  const {return Build1([&](size_t i,size_t j){return its.Overlap(i,j);});}
    rsmat_t KineticMatrix()                  const {return Build1([&](size_t i,size_t j){return its.Grad2  (i,j);});}
    rsmat_t NuclearMatrix(const Structure* cl) const {return Build1([&](size_t i,size_t j){return its.Nuclear(i,j,cl);});}

    // --- isM_DFT (3-centre) ---
    ERI3<double> OverlapThreeC_Matrix(const Matrix_Adapter& fit) const
    { return Build3(fit, [](const E& a,size_t ia,size_t ib,const E& c,size_t ic){return a.OverlapThreeC  (ia,a,ib,c,ic);}); }
    ERI3<double> RepulsionThreeC_Matrix(const Matrix_Adapter& fit) const
    { return Build3(fit, [](const E& a,size_t ia,size_t ib,const E& c,size_t ic){return a.RepulsionThreeC(ia,a,ib,c,ic);}); }

    // --- isM_HF (4-centre).  The Direct loop is trivial, but ExchangeMatrix has to reproduce
    // Orbital_HF_IBS::MakeExchange's symmetry packing VERBATIM -- that an opaque matrix wrapper must
    // duplicate this is exactly the friction the (later) block-granularity category removes. ---
    ERI4 DirectMatrix(const Matrix_Adapter& p) const
    {
        const E& a=its; const E& c=p.its; size_t Na=a.size(), Nc=c.size();
        ERI4 J(Na,Nc);
        for (size_t ia:iv_t(0,Na)) for (size_t ib:iv_t(ia,Na))
        {
            rsmat_t& Jab=J(ia,ib);
            for (size_t ic:iv_t(0,Nc)) for (size_t id:iv_t(ic,Nc)) Jab(ic,id)=a.FourC(ia,a,ib,c,ic,c,id);
        }
        return J;
    }
    ERI4 ExchangeMatrix(const Matrix_Adapter& pb) const
    {
        const E& a=its; const E& b=pb.its; size_t Na=a.size(), Nb=b.size();
        ERI4 K(Na,Nb);
        for (size_t ia:iv_t(0,Na)) for (size_t ib:iv_t(0,Nb)) for (size_t ic:iv_t(ia,Na))
        {
            rsmat_t& Kac=K(ia,ic);
            for (size_t id:iv_t(0,Nb))
            {
                double v=a.FourC(ia,b,ib,a,ic,b,id);
                if (ib==id) Kac(ib,id)=v; else if (ib<id) Kac(ib,id)+=0.5*v; else Kac(id,ib)+=0.5*v;
            }
        }
        return K;
    }
};

static_assert( isM_1E_DFT_HF_Evaluator<Matrix_Adapter<NRE>>, "Matrix_Adapter is a full matrix evaluator");
static_assert(!is1E_Evaluator        <Matrix_Adapter<NRE>>, "a matrix evaluator does NOT provide scalar kernels");

// Compile-check the matrix branch of all three mixins (PG_Cart already covers the scalar branches).
template class BasisSet::Molecule::Orbital_1E_IBS <Matrix_Adapter<NRE>>;
template class BasisSet::Molecule::Orbital_DFT_IBS<Matrix_Adapter<NRE>>;
template class BasisSet::Molecule::Orbital_HF_IBS <Matrix_Adapter<NRE>>;

static Molecule* MakeWater()
{
    Molecule* w = new Molecule();
    w->Insert(new Atom(8, 0, rvec3_t(0, 0,      0.117)));
    w->Insert(new Atom(1, 0, rvec3_t(0, 0.757, -0.467)));
    w->Insert(new Atom(1, 0, rvec3_t(0,-0.757, -0.467)));
    return w;
}

TEST(M_MEvaluator, matrix_1E_matches_scalar)
{
    Molecule* h2o = MakeWater();
    rvec_t exps{1.0, 0.25};
    Orbital_IBS ibs(exps, 1, h2o);                           // scalar PG evaluator (s + p, 2 exponents)
    Matrix_Adapter<NRE> mev(static_cast<const NRE&>(ibs));   // its matrix-delivery view
    ASSERT_EQ(mev.size(), ibs.GetNumFunctions());

    const BasisSet::Orbital_1E_IBS<double>& bs1e = ibs;
    rsmat_t So=mev.OverlapMatrix(), Ko=mev.KineticMatrix(), Vo=mev.NuclearMatrix(h2o);
    for (size_t i=0;i<mev.size();++i)
        for (size_t j=i;j<mev.size();++j)
        {
            EXPECT_NEAR(So(i,j), bs1e.Overlap()(i,j),    1e-12) << "overlap ("<<i<<","<<j<<")";
            EXPECT_NEAR(Ko(i,j), bs1e.Kinetic()(i,j),    1e-12) << "kinetic ("<<i<<","<<j<<")";
            EXPECT_NEAR(Vo(i,j), bs1e.Nuclear(h2o)(i,j), 1e-12) << "nuclear ("<<i<<","<<j<<")";
        }
}

TEST(M_MEvaluator, matrix_3C_4C_match_scalar)
{
    Molecule* h2o = MakeWater();
    rvec_t exps{1.0, 0.25};
    Orbital_IBS ibs(exps, 1, h2o);
    Matrix_Adapter<NRE> mev(static_cast<const NRE&>(ibs));

    // 4-centre (HF): partner = self.  Reference = the public, energy-validated scalar Direct/Exchange.
    const BasisSet::Orbital_HF_IBS<double>& hf = ibs;
    EXPECT_LT(fnorm(mev.DirectMatrix(mev),   hf.Direct(ibs)),   1e-10);
    EXPECT_LT(fnorm(mev.ExchangeMatrix(mev), hf.Exchange(ibs)), 1e-10);

    // 3-centre (DFT): a real Coulomb-fit basis from the orbital IBS; reference = public Overlap3C/Repulsion3C.
    BasisSet::Fit_IBS* fit = ibs.CreateCDFitBasisSet(h2o);
    Matrix_Adapter<NRE> mfit(dynamic_cast<const NRE&>(*fit));
    const BasisSet::Orbital_DFT_IBS<double>& dft = ibs;
    EXPECT_LT(fnorm(mev.OverlapThreeC_Matrix(mfit),   dft.Overlap3C(*fit)),   1e-10);
    EXPECT_LT(fnorm(mev.RepulsionThreeC_Matrix(mfit), dft.Repulsion3C(*fit)), 1e-10);
    delete fit;
}
