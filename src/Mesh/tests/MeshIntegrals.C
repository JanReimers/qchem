// File: src/Mesh/tests/MeshIntegrals.C
#include "gtest/gtest.h"
#include <nlohmann/json.hpp>
// #include <iostream>
import qchem.Mesh.Factory;
import qchem.VectorFunction;
import qchem.Blaze;

class MeshIntegrals : public ::testing::Test
{
public:
};

class Exp : public virtual VectorFunction<double>
{
public:
    virtual size_t  GetVectorSize() const {return 1;};

    virtual rvec_t operator()  (const rvec3_t& r) const
    {
        return {exp(-2*norm(r))};
    }
    virtual rvec3vec_t Gradient(const rvec3_t& r) const
    {
        double mr=norm(r);
        if (mr==0.0) return {rvec3_t(0,0,0)};
        rvec3_t rhat=r/mr;
        return {-2*rhat*exp(-2*mr)};
    }

};

TEST_F(MeshIntegrals,MHL)
{
    nlohmann::json js={{"N",200},{"m",2},{"alpha",3.0}};
    RadialMesh* m=MeshF::Factory(qchem::RadialType::MHL,js);
    Exp f;

    size_t n=f.GetVectorSize();
    rvec_t I(n,0.0);

    int iw=0;
    for (auto rw:*m)
    {
        rvec_t fr=f(rvec3_t(r(rw),0,0));
        for (size_t i=0; i<n; i++)
            I[i]+=fr[i]*w(rw);
        iw++;
    }

    EXPECT_NEAR(I[0],0.25,4e-15);
}