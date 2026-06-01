// Unit tests for libxc presence and basic initialization using Google Test
#include <src/xc.h>
#include <gtest/gtest.h>
#include <iostream>

using std::cout;
using std::endl;

TEST(Libxc, Version)
{
    int vmajor, vminor, vmicro;
    xc_version(&vmajor, &vminor, &vmicro);
    cout << "Libxc version: " << vmajor << "." << vminor << "." << vmicro << endl;
    EXPECT_GE(vmajor, 7) << "libxc not available (xc_version returned " << vmajor << ")";
}

TEST(Libxc, CommonFunctionals)
{
    int ids[] = { XC_FAMILY_LDA, XC_FAMILY_GGA, XC_FAMILY_LCA };
    double rho[5] = {0.1, 0.2, 0.3, 0.4, 0.5};
    double sigma[5] = {0.2, 0.3, 0.4, 0.5, 0.6};
    double exc[5];
    xc_func_type func;
    for (int id : ids) {
        EXPECT_EQ(xc_func_init(&func,id, XC_UNPOLARIZED),0);
        ASSERT_NE(func.info, nullptr) << "no info for functional id " << id;
        EXPECT_EQ(func.info->number, id) << "mismatch id for functional (requested " << id << ", got " << (func.info ? func.info->number : -1) << ")";
        switch(func.info->family)
        {
            case XC_FAMILY_LDA:
                xc_lda_exc(&func, 5, rho, exc);
                break;
            case XC_FAMILY_GGA:
                xc_gga_exc(&func, 5, rho, sigma, exc);
                break;
            case XC_FAMILY_LCA:
                xc_lda_exc(&func, 5, rho, exc);
                break;
            }
            /* Print out density and energy density per particle */
            cout << func.info->family << " " << func.info->name << endl;
            /* Print out references */
            for(size_t ii = 0; func.info->refs[ii] != NULL; ii++)
                cout << func.info->refs[ii]->ref;
            for(size_t i=0; i<5; i++){
                cout << rho[i] << " " << exc[i] << endl;
        } //swtich
        xc_func_end(&func);
    }
}

