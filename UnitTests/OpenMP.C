// File: OpenMP.C  

#include "gtest/gtest.h"
#include <omp.h>
#include <chrono>
import qchem.BasisSet.ERI4;

class OpenMPTests : public ::testing::Test
{
public:
    OpenMPTests()
    {
        StreamableObject::SetToPretty();
    }
    auto now() {return std::chrono::high_resolution_clock::now();}
    
};


ERI4::SMat ContractSerial(const ERI4::SMat& Sab, const ERI4& gabcd)
{
    ERI4::SMat Scd(gabcd(1,1).GetLimits());
    Fill(Scd,0.0);
    for (auto ia:Sab.rows())
        Scd+=gabcd(ia,ia)*Sab(ia,ia);

    for (auto ia:Sab.rows())
        for (auto ib:Sab.cols(ia+1))
            Scd+=2*gabcd(ia,ib)*Sab(ia,ib);

    return Scd;
}
ERI4::SMat ContractParallel(const ERI4::SMat& Sab, const ERI4& gabcd)
{
    ERI4::SMat Scd(gabcd(1,1).GetLimits());
    Fill(Scd,0.0);
    #pragma omp parallel for 
    for (size_t ia=1;ia<Sab.GetNumRows();ia++)
    {
        ERI4::SMat s=gabcd(ia,ia)*Sab(ia,ia);
        // #pragma omp critical
        {
            Scd+=s;
        }
    }
        

    #pragma omp parallel for collapse(2) 
    for (size_t ia=1;ia<Sab.GetNumRows();ia++)
        for (size_t ib=ia+1;ib<Sab.GetNumCols();ib++)
        {
            ERI4::SMat s=2*gabcd(ia,ib)*Sab(ia,ib);
            // #pragma omp critical
            {
                Scd+=s;
            }

        }
    return Scd;
}


TEST_F(OpenMPTests,Dot)
{
    
    std::cout << "n   threads=" << omp_get_num_threads() << std::endl;
    std::cout << "max threads=" << omp_get_max_threads() << std::endl;
    #pragma omp parallel
    {
        #pragma omp single
        std::cout << "n   threads=" << omp_get_num_threads() << " #=" << omp_get_thread_num() <<std::endl;
    }
    size_t N=200;
    ERI4 g(N,N);
    for (size_t ia=1;ia<N;ia++)
        for (size_t ib=ia;ib<N;ib++)
            FillRandom(g(ia,ib));


    double tp,ts;
    SMatrix<double> A(N),Bs(N),Bp(N);
    FillRandom(A);


    {
        auto start_time = now();
        Bp=ContractParallel(A,g);
        auto end_time = now();
        std::chrono::duration<double> serial_duration  = end_time - start_time;
        tp=serial_duration.count()*1e3;
        std::cout << "Parallel Elapsed time = " << tp << "(ms)" << std::endl;
    }

    {
        auto start_time = now();
        Bs=ContractSerial(A,g);
        auto end_time = now();
        std::chrono::duration<double> serial_duration  = end_time - start_time;
        ts=serial_duration.count()*1e3;
        std::cout << "Serial   Elapsed time = " << ts << "(ms)" << std::endl;
    }
     
    std::cout.precision(3);
    // std::cout << Bp-Bs << std::endl;
    std::cout << Max(fabs(Bs-Bp))/N << " ts/tp=" << ts/tp << std::endl;
    EXPECT_NEAR(Max(fabs(Bs-Bp))/N,0.0,2e-12);
    // EXPECT_EQ(Sum(fabs(Bs-Bp)),0.0);
    EXPECT_GE(ts/tp,3.0);
    {
        auto start_time = now();
        Bp=ContractParallel(A,g);
        auto end_time = now();
        std::chrono::duration<double> serial_duration  = end_time - start_time;
        tp=serial_duration.count()*1e3;
        std::cout << "Parallel Elapsed time = " << tp << "(ms)" << std::endl;
    }
    EXPECT_GE(ts/tp,3.0);

}