// FIle: SCFAccelerator.C  Interface for an accelerator alogrithm


#include "Imp/SCFAccelerator.H"
#include <LASolver.H>

SCFIrrepAccelerator_DIIS::SCFIrrepAccelerator_DIIS(const SMat& _S) : S(_S) {};
SCFIrrepAccelerator_DIIS::~SCFIrrepAccelerator_DIIS() {};
void SCFIrrepAccelerator_DIIS::Init(const LASolver<double>* las)
{
    itsLaSolver=las;
}

SCFIrrepAccelerator::SMat SCFIrrepAccelerator_DIIS::Project(const SMat& F, const SMat& DPrime)
{
    Mat Fproj_nonsym=itsLaSolver->Transform(F);
    double del=MakeSymmetric(Fproj_nonsym);
    assert(del<1e-15);
    SMat Fproj(Fproj_nonsym); //Fproj = Fprime = Vd*F*V
    if (Max(fabs(DPrime))==0.0) return Fproj;
    Mat E=Fproj*DPrime-DPrime*Fproj;// Make the communtator
    double En=FrobeniusNorm(E);
    if (En<1.0)
    {
        // std::cout << "||E|| = " << FrobeniusNorm(E) << std:: endl;
        itsEs.push_back(E);
        itsF_Primes.push_back(Fproj);
        Purge(8);
        if (itsEs.size()>1  && En>1e-6) 
            Fproj=Solve();
        else
            if (itsEs.size()>1 ) std::cout << "Skip projection ||E|| = " << En << std::endl;
    }
    return Fproj;
}

template <class T> const SMatrix<T>& operator+=(SMatrix<T>& a, const SMatrix<T>& b)
{
    if (a.size()==0) 
    {
        a.SetLimits(b.GetLimits());
        Fill(a,0.0);
    }
    return ArrayAdd(a,b);
}

#include "oml/numeric/LapackLinearSolver.H"
#include "oml/vector.h"
SCFIrrepAccelerator::SMat SCFIrrepAccelerator_DIIS::Solve()
{
    index_t N=itsEs.size()+1;
    Mat B(N,N);
    index_t i=1;
    for (auto ei:itsEs)
    {
        index_t j=1;
        for (auto ej:itsEs)
            B(i,j++)=Dot(ei,ej);
        i++;
    }
    for (index_t i=1;i<=N-1;i++)
    {
        B(N,i)=B(i,N)=-1.0;
    }
    B(N,N)=0.0;
    Vector<double> v(N);
    Fill(v,0.0); 
    v(N)=-1.0;
    auto solver=new LapackLinearSolver<double>();
    Vector<double> C=solver->Solve(B,v);
    delete solver;
    Vector<double> del=B*C-v;
    // StreamableObject::SetToPretty();
    // std:: cout << B << C << v << del <<std::endl;
    // std:: cout << "del,[F,D] = " << sqrt(del*del) << " " << C(N) << std::endl;
    assert(sqrt(del*del)<1e-14);
    Vector<double> c=C.SubVector(N-1);
    // std::cout << " sum(c)-1 = " << Sum(c)-1.0 << std::endl;
    assert(fabs(Sum(c)-1.0)<1e-13);
    SMat Fproj;
    i=1;
    for (auto f:itsF_Primes) 
        Fproj+=SMat(c(i++)*f);
    return Fproj;

}
void SCFIrrepAccelerator_DIIS::Purge(int N)
{
    assert(itsEs.size()==itsF_Primes.size());
    if (itsEs.size()>N)
    {
        itsEs.pop_front();
        itsF_Primes.pop_front();
    }
}


#include <Irrep_BS.H>

SCFAccelerator_DIIS::SCFAccelerator_DIIS() {};
SCFAccelerator_DIIS::~SCFAccelerator_DIIS() {};
SCFIrrepAccelerator* SCFAccelerator_DIIS::Create(const TOrbital_IBS<double>* bs) const
{
    return new SCFIrrepAccelerator_DIIS(bs->Overlap());
}

