#include "Imp/Containers/ERI4.H"

template <class T,template<class> class M> ERI4T<T,M>::ERI4T(size_t Nab, size_t Ncd) : itsData(Nab,Nab)
{
    M<T> Jcd(Ncd,Ncd);
    Fill(Jcd,0.0);
    Fill(itsData,Jcd);
}

template <class T,template<class> class M> size_t ERI4T<T,M>::size() const
{
    size_t Nab=itsData.size();
    size_t Ncd=itsData(1,1).size();
    return Nab*Ncd;
}

template class ERI4T<double,SMatrix>;
template class ERI4T<double, Matrix>;

ERI4::SMat operator*(const ERI4& gabcd,const ERI4::SMat& Scd)
{
    //std::cout << "gabcd=" << gabcd.GetLimits() << " Scd=" << Scd.GetLimits() << std::endl;
    ERI4::SMat Sab(gabcd.GetLimits());
    for (auto ia:Sab.rows())
        for (auto ib:Sab.cols(ia))
            Sab(ia,ib)=ERI4::contract(gabcd(ia,ib),Scd); //Dot(DirectMultiply(A,B))
    return Sab;
}

// Profiling hot loop
ERI4::SMat operator*(const ERI4::SMat& Sab, const ERI4& gabcd)
{
    ERI4::SMat Scd(gabcd(1,1).GetLimits());
    Fill(Scd,0.0);
    for (auto ia:Sab.rows())
    {
        Scd+=gabcd(ia,ia)*Sab(ia,ia);
        for (auto ib:Sab.cols(ia+1))
            Scd+=2*gabcd(ia,ib)*Sab(ia,ib);
    }
    return Scd;
}

double ERI4::contract(const ERI4::SMat& A,const ERI4::SMat& B)
{
    //std::cout << "ERI4::contract " << A.GetLimits() << " " << B.GetLimits() << std::endl;
    assert(A.GetLimits()==B.GetLimits());
    double ret=Dot(A.GetDiagonal(),B.GetDiagonal());
    for (auto ia:A.rows())
        for (auto ib:A.cols(ia+1))
            ret+=2*A(ia,ib)*B(ia,ib);
    return ret;         
}

