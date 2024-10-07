#include "ERI4.H"

ERI4::ERI4(int n)
    : itsN()
    , itsData()
{
    SetSize(n);
};
ERI4::ERI4(int n,double _fill)
    : itsN()
    , itsData()
{
    SetSize(n,_fill);
};

ERI4::ERI4()
    : itsN(0)
    , itsData()
{};

ERI4::~ERI4()
{
    
}

//Clear out all data.
void ERI4::Empty()
{
    itsData.SetLimits(0);
    itsN=0;
}

void ERI4::SetSize(int n,double _fill)
{
    itsN=n;
    itsData.SetLimits(0,GetIndex(n,n,n,n,n),false);
    Fill(itsData,_fill);
//    for (index_t i=0;i<itsData.size();i++) itsData[i]=_fill;
}

std::ostream&    ERI4::Write(std::ostream& os) const
{
    return os << itsN << " " << itsData << std::endl;
}

std::istream&    ERI4::Read (std::istream& is)
{
    return is >> itsN >> itsData;
}

#include <iomanip>
void ERI4::Dump(std::ostream& os) const
{
    os.precision(4);
    for (int ia=1; ia<=itsN; ia++)
        for (int ib=ia; ib<=itsN; ib++)
        {
            SMatrix<double> Jcd(itsN,itsN);
            for (int ic=1; ic<=itsN; ic++)
                for (int id=ic; id<=itsN; id++)
                {
                    Jcd(ic,id)=(*this)(ia,ib,ic,id);
                }
            os << "J(" << ia << "," << ib<<")="<< std::endl<< std::setw(7) << Jcd << std::endl;
        }
}

void ERI4::DumpExchange(std::ostream& os) const
{
    os.precision(4);
    for (int ia=1; ia<=itsN; ia++)
        for (int ib=ia; ib<=itsN; ib++)
        {
            SMatrix<double> Kcd(itsN,itsN);
            for (int ic=1; ic<=itsN; ic++)
                for (int id=ic; id<=itsN; id++)
                {
                    Kcd(ic,id)=(*this)(ia,id,ic,ib);
                }
            os << "K(" << ia << "," << ib<<")="<< std::endl<< std::setw(7) << Kcd << std::endl;
        }
}


#include "ERI4.H"

ERI4view::ERI4view(ERI4& eril,int start_ab, int start_cd)
 : itsERI4(eril)
 , itsStart_a(start_ab)
 , itsStart_b(start_ab)
 , itsStart_c(start_cd)
 , itsStart_d(start_cd)
{
    assert(itsStart_a>0);
    assert(itsStart_b>0);
    assert(itsStart_c>0);
    assert(itsStart_d>0);
    //ctor
}

ERI4view::ERI4view(ERI4& eril,int sa,int sb,int sc, int sd)
 : itsERI4(eril)
 , itsStart_a(sa)
 , itsStart_b(sb)
 , itsStart_c(sc)
 , itsStart_d(sd)
{
    assert(itsStart_a>0);
    assert(itsStart_b>0);
    assert(itsStart_c>0);
    assert(itsStart_d>0);
    //ctor
}

