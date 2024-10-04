#include "ERIList1.H"

ERIList1::ERIList1(int n)
    : itsN()
    , itsData(0)
{
    SetSize(n);
};
ERIList1::ERIList1(int n,double _fill)
    : itsN()
    , itsData(0)
{
    SetSize(n,_fill);
};

ERIList1::ERIList1()
    : itsN(0)
    , itsData(0)
{};

ERIList1::~ERIList1()
{
    
}

//Clear out all data.
void ERIList1::Empty()
{
    itsData=cow_array<double>(0);
    itsN=0;
}

void ERIList1::SetSize(int n,double _fill)
{
    itsN=n;
    itsData=cow_array<double>(GetIndex(n,n,n,n,n)+1);
    for (index_t i=0;i<itsData.size();i++) itsData[i]=_fill;

    // Fill?
}

std::ostream&    ERIList1::Write(std::ostream& os) const
{
//    return os << itsN << " " << itsData << std::endl;
    return os;
}

std::istream&    ERIList1::Read (std::istream& is)
{
//    return is >> itsN >> itsData;
    return is;
}

#include <iomanip>
void ERIList1::Dump(std::ostream& os) const
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

void ERIList1::DumpExchange(std::ostream& os) const
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
