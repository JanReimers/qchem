#include "ERIList.H"
#include "Misc/stl_io.h"

ERIList::ERIList(int n)
    : itsN()
    , itsData()
{
    SetSize(n);
};

ERIList::ERIList()
    : itsN(0)
    , itsData()
{};

ERIList::~ERIList()
{
}

//Clear out all data.
void ERIList::Empty()
{
    itsData.clear();
    itsN=0;
}

void ERIList::SetSize(int n)
{
    itsN=n;
    itsData.resize(GetIndex(n,n,n,n,n)+1,0);
}

std::ostream&    ERIList::Write(std::ostream& os) const
{
    return os << itsN << " " << itsData << std::endl;
}

std::istream&    ERIList::Read (std::istream& is)
{
    return is >> itsN >> itsData;
}

#include <iomanip>
void ERIList::Dump(std::ostream& os) const
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

void ERIList::DumpExchange(std::ostream& os) const
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
