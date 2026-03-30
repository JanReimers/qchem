module;
module qchem.BasisSet.Internal.ERI4T;
import qchem.Blaze;

template <> ERI4T<double,mat_t>::ERI4T(size_t Nab, size_t Ncd) : itsData(Nab,Nab)
{
    mat_t<double> Jcd(Ncd,Ncd,0.0);
    itsData=Jcd;
}
template <> ERI4T<double,smat_t>::ERI4T(size_t Nab, size_t Ncd) : itsData(Nab)
{
     smat_t<double> Jcd=zero<double>(Ncd);
     for (auto i:iv_t(0,Nab))
        for (auto j:iv_t(i,Nab)) itsData(i,j)=Jcd;
}

template <class T,template<class> class M> size_t ERI4T<T,M>::size() const
{
    size_t Nab=itsData.rows();
    size_t Ncd=itsData(0,0).rows();
    return Nab*Nab*Ncd*Ncd;
}


template class ERI4T<double,smat_t>; //g++ 15.2 BUG internal compiler error: in import_export_decl, at cp/decl2.cc:3563
template class ERI4T<double, mat_t>; //g++ 15.2 BUG internal compiler error: in import_export_decl, at cp/decl2.cc:3563
