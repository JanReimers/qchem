module;
// #include <cassert>
module qchem.BasisSet.Internal.ERI4T;

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


template class ERI4T<double,SMatrix>; //g++ 15.2 BUG internal compiler error: in import_export_decl, at cp/decl2.cc:3563
template class ERI4T<double, Matrix>; //g++ 15.2 BUG internal compiler error: in import_export_decl, at cp/decl2.cc:3563
