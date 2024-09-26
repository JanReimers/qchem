// File: MatrixList.C  A list of matricies which are mostly stored on disk.



#include "Misc/MatrixList.H"
#include "Misc/stl_io.h"
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>

template <class T> MatrixList<T>::MatrixList()
    : itsMatrixIndex(-1)
{};

template <class T> MatrixList<T>::MatrixList(const MatrixList<T>& m)
    : itsFileNames(m.itsFileNames)
    , itsMatrixIndex(-1)
{};

template <class T> MatrixList<T>::~MatrixList()
{
    Empty();
}

template <class T> void MatrixList<T>::Add(const SMat& m)
{
    char fname[L_tmpnam+100];
    char* ret=tmpnam(fname);
    (void)*ret; //Avoid unused warning.
    std::string t(fname);
    itsFileNames.push_back(t);
    std::ofstream mfile(itsFileNames.back().c_str());
    StreamableObject::Mode  mode=SetToBinary();
    mfile << m;
    StreamableObject::SetOutputMode(mode);
}

template <class T> const typename MatrixList<T>::SMat& MatrixList<T>::operator[](index_t i) const
{
    LoadMatrix(i);
    return itsMatrix;
}

template <class T> void MatrixList<T>::Clear() const
{
    itsMatrix.SetLimits(MatLimits(0,0));
    itsMatrixIndex=-1;
}

template <class T> void MatrixList<T>::Empty()
{
    Clear();

    for (std::vector<std::string>::iterator i(itsFileNames.begin()); i!=itsFileNames.end(); i++) unlink(i->c_str());
    itsFileNames.clear();
}

template <class T> std::ostream& MatrixList<T>::Write(std::ostream& os) const
{
    os << size() << " ";
    if (Pretty())
        os << "MatrixList with files " << itsFileNames << std::endl;
    else
        for (index_t i=0; i<size(); i++) os << (*this)[i];

    return os;
}

template <class T> std::istream& MatrixList<T>::Read (std::istream& is)
{
    Empty();
    index_t N;
    is >> N;
    is.get(); //mop up the space.
    for (index_t i=0; i<N; i++)
    {
        SMat m;
        is >> m;
        Add(m);
    }
    return is;
}

template <class T> MatrixList<T>* MatrixList<T>::Clone() const
{
    return new MatrixList(*this);
}


template <class T> void MatrixList<T>::LoadMatrix(index_t i) const
{
//    std::cout << "LoadMatrix i=" << i << " itsMatrixIndex=" << itsMatrixIndex << std::endl;
//    std::cout << itsFileNames.size() << std::endl;
    if (i!=itsMatrixIndex)
    {
//        std::cout << "Filename=" << itsFileNames[i] << std::endl;
        std::ifstream mfile(itsFileNames[i].c_str());
        if(!mfile)
        {
            std::cerr << "MatrixList::op(i) could not open matrix file '" << itsFileNames[i] << "'" << std::endl;
            exit(-1);
        }
        itsMatrix.SetLimits(MatLimits(0,0));
        mfile >> itsMatrix;
        itsMatrixIndex=i;
    }
}

template class MatrixList<double>;
template class MatrixList<std::complex<double> >;
