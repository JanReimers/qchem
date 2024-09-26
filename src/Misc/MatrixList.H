// File: MatrixList.H  A list of matricies which are mostly stored on disk.
#ifndef _MatrixList_H_
#define _MatrixList_H_



#include "Misc/pmstream.h"
#include "oml/smatrix.h"
#include "oml/matrix.h"
#include <vector>
#include <string>

//---------------------------------------------------------------------------------
//
//  This is a read only list of matricies.  Use disk files to store all matricies.
//
template <class T> class MatrixList
    : public virtual PMStreamableObject
{
public:
    typedef SMatrix<T> SMat;

    MatrixList();

    virtual ~MatrixList();

    void    Add(const SMat&);
    index_t size() const
    {
        return itsFileNames.size();
    }

    const SMat& operator[](index_t i) const;

    void Clear() const; //Release all memory used by buffered matricies.
    void Empty()      ; //Do a Clear() and also wipe all disk data.

    virtual std::ostream&    Write(std::ostream&) const;
    virtual std::istream&    Read (std::istream&)      ;
    virtual MatrixList* Clone(        ) const;

private:
    using StreamableObject::SetToBinary;
    using StreamableObject::Pretty;

    void LoadMatrix (index_t i) const;
    MatrixList(const MatrixList&);
    MatrixList& operator=(const MatrixList&);

    std::vector<std::string> itsFileNames;
    mutable SMat             itsMatrix;
    mutable index_t          itsMatrixIndex;
};
#endif // _MatrixList_H_
