module;
#include <cassert>
#include <iostream>
#include <iomanip>

export module qchem.Matrix2D;
export import qchem.Vector2D;
import qchem.Math;

//------------------------------------------------------------------------
//
//  Hard coded 2x2 matrix, mirroring Matrix3D.  Default construction is the identity.
//
export template <class T> class Matrix2D
{
public:
   T M11,M12;
   T M21,M22;

   Matrix2D() : M11(1),M12(0),M21(0),M22(1) {}
   Matrix2D(const T& m11,const T& m12,const T& m21,const T& m22) : M11(m11),M12(m12),M21(m21),M22(m22) {}
   Matrix2D(const Matrix2D<T>& m) : M11(m.M11),M12(m.M12),M21(m.M21),M22(m.M22) {}

   Matrix2D& operator =(const Matrix2D<T>& m) {M11=m.M11;M12=m.M12;M21=m.M21;M22=m.M22;return *this;}

   Vector2D<T> GetRow(size_t r) const {return Vector2D<T>((*this)(r,1),(*this)(r,2));}
   Vector2D<T> GetCol(size_t c) const {return Vector2D<T>((*this)(1,c),(*this)(2,c));}

   const T& operator()(size_t i, size_t j) const {return (&M11)[2*(i-1)+j-1];}
         T& operator()(size_t i, size_t j)       {return (&M11)[2*(i-1)+j-1];}
};

export {

template <class T> inline bool operator==(const Matrix2D<T>& a,const Matrix2D<T>& b)
{ return a.M11==b.M11 && a.M12==b.M12 && a.M21==b.M21 && a.M22==b.M22; }
template <class T> inline bool operator!=(const Matrix2D<T>& a,const Matrix2D<T>& b) {return !(a==b);}

template <class T> inline Matrix2D<T> Transpose(const Matrix2D<T>& m) {return Matrix2D<T>(m.M11,m.M21,m.M12,m.M22);}
template <class T> inline Matrix2D<T> operator~(const Matrix2D<T>& m) {return Transpose(m);}

//
//  Algebra
//
template <class T1, class T2> inline
auto operator+(const Matrix2D<T1>& a,const Matrix2D<T2>& b) {return Matrix2D(a.M11+b.M11,a.M12+b.M12,a.M21+b.M21,a.M22+b.M22);}
template <class T1, class T2> inline
auto operator-(const Matrix2D<T1>& a,const Matrix2D<T2>& b) {return Matrix2D(a.M11-b.M11,a.M12-b.M12,a.M21-b.M21,a.M22-b.M22);}
template <class T1, class T2> inline
Matrix2D<T1>& operator+=(Matrix2D<T1>& a,const Matrix2D<T2>& b) {a.M11+=b.M11;a.M12+=b.M12;a.M21+=b.M21;a.M22+=b.M22;return a;}
template <class T1, class T2> inline
Matrix2D<T1>& operator-=(Matrix2D<T1>& a,const Matrix2D<T2>& b) {a.M11-=b.M11;a.M12-=b.M12;a.M21-=b.M21;a.M22-=b.M22;return a;}

template <class T> inline Matrix2D<T> operator*(const Matrix2D<T>& m,const T& c) {return Matrix2D<T>(m.M11*c,m.M12*c,m.M21*c,m.M22*c);}
template <class T> inline Matrix2D<T> operator*(const T& c,const Matrix2D<T>& m) {return m*c;}
template <class T> inline Matrix2D<T> operator/(const Matrix2D<T>& m,const T& c) {return Matrix2D<T>(m.M11/c,m.M12/c,m.M21/c,m.M22/c);}
template <class T> inline Matrix2D<T> operator-(const Matrix2D<T>& m) {return Matrix2D<T>(-m.M11,-m.M12,-m.M21,-m.M22);}

//
//  Matrix algebra
//
template <class T1, class T2> inline
auto operator*(const Matrix2D<T1>& a, const Matrix2D<T2>& b)
{
   return Matrix2D
   (
      a.M11*b.M11 + a.M12*b.M21,  a.M11*b.M12 + a.M12*b.M22,
      a.M21*b.M11 + a.M22*b.M21,  a.M21*b.M12 + a.M22*b.M22
   );
}

template <class T1, class T2> inline
auto operator*(const Matrix2D<T2>& m, const Vector2D<T1>& v) {return Vector2D(m.M11*v.x + m.M12*v.y, m.M21*v.x + m.M22*v.y);}
template <class T1, class T2> inline
auto operator*(const Vector2D<T1>& v, const Matrix2D<T2>& m) {return Vector2D(m.M11*v.x + m.M21*v.y, m.M12*v.x + m.M22*v.y);}

//! Outer product  a (x) b :  (a b^T)_{ij} = a_i b_j.
template <class T> inline Matrix2D<T> Outer(const Vector2D<T>& a, const Vector2D<T>& b)
{ return Matrix2D<T>(a.x*b.x, a.x*b.y, a.y*b.x, a.y*b.y); }

template <class T> inline T Determinant(const Matrix2D<T>& m) {return m.M11*m.M22 - m.M12*m.M21;}

template <class T> inline Matrix2D<T> Invert(const Matrix2D<T>& m)
{
   T d=(T)(1.0/Determinant(m));
   return Matrix2D<T>(m.M22*d, -m.M12*d, -m.M21*d, m.M11*d);
}

template <class T> std::ostream& operator<<(std::ostream& os,const Matrix2D<T>& a)
{
   std::streamsize prec=os.precision(), wid=os.width();
   os << std::setw(0);
   os << "[ " << std::setw(wid) << std::setprecision(prec) << a.M11 << " " << std::setw(wid) << std::setprecision(prec) << a.M12 << " ]" << std::endl;
   os << "[ " << std::setw(wid) << std::setprecision(prec) << a.M21 << " " << std::setw(wid) << std::setprecision(prec) << a.M22 << " ]" << std::endl;
   return os;
}

} //export block
