module;
#include <cassert>
#include <iostream>
#include <iomanip>
#include <complex>

export module qchem.Vector2D;
export import qchem.Math;

//------------------------------------------------------------------------
//
//  Light weight 2D vector, mirroring Vector3D.  Most algebraic operators are
//  overloaded so a Vector2D behaves like an intrinsic type.  \c * is the dot product.
//
export template <class T> class Vector2D
{
public:
  Vector2D(                 ): x(0),y(0) {}
  Vector2D(T _x,T _y        ): x(_x),y(_y) {}
  Vector2D(const Vector2D& v): x(v.x),y(v.y) {}
  template <class T1> Vector2D(const Vector2D<T1>& v) : x(v.x),y(v.y) {}

  Vector2D& operator =(const Vector2D& v) {x=v.x;y=v.y;return *this;}
  template <class T1> Vector2D& operator =(const Vector2D<T1>& v) {x=v.x;y=v.y;return *this;}

  T x; //!< \a x coordinate.
  T y; //!< \a y coordinate.
};

export {

template <class T1, class T2> inline
auto operator +(const Vector2D<T1>& a,const Vector2D<T2>& b) {return Vector2D(a.x+b.x,a.y+b.y);}

template <class T1, class T2> inline
auto operator -(const Vector2D<T1>& a,const Vector2D<T2>& b) {return Vector2D(a.x-b.x,a.y-b.y);}

//! Dot product.
template <class T1, class T2> inline
auto operator *(const Vector2D<T1>& a,const Vector2D<T2>& b) {return a.x*b.x+a.y*b.y;}

template <class T> inline
Vector2D<T> operator *(const Vector2D<T>& a,const T F) {return Vector2D<T>(a.x*F,a.y*F);}

template <class T1, class T2> inline
auto operator *(const T1 F,const Vector2D<T2>& a) {return Vector2D(a.x*F,a.y*F);}

template <class T> inline
Vector2D<T> operator /(const Vector2D<T>& a,const T F) {return Vector2D<T>(a.x/F,a.y/F);}

template <class T> inline Vector2D<T> operator -(const Vector2D<T>& a) {return Vector2D<T>(-a.x,-a.y);}

template <class T1, class T2> inline
Vector2D<T1>& operator +=(Vector2D<T1>& a,const Vector2D<T2>& b) {a.x+=b.x;a.y+=b.y;return a;}

template <class T1, class T2> inline
Vector2D<T1>& operator -=(Vector2D<T1>& a,const Vector2D<T2>& b) {a.x-=b.x;a.y-=b.y;return a;}

//
//  Magnitude.
//
template <class T> inline T          norm     (const Vector2D<T>& a) {return (T)sqrt(a*a);}
template <class T> inline Vector2D<T> normalize(const Vector2D<T>& a) {return a/norm(a);}

//
//  Complex helpers (mirroring Vector3D).
//
template <class T> inline Vector2D<std::complex<T> > conj(const Vector2D<std::complex<T> >& v)
{ return Vector2D<std::complex<T> >(conj(v.x),conj(v.y)); }
template <class T> inline Vector2D<T> real(const Vector2D<std::complex<T> >& v) {return Vector2D<T>(real(v.x),real(v.y));}
template <class T> inline Vector2D<T> imag(const Vector2D<std::complex<T> >& v) {return Vector2D<T>(imag(v.x),imag(v.y));}
inline const Vector2D<double>& conj(const Vector2D<double>& v) {return v;}

template <class T> std::ostream& operator<<(std::ostream& os,const Vector2D<T>& v)
{
    std::streamsize prec=os.precision();
    std::streamsize wid =os.width();
    os << std::setw(0) << "(";
    os << std::setw(wid) << std::setprecision(prec) << v.x << ",";
    os << std::setw(wid) << std::setprecision(prec) << v.y << ")";
    return os;
}

} //export block
