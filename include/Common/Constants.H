// File: Common/Constants.H  Some commonly used math and physics constants.
#ifndef _Common_Constants_H
#define _Common_Constants_H

extern const double  Pi  ;
extern const double  Pi12;
extern const double  Pi32;
extern const double  Pi52;

extern const double  c_light; // speed of light in atomic units

inline double Square(double x)
{
    return x*x;
}
inline double Cube  (double x)
{
    return x*x*x;
}

#endif //_Common_Constants_H
