// File: Common/Constants.C  Some commonly used math and physics constants.
module;
#include <cmath>
  
export module Common;

export 
{
const double  Pi   = M_PI;
const double  Pi12 = sqrt(Pi);
const double  Pi32 = Pi*Pi12;
const double  Pi52 = Pi*Pi32;
const double  c_light = 137.035999139; // speed of light in atomic units

inline double Square(double x)
{
    return x*x;
}
inline double Cube  (double x)
{
    return x*x*x;
}

} //export