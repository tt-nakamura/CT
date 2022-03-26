// linear interpolation
// W. H. Press, et al, "Numerical Recipes", section 3.4, 3.6

#include "nr.h"

void locate(Vec_I_DP &xx, const DP x, int &j)
{
    int ju,jm,jl;
    bool ascnd;
    
    int n=xx.size();
    jl=-1;
    ju=n;
    ascnd=(xx[n-1] >= xx[0]);
    while (ju-jl > 1) {
        jm=(ju+jl) >> 1;
        if (x >= xx[jm] == ascnd)
            jl=jm;
        else
            ju=jm;
    }
    if (x == xx[0]) j=0;
    else if (x == xx[n-1]) j=n-2;
    else j=jl;
}

double interp(double x1, const Vec_DP& x, const Vec_DP& y,
              double fill_value)
// linear interpolation in one dimension
// input:
//   x1 = evaluation point
//   x,y = data points
// return fill_value if x1 is out of x
{
    int i;
    double u;
    locate(x,x1,i);
    if(i<0 || i>=x.size()-1) return fill_value;
    u = (x1 - x[i])/(x[i+1] - x[i]);
    return (1-u)*y[i] + u*y[i+1];
}

double interp2d(double x1, double y1,
                const Vec_DP& x, const Vec_DP& y, const Mat_DP& z,
                double fill_value)
// linear interpolation in two dimensions
// input:
//   x1,y1 = evaluation point
//   x,y,z = data points
// return fill_value if (x1,y1) is out of (x,y)
{
    int i,j;
    double u,v;
    locate(x,x1,i);
    locate(y,y1,j);
    if(i<0 || i>=x.size()-1 || j<0 || j>=y.size()-1)
        return fill_value;
    u = (x1 - x[i])/(x[i+1] - x[i]);
    v = (y1 - y[j])/(y[j+1] - y[j]);
    return (1-u)*(1-v)*z[i][j] + u*(1-v)*z[i+1][j]
    + u*v*z[i+1][j+1] + (1-u)*v*z[i][j+1];
}
