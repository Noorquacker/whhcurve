#include <math.h>
#include <iostream>

#include "twave1d.h"
#include "fwave1d.h"

using namespace std;

double hp(double f);
double lp(double f);

REAL_FUNCTION(sin)
REAL_FUNCTION(cos)
REAL_FUNCTION(fabs)
REAL_FUNCTION(hp)
REAL_FUNCTION(lp)

// this can't have a main() and live in the library directory. until we move
// it somewhere else, change the name around a bit
int main_test_wave1d(int,char**)
{
    double PI = 3.141592654;
    TReal1D x(65536);
    x.mint(-0.5);
    x.maxt(0.5);
    x = x.t();
    TReal1D sin2x = fabs(x);
    FReal1D fft;
    fft = sin2x;
    FReal1D f = fft.f();

    for(int i=0;i<fft.dim();i+=2)
        fft[i] *= lp(f[i]);


/*    cout << f[0] << "\t" << abs(fft[0]) << endl;
    for(int i=2;i<fft.dim();i+=2)
    {
        cout << f[i] << "\t" << sqrt(fft[i]*fft[i]+fft[i+1]*fft[i+1]) << endl;
    }
    cout << f[1] << "\t" << abs(fft[1]) << endl;
*/

    TReal1D x2;
    x2 = fft;
    
    for(int i=0;i<x.dim();i++)
        cout << x[i] << "\t" << sin2x[i] << "\t" << x2[i] << endl;

    cout << endl;


    return 0;
}

double hp(double f)
{
    return f>3 ? 1 : 0;
}
double lp(double f)
{
    return (f<3 ? 1 : 0);
}

