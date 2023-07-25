#include <iostream>
#include <math.h>

#define TNT_NR_DIM_CHECK

#include <math.h>
#include "gaussianfit.h"
#include "tnt/tnt.h"

using namespace std;
using namespace TNT;

void nrgaussj(double **a,int n,double **b,int m);

// this can't have a main() and live in the library directory. until we move
// it somewhere else, change the name around a bit
int main_gaussianfittest(int,char**)
{
    Array1D<double> x(20);
    Array1D<double> y(20);
    Array1D<double> sig(20);
    Array1D<double> a(3);
    Array1D<double> dydt(3);
    a[0] = 7.0; a[1] = 1.0; a[2] = 1.5;
    for(int i=0; i<20; i++)
    {
        x[i] = 1.0*(i-10);
        y[i] = GaussianFit<double>::gaussianfunc(x[i],a,dydt);
        sig[i] = 1.0;
    }
    //y[4] = 3;
    a[0] = 5.0; a[1] = 3; a[2] = 3.0;
    try {
    GaussianFit<double> f(20,3);
    cout << f.fit(x,y,sig,a);
    cout << f.stddev();
    
    } catch(const char* foo) { cout << "GFT error: " << foo << endl; }
    
    return 0;
}
