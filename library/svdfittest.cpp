#include <iostream>
#include <math.h>
#include "polyfit.h"
#include "tnt/tnt.h"

using namespace std;
using namespace TNT;

// this can't have a main() and live in the library directory. until we move
// it somewhere else, change the name around a bit
int main_svdfittest(int,char**)
{
    Array1D<double> x(20);
    Array1D<double> y(20);
    for(unsigned i=0; i<20; i++)
    {
        x[i] = 1.0*i;
        y[i] = 7.0*sin(i/30.0*2*3.14159);

        cout << x[i] << "\t" << y[i] << endl;        
    }
    
    PolyFit<double> f(20,3);
    try{
    f.fit(x,y);

    cout << "f(x) = " << f[0] << " + " << f[1] << "*x + " << f[2] <<  "*x**2 + " /*<< f[3] << "*x**3 + " << f[4] << "*x**4"*/ << endl;
    
    }catch(const char* foo) { cout << foo << endl; }
    
    return 0;
}
