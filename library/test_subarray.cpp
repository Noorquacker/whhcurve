#include <math.h>
#include <iostream>

#include "tnt/tnt_array1d.h"

using namespace std;
using namespace TNT;

// this can't have a main() and live in the library directory. until we move
// it somewhere else, change the name around a bit
int main_test_subarray(int,char**)
{
    Array1D<double> a(10);
    for(int i=0; i<a.dim(); i++) a[i] = 10.0*i;
  
    for(int i=0; i<a.dim(); i++)
        cout << a[i] << "\t";
    cout << endl;

    Array1D<double> b = a.subarray(11, 40);
    
    for(int i=0; i<b.dim(); i++)
        cout << b[i] << "\t";
    cout << endl;

    return 0;
}

