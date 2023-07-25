#include "wave1d_operators.h"
#include "textwavedata.h"
#include <math.h>

#define PI 3.141592654

// this can't have a main() and live in the library directory. until we move
// it somewhere else, change the name around a bit
int main_w1d_testops(int,char**)
{
    try {
    TextWaveData<double> wd("wavein.txt");
    wd.read(2);

    TReal1D u = wd[0];
    TReal1D v = wd[1];

    u[3] = PI;
    
    wd.openw("waveout.txt");
    wd.write();


    for(unsigned i=0; i<u.dim(); i++)
        cout << u[i] << "\t" << v[i] << endl;
    
    } catch(const char* errmsg)
    {
        cerr << errmsg << endl;
    }
    
    return 0;    
}
