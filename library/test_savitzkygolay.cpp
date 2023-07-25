#include <math.h>
#include <iostream>
#include <string>

#include "twave1d.h"
#include "fwave1d.h"

#include "savitzkygolay.h"

using namespace std;

// this can't have a main() and live in the library directory. until we move
// it somewhere else, change the name around a bit
int main_test_savitzkygolay(int,char**)
{
try {

    TReal1D w(1400);
	TReal1D field(w.dim());
	TReal1D deriv(w.dim());
	TReal1D deriv2(w.dim());

    TReal1D freq, w1, w2;

	double dummy;
	string sdummy;
	getline(cin,sdummy);
    for(int i=0; i<w.dim(); i++)
    {
        cin >> dummy >> field[i] >> w[i] >> deriv[i] >> deriv2[i];
    }
	w.dt(field[1]-field[0]);
	freq = w;
	w1 = w;
	w2 = w;

    SavitzkyGolay<double> sg(6,60,0);
	w = sg.filter(freq);
	sg.derivative_order(1);
	w1 = sg.filter(freq);
	sg.derivative_order(2);
	w2 = sg.filter(freq);

	
    for(int i=0; i<w.dim(); i++)
    {
//        cout << w[i] << "\t" << x[i] << "\t" << y[i] << endl;
		cout << field[i] << "\t" << freq[i] << "\t" << w[i] << "\t" << w1[i] << "\t" << w2[i] << endl;
    }

    
}
catch(const char* err)
{
    cerr << err << endl;    
}
    return 0;
}
