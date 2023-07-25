#include <math.h>
#include <iostream>

#include "textwavedata.h"

using namespace std;

// this can't have a main() and live in the library directory. until we move
// it somewhere else, change the name around a bit
int main_test_textwavedata(int,char**)
{
	TextWaveData<double> wd(cin);
	wd.read();
	cerr << "read finished" << endl;
	
	for(int i=0; i<(wd[0].dim()); i++)
	{
		cout << wd[0][i] << endl;
	}
    
	return 0;
}
