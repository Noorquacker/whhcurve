#include <iostream>
#include <fstream>
#include <string>

#include "tnt/tnt.h"
#include "spline.h"
#include "twave1d.h"
#include "fwave1d.h"
#include "filters.h"
#include "wave1d_operators.h"


#include "textwavedata.h"

using namespace TNT;
using namespace std;

// this can't have a main() and live in the library directory. until we move
// it somewhere else, change the name around a bit
int main_test_spline(int argc,char* argv[])
{
    try
    {
    string sShotNumber = argv[1];
    string sFilename = string("/Users/andy/DataAnalysis/20071009_probe_vibration_analysis/C3_") + sShotNumber + string(".txt");
    string outFilename = string("/Users/andy/DataAnalysis/20071009_probe_vibration_analysis/fft_C3_") + sShotNumber + string(".txt");
    
    TextWaveData<double> wd;
    wd.openr(sFilename.c_str());
    wd.read();

    CubicSpline<double> csf(wd[2],wd[0]);
    
    TReal1D freq = csf.toWave(65536);

/*    TReal1D ysmooth = filterHP<double>(ycs,0.3,0.3);
    ysmooth = filterLP<double>(ysmooth,1,0.5);
    ysmooth *= 10.0;
*/
//    TReal1D dfreq_dfield = qdDeriv<double>(freqvsfield);

/*    cout << "s" << sShotNumber << "upField_cs\t" <<
        "s" << sShotNumber << "upFreq_cs\t" << 
        "s" << sShotNumber << "updFreqdField" << endl;
*/
    FReal1D ffreq;
    ffreq.isPeriodic(false);
    ffreq = freq;

    ofstream fout(outFilename.c_str());
    
    fout << "freqfft" << sShotNumber << "\tmagfft" << sShotNumber << endl;
    for(unsigned i=2; i<ffreq.dim(); i+=2)
    {
        fout << ffreq.f(i) << "\t" << sqrt(ffreq[i]*ffreq[i]+ffreq[i+1]*ffreq[i+1]) << endl;
    }
  
/*    freq = ffreq;
    
    cout << "tim" << sShotNumber << "\tfreq" << sShotNumber << endl;
    for(unsigned i=0; i<freq.dim(); i++)
    {
        cout << freq.t(i) << "\t" << freq[i] << endl;
    }
*/

    }//try
    catch(const char* err)
    {
        cerr << err << endl;
    }
    return 0;
}
