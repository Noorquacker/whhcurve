#include <iostream>
#include <fstream>
#include <string>

#include "tnt/tnt.h"
#include "twave1d.h"
#include "fwave1d.h"
#include "wave1d_operators.h"

#include "binarywavedata.h"
#include "textwavedata.h"
#include "lecroy_wavedata.h"

using namespace TNT;
using namespace std;

// this can't have a main() and live in the library directory. until we move
// it somewhere else, change the name around a bit
int main_filter_rowing_wav(int argc,char* argv[])
{
    try
    {
    BinaryWaveData<signed short int> wd;
    wd.openr("morningrow/3_downwind500_raw.wav");
    wd.read();

    TReal1D w((wd[0]).dim());
    Array1D<signed short int> wd_data = wd[0];
    for(unsigned i=0; i<w.dim(); i++)
        w[i] = (double) wd_data[i];
        
    w = w.subwave(22,1300000);	// get rid of the WAV header (44 bytes = 22 words) and end of wave
    w.startt(0);
    w.dt(1.0/8000);

    FReal1D freq;
    freq.isPeriodic(true);
    freq = w;

    for(unsigned i=0; i<freq.dim(); i++)
    {
        // this part to reverse the effects of the dB/dt
//        if(freq.f(i) > 25 && freq.f(i) < 50)
//            freq[i] *= 25.0/freq.f(i);
        
        // this one to filter all the high frequency garbage in two stages
        if(freq.f(i) >= 50)
            freq[i] *= exp( -1.0/50*log(2.0)/0.3*(freq.f(i)-50) );
//        if(freq.f(i) >= 60)
//            freq[i] *= exp( -1.0/60*log(2.0)/0.3*(freq.f(i)-60) );
        
        
        // and then high pass
        if(freq.f(i)==0)
            freq[i] *= 0;
        if(freq.f(i) < 20)
            freq[i] *= exp( -1.0/20*log(2.0)/0.5*(20-freq.f(i)) );
    }

    w = freq;	// reverse transform

    // apply a compressor
    for(unsigned i=0; i<w.dim(); i++)
    {
        if(w[i] > 0)
            w[i] = pow(w[i],0.33);
        else if(w[i] < 0)
            w[i] = -pow(-w[i],0.33);
    }
    
    // rescale output to 95% of full range
    double absmax = 0;
    for(unsigned i=0; i<w.dim(); i++)
    {
        if(fabs(w[i]) > absmax)
            absmax = fabs(w[i]);
    }
    for(unsigned i=0; i<w.dim(); i++)
    {
        w[i] *= (127*0.95/absmax);
    }
    

    LecroyWaveData<signed char> owd(w.dim());
    owd.openw("morningrow/downwind500_8kHz.trc");
    owd.dt(w.dt());
    owd.newHeader();
    
    // this line is a mess.
    signed char* owd_data = (Wave1D<signed char,signed char>(owd[0])).data();

    // copy data into output buffer with rounding
    for(unsigned i=0; i<w.dim(); i++)
        owd_data[i] = (signed char) (w[i]+0.5);
    
    owd.write();
    
    
    }
    catch(const char* err)
    {
        cerr << err << endl;
    }
    return 0;
}
