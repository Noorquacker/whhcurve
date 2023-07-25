#include <iostream>

#include "peakfinder.h"
#include "../anyoption/anyoption.h"

REAL_FUNCTION(sin)
const double PI = 3.141592654;

using namespace std;

template <class T>
ostream& operator<<(ostream& out, const PeakData<T> p)
{
    out << "Peak: t=" << p.time() << ", a=" << p.amp() << ", d=" << p.direction() << endl;
    return out;
}
/*
int main(int argc, char** argv)
{
    try{

    double starttime = 0.000;
    double endtime = 0.0;
    bool alltime = true;
    double digitalfreq = 432e3;
    double startphase = 0.0;
    unsigned int datasize = 1000000;

    AnyOption opt;
    opt.addUsage("");
    opt.addUsage("Usage: ");
    opt.addUsage("");
    opt.addUsage(" -h  --help          Print this help ");
    opt.addUsage(" -s  --start  0.004  Start at this time (default: 0.0)");
    opt.addUsage(" -e  --end 0.008     End at this time (default: all)");
    opt.addUsage(" -f  --freq 432e8    Digital seed frequency");
    opt.addUsage(" -p  --phase 3.14159 Starting phase of wave");
    opt.addUsage(" -d  --data 1000     Number of data points in wave");
    opt.addUsage("");

    opt.setFlag("help", 'h' );
    opt.setOption( "start", 's' );
    opt.setOption( "finish", 'f' );
    opt.setOption( "freq", 'f' );
    opt.setOption( "phase", 'p' );
    opt.setOption( "data", 'd' );

    opt.processCommandArgs(argc, argv);

    if( opt.getFlag( "help" ) || opt.getFlag( 'h' ) )
    {
        opt.printUsage();
        exit(0);
    }
	if( opt.getValue( 's' ) != NULL  || opt.getValue( "start" ) != NULL  )
        starttime = atof(opt.getValue('s'));
	if( opt.getValue( 'e' ) != NULL  || opt.getValue( "end" ) != NULL  )
        endtime = atof(opt.getValue('e'));
	if( opt.getValue( 'f' ) != NULL  || opt.getValue( "freq" ) != NULL  )
        digitalfreq = atof(opt.getValue('f'));
	if( opt.getValue( 'p' ) != NULL  || opt.getValue( "phase" ) != NULL  )
        startphase = atof(opt.getValue('p'));
	if( opt.getValue( 'p' ) != NULL  || opt.getValue( "phase" ) != NULL  )
        startphase = atof(opt.getValue('p'));
	if( opt.getValue( 'd' ) != 0  || opt.getValue( "data" ) != 0  )
        datasize = atoi(opt.getValue('d'));

    TReal1D wave(datasize);
    wave.dt(5e-8);
    wave.startt(-starttime);

    TReal1D time;
    time = wave.t();
    
    wave = sin(2.0*PI*digitalfreq*time+startphase);
    
/*    for(unsigned long j=0; j<wave.dim(); j+=5)
    {
        cout << j << " ";
        for(double k=0; k<wave[j]+1.2; k+=0.05)
            cout << " ";
        cout << "*" << endl;        
    }



// cout << "wave: " << wave.subarray(0,14);

    PeakFinder<double> p(wave);
cerr << "Finding peaks from wave(" << wave.dim() << ")...";
    p.findPeaks();
cerr << "done." << endl;
    list<PeakData<double> > pl;  // Peak List
    list<PeakData<double> >::iterator pli, plj; // Peak List Iterators
    
    pl = p.getPeaks();
cerr << "Peaks found: " << pl.size() << endl;
    
    // iterate over all peaks in the list. pli will track one peak
    // ahead of plj
    pli = plj = pl.begin(); pli++;
    while(pli != pl.end())
    {
        cout << (*pli).freq(*plj) << "\t" << (0.5*(fabs((*pli).amp())+fabs((*plj).amp()))) << 
        "\t" << 
        (0.5*
        ( (*pli).time() + (*plj).time() )
        ) << endl;
        
        pli++; plj++;
    }
    
    
    
    }   // default exception handler spits it to screen.
    catch(const char* foo)
    {
        cout << foo << endl;
    }
    
    return 0;    
}
*/