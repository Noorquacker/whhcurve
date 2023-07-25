#ifndef _PEAKFINDER_H
#define _PEAKFINDER_H

// STL linked lists
#include <list>
#include "peakdata.h"
#include "../tntnr.h"
#include "../wave1d.h"
#include "../twave1d.h"
#include "../fwave1d.h"
#include "../wave1d_operators.h"
#include "../tnt/tnt.h"
#include "../polyfit.h"

using namespace std;

extern const char* FITPEAK_INDEX_OUT_OF_BOUNDS;

/** PeakFinder
Implements algorithm to find (in floating point) peaks of a given wave. DO NOT
PASS BY VALUE. Destroying copies of object will attempt to free memory multiple
times.
*/
template <class T>
class PeakFinder
{
    public:
    PeakFinder(const TWave1D<T,T>& w) { init(w); }
    ~PeakFinder()
    {
        delete polyfitter;
    }
    
    void findPeaks()
    {
        unsigned long n;
        T freq;
        typename list<PeakData<T> >::const_reverse_iterator it;
        PeakData<T> p1,p2;

        try
        {
            n = addForwardPeak(INITIAL_POINTS,INITIAL_POINTS);
    
cout << "addForwardPeak() returned " << n << endl;
    
            bool first = true;
        
            // IMPORTANT! Exit condition is taken care of by the try/catch
            // exception blocks. We don't know exactly *which* function call
            // might generate the signal to get out, so exit condition is
            // the thrown exception of the proper type and value. This isn't
            // a standard control structure, but it's convenient. Look up
            // Exception Handling in a C++ reference.
            while(1)
            {
                if(first)
                {
cout << "calling addForwardPeak(" << n << ") from first" << endl;
                    n = addForwardPeak(n,INITIAL_POINTS);
                    first = false;

                    it = peaks.rbegin();
                    p1 = (*it);
//cout << "p1: " << p1 << endl;                
                    it++;
                    p2 = (*it);
//cout << "p2: " << p2 << endl;                
                    freq = p1.freq(p2);
cout << "Freq from first loop: " << freq << endl;
                }
                else
//                    n = addForwardPeak(n,fitSize(freq));
                    addForwardPeak(n,fitSize(freq));
                                
                // recalculate frequency between last two peaks.
                it = peaks.rbegin();
                p1 = (*it);
//cout << "p1: " << p1 << endl;                
                it++;
                p2 = (*it);
//cout << "p2: " << p2 << endl;                
                freq = p1.freq(p2);
cout << "Freq from subsequent loop: " << freq << endl;
                
                n = wave.index(2*p1.time()-p2.time())-fitSize(freq)/2;
            }
        }
        catch(const char* err)
        {
            // technique: if it's the error we expect, everything is good,
            // and the unroll mechanism has gotten us out of all functions
            // and loops.
            //
            // if the error isn't what we expect, something is wrong, so
            // rethrow it.
            if(err != FITPEAK_INDEX_OUT_OF_BOUNDS)
                throw err;
        }
        // through the virtues of the exception-unrolling mechanism, we're
        // already out of our loops and done!
    }
    
    /**
    Adds the next peak found in the data to the list of peaks.
    */
    unsigned long addForwardPeak(unsigned long start, const unsigned fitpoints)
    {
cout << "Starting addForwardPeak() with start=" << start << ", fitpoints=" << fitpoints << endl;        
        // test current starting point, then move forward half the
        // width until we find peak
        PeakData<T> p, p2;
        unsigned long newstart;

        // first, find it roughly
        while(p.direction()==0)
        {
cout << "Looping, looking for a maximum." << endl;

            // first check for end of array (yes, this happens a few times in
            // the algo to cover all possible runover cases)
            if(wave.dim() <= start+fitpoints)
                throw FITPEAK_INDEX_OUT_OF_BOUNDS;
            
            // speedup: if the max is at an endpoint, ignore this bit of data            
            if(!extBetweenEndPoints(wave.subarray(start,start+fitpoints-1)))
            {
                start += fitpoints - 1;
                continue;
            }
cout << "Checking for a peak at " << start << ", " <<  fitpoints << endl;
            
            p = fitPeak((*polyfitter),start,fitpoints);
            start += fitpoints/2;    // increment in anticipation of another
        }
        start -= fitpoints/2;    // put it back after the last one

        // put new starting point at found peak and find peak from symmetric
        // data
cout << "    p.time()=" << p.time() << endl;
        start = wave.index(p.time()) - fitpoints/2;
cout << "    start=" << start << endl;         
        p = fitPeak((*polyfitter),start,fitpoints);
        
        // Alright, this gives us one peak using so-called reasonably symmetric data. Now get another different
        // one:
        if(p.time() > wave.t(start+fitpoints/2)) // if peak is AFTER our window-center
        {
            newstart = wave.index(p.time());	// first set new starting point to our peak time.
            if(newstart == start+fitpoints/2)	// but if it's the same as our old one, increment it
                newstart++;
        }
        else if(p.time() < wave.t(start+fitpoints/2)) // if peak is BEFORE our window center
        {
            newstart = wave.index(p.time());	// first set new starting point to our peak time.
            if(newstart == start+fitpoints/2)	// but if it's the same as our old one, decrement it
                newstart--;
        }
        newstart -= fitpoints/2;	// convert from anticipated peak to starting point
//cout << "Peak 1: " << p << endl;        
        
        p2 = fitPeak((*polyfitter),newstart,fitpoints);
//cout << "Peak 2: " << p2 << endl; 

        // Weight and average them according to the distance from the peak to the center of the window used.
        // Note that the weights are SWITCHED, so the closer point gets the higher weighting.
        T weight1, weight2;
        weight2 = abs<T>(wave.t(start+fitpoints/2)-p.time());
        weight1 = abs<T>(wave.t(newstart+fitpoints/2)-p2.time());

//cout << "Weight 1: " << weight1 << ", weight 2: " << weight2 << endl;

        // weighted average
        PeakData<T> avgpeak;
        avgpeak.time((p.time()*weight1+p2.time()*weight2)/(weight1+weight2));
        avgpeak.amp((p.amp()*weight1+p2.amp()*weight2)/(weight1+weight2));

//cout << "Averaged peak: " << avgpeak << endl;

        // add to list
cout << "Adding peak: " << p << endl;        
        peaks.push_back(avgpeak);
        
        // give back an index at beginning of the next window.
cout << "addForwardPeak() returns " << start+fitpoints/2+fitpoints/4 << endl;
        return wave.index(avgpeak.time())+fitpoints/4;
    }

// haven't implemented this yet
//    unsigned long addNearestPeak(unsigned long start, unsigned fitpoints);
    
    const list<PeakData<T> >& getPeaks() { return peaks; }
    void clearPeaks() { peaks.clear(); }
    
    protected:
    unsigned INITIAL_POINTS;
    unsigned MAX_POINTS;
    unsigned FIT_FUNCTIONS;

    Array1D<T> timearray;
    TWave1D<T,T> wave;
    list<PeakData<T> > peaks;
    
    PolyFit<T> *polyfitter;
    PolyFit<T> *nearestfitter;
    
    /**
    Standard Andy callable constructor
    */
    void init(const TWave1D<T,T>& w)
    {
        INITIAL_POINTS = 7;
        MAX_POINTS = 1001;
        FIT_FUNCTIONS = 4;

        wave = w;
        
        // timearray always starts at zero, so fits aren't
        // horribly ill-conditioned due to a large constant factor
        timearray = Array1D<T>(MAX_POINTS);
        for(unsigned i=0; i< timearray.dim(); i++)
            timearray[i] = i*w.dt();

        polyfitter = new PolyFit<T>(MAX_POINTS,FIT_FUNCTIONS);
    }

    /**
    Fits one peak to the data using the given fitting object and data window.
    If unsuccessful, returns a null peak (peak.direction() == 0).
    */
    PeakData<T> fitPeak(PolyFit<T>& fitter, unsigned long start, unsigned fitpoints)
    {
//cout << "fitPeak: start=" << start << ", fitpoints=" << fitpoints << endl;

        // first check for end of array
        if(wave.dim() <= start+fitpoints)
            throw FITPEAK_INDEX_OUT_OF_BOUNDS;

        // need to put a check here to make sure we don't let ourselves
        // ask for a fit that's bigger than our fitter can handle
        Array1D<T> x = timearray.subarray(0,fitpoints-1);
        Array1D<T> y = wave.subarray(start,start+fitpoints-1);

//cout << "x: " << x << endl;        
//cout << "y: " << y << endl;        
        fitter.fit(x,y);
//cout << "Fit: " << fitter.params() << endl;
        PeakData<T> p = polyExtremum(fitter.params(),start,fitpoints);

cout << "fitPeak() returning " << p << "...whose real time is " << (p.time()+wave.dt()*start) << endl;

        // finally, shift time back to full wave scale            
        p.time(p.time()+wave.dt()*start);
        
        return p;
        
    }
    
    /**
    Computes the extremum (min/max) of a polynomial nearest to the guess g
    using bisection, and returns result as a PeakData<T> peak. If peak falls
    outside interval, returns a null PeakData<T>.
    */    
    PeakData<T> polyExtremum(const Array1D<T>& a, unsigned long start, unsigned size)
    {
        const T ABS_TOL = 1e-10;

        // technique here is to take the derivative of the
        // polynomial, then find the roots
        
        // polynomial derivative
        Array1D<T> dadt(a.dim()-1);
        for(unsigned long i=0; i<dadt.dim(); i++)
        {
            dadt[i] = (i+1)*a[i+1];
        }
//cout << "***polyExtremum:***" << endl;
//cout << "a: " << a << endl;
//cout << "dadt: " << dadt << endl;
        
        PeakData<T> q;

        T xlow = 0;
        T low = fpoly(dadt,xlow);
        T xhigh = wave.dt()*size;
        T high = fpoly(dadt,xhigh);

//cout << "xlow=" << xlow << ", xhigh=" << xhigh << endl;
//cout << "low=" << low << ", high=" << high << endl;

        // find the peak of a quadratic quickly without searching
        if(dadt.dim() == 2)
        {
            q.time(-dadt[0]/dadt[1]);
            q.amp(fpoly(a,q.time()));
            q.direction(dadt[1]==0 ? 0 : dadt[1]>0 ? -1 : 1);
//cout << "Quad peak: " << q.time() << endl;
            if(q.time() > xlow && q.time() < xhigh)
                return q;
            else
                return PeakData<T>(0,0,0);
        }
        
        // unlikely solutions first
        if(abs<T>(low) <= ABS_TOL)
        {
            q.time(xlow);
            q.amp(fpoly(a,xlow));
            if(high>ABS_TOL)
                q.direction(1);
            else if(high < -ABS_TOL)
                q.direction(-1);
            return q;
        }
        if(abs<T>(high) <= ABS_TOL)
        {
            q.time(xhigh);
            q.amp(fpoly(a,xhigh));
            if(low>ABS_TOL)
                q.direction(-1);
            else if(low < -ABS_TOL)
                q.direction(1);
            return q;
        }
//cout << "dadt: " << dadt;
        q.time(bisectRoot(dadt,xlow,xhigh,low,high,ABS_TOL));
cout << "polyExtremum: bisectRoot returns " << q.time() << endl;
        q.amp(fpoly(a,q.time()));
        if(low>ABS_TOL)
            q.direction(-1);
        else if(low < -ABS_TOL)
            q.direction(1);
        
        return q;
    }
    /**
    Bisection root-finder:
    xhigh and xlow define bracket, high and low are the function values
    at xhigh and xlow respectively. high and low MUST have opposite sign.
    */
    T bisectRoot(const Array1D<T>& a, T xlow, T xhigh, T low, T high, T abstol)
    {
//cout << "**bisectRoot()** with abstol " << abstol << endl;
//cout << "x: " << xlow << ", " << xhigh << endl;
//cout << "y: " << low << ", " << high << endl;

        T xmid = (xlow+xhigh)/2;
        T mid = fpoly(a,xmid);
//cout << "xmid: " << xmid << ", mid: " << mid << endl;
        
        // trivial case
        if(abs<T>(mid) <= abstol || (xhigh-xlow) <= abstol)
            return xmid;
        if(mid * low < 0)
            bisectRoot(a,xlow,xmid,low,mid,abstol);
        if(mid * high < 0)
            bisectRoot(a,xmid,xhigh,mid,high,abstol);
    }
    /**
    Computes the value of a polynomial described by array a at point x. There is a much better algorithm for
    this that is linear in polynomial order.
    */
    T fpoly(const Array1D<T>& a, const T x)
    {
//cout << "fpoly(a,x): x=" << x << " a:" << a << endl;
        T q = 0.0;
        for(unsigned i=0; i<a.dim(); i++)
        {
            T z = a[i];
            for(unsigned j=0; j<i; j++)
                z *= x;
            q += z;             
        }
//cout << "fpoly returns " << q << endl;
        return q;
    }
    /**
    Computes the desired size of the fit from a frequency guess. Currently
    aims to use the top 1/3 of each peak, so 1/6 of the period.
    */
    unsigned fitSize(const T freq)
    {
        // number of points per period /6, rounded
        //unsigned int p = (unsigned int) (1.0 / (wave.dt() * freq) / 6.0 + 0.5);
        unsigned int p = (unsigned int) (1.0 / (wave.dt() * freq) / 4.0 + 0.5);
        
        // round p up to next odd integer
        if(p%2 == 0) p++;
        
        return p;
    }

    bool extBetweenEndPoints(const Array1D<T> a)
    {
        T max,min;
        unsigned long mini,maxi;
        
        max = min = a[0];
        maxi = mini = 0;
        for(unsigned long i=0; i<a.dim(); i++)
        {

            if(min > a[i])
            {
                min = a[i];
                mini = i;
            }
            if(max < a[i])
            {
                max = a[i];
                maxi = i;
            }
        }
        if( (mini > 0 && mini < a.dim()-1) || (maxi > 0 && maxi < a.dim()-1) )
        {
            return true;
        }
        // else
        return false;
    }
};

#endif //_PEAKFINDER_H
