#ifndef _DEMOD_H
#define _DEMOD_H

// comment out to do without multithreading with OpenMP
//#define _OPENMP

#ifdef _OPENMP
	#include <omp.h>
#endif

#ifdef __OBJC__ 
	#include "fwave1d.h"
	#include "twave1d.h"

	#include <algorithm>
	#include "tntnr.h"
	#include "textwavedata.h"
	#include "lecroy_wavedata.h"
	#include "savitzkygolay.h"
#else
	#include "textwavedata.h"
	#include "lecroy_wavedata.h"
	#include "twave1d.h"
	#include "fwave1d.h"
	#include "savitzkygolay.h"
	#include "tntnr.h"
#endif

#define	FOURIER_DEMOD_UNINITIALIZED	0
#define FOURIER_DEMOD_INITIALIZED	1
#define FOURIER_DEMOD_UNINITIALIZED_ERROR	"FourierDemodulator (or a subclass) has not been initialized."

typedef struct
{
	double mag;
	double a;
	double b;
} fint_ret_t;

typedef struct
{
	double freq;
	double mag;
} find_fft_peak_ret_t;

// return structure for brent_four
typedef struct 
{
	double freq;
	double a;
	double b;
} brent_four_ret_t;

typedef struct
{
	TWave1D<double,double> time;
	TWave1D<double,double> freq;
	TWave1D<double,double> amp;
	TWave1D<double,double> phase;
	TWave1D<double,double> jitter;
} tfapj_demodulation_t;

//Haven't gotten it to run without this here at the moment. -Steve
double percent_done;
bool cancel_demod = false;

template <class T, class Tbigger>
class FourierDemodulator
{
	public:
	FourierDemodulator(const TWave1D<T,double> wave_, const int window_size_) { create(wave_,window_size_); }
	FourierDemodulator(const FourierDemodulator<T,Tbigger>& old) { create(old); }
	FourierDemodulator() {instance_state = FOURIER_DEMOD_UNINITIALIZED;}
	
	void create(const TWave1D<T,double> wave_, int window_size_)
	{
		wave = wave_;
		windowsize(window_size_);
		percent_done = 0;
		cancel_demod = false;
		
		jitterthreshold(defaultjitterthreshold());
		jitteraveraging(defaultjitteraveraging());
		correctionlevel(defaultcorrectionlevel());
		
		definewindow();
		instance_state = FOURIER_DEMOD_INITIALIZED;
	}
	
	void create(const FourierDemodulator<T,Tbigger>& old)
	{
		if(old.instance_state == FOURIER_DEMOD_UNINITIALIZED)
		{
			instance_state = FOURIER_DEMOD_UNINITIALIZED;
			return;
		}
		
		// old is good. copy it.
		wave = old.wave;
		windowsize(old.windowsize());
		win = old.win;
		jitterthreshold(old.jitterthreshold());
		jitteraveraging(old.jitteraveraging());
		correctionlevel(defaultcorrectionlevel());

		percent_done = 0;
		cancel_demod = false;
		instance_state = FOURIER_DEMOD_INITIALIZED;
	}
	
	FourierDemodulator<T,Tbigger>& operator = (const FourierDemodulator<T,Tbigger>& old) { create(old); return (*this); }
	
	int windowsize(const int s) { return window_size = (s > 0 ? s : defaultwindowsize()); }
	int windowsize() const { return window_size; }
	int defaultwindowsize() { return wave.dim() / 1024; }

	int jitteraveraging(const int s) { return JITTER_AVERAGING_LENGTH = (s > 0 ? s : defaultjitteraveraging()); }
	int jitteraveraging() const { return JITTER_AVERAGING_LENGTH; }
	
	double jitterthreshold(const double j) { return JITTER_THRESHOLD = abs<double>(j); }
	double jitterthreshold() const { return JITTER_THRESHOLD; }

	int correctionlevel(int cr) { return CORRECTION_LEVEL = (cr >=0 && cr <= 2 ? cr : CORRECTION_LEVEL); }
	int correctionlevel() const { return CORRECTION_LEVEL; }
	int defaultcorrectionlevel() const { return 2; }

	// rationale behind jitter analysis is tricky:
	// N is point spacing (hann window: window_size/2)
	// n is number of full oscillations
	// dt is sampling interval
	// 1.0/2/dt is nyquist frequency (worst case scenario)
	// f0 is frequency
	// df is maximum allowed jitter
	//
	// n = N/(2*dt*f0)
	// likewise, n-0.5 = N/(2*dt*(f0 + df))
	//
	// combining,
	// N/(2*dt*f0) - 0.5 = N/(2*dt*(f0 + df)
	//
	// assuming the nyquist frequency for f0,
	// N - 0.5 = N / ( 1 + 2*dt*df )
	// 1 + 2*dt*df = N/(N-0.5)
	// 2*dt*df = N/(N-0.5) - 1
	// df = 1/2/dt * (N/(N-0.5) - 1)
	// df = 1/2/dt * (window_size/(window_size-1)-1)
	//
	// But the data points are not independent with a windowed wave. We see half the data our
	// neighbor sees each point. So as an approximation, cut our df in half. We have already
	// made a very conservative choice in using the Nyquist frequency.
	double defaultjitterthreshold() const { return 1.0/(4*wave.dt())*(window_size/(double)(window_size-1)-1); }
	int defaultjitteraveraging() const { return 9; }

	
	// prime movers
	virtual tfapj_demodulation_t demod();	// runs the algorithm with current options
	
	// track our progress
	double progress_bar() const { return percent_done; }
	void cancel() { cancel_demod = true; }
	void reset_cancellation() { cancel_demod = false; }
	
	// some constants
	static const double PI;
	
	protected:
	TWave1D<T,double> wave;

	TWave1D<double,double> freqs;
	TWave1D<double,double> amps;
	TWave1D<double,double> phases;
	TWave1D<double,double> jitters;
	TWave1D<double,double> jitters_smth;
	TWave1D<int,double> jitter_free_points;
	
	Array1D<Tbigger> win;

	//Declared up top
	//bool cancel_demod;
	
	// window function
	TWave1D<Tbigger,double> applywindow(const TWave1D<T,double>& data);
	
	// brent's method
	brent_four_ret_t brent_four(double ax, double bx, double cx, double tol, const TWave1D<Tbigger,double>& data);
	
	// fixed and floating point fourier integrals
	virtual fint_ret_t four(const TWave1D<Tbigger,double>& data, double freq) = 0;
	
	// identify peak in fft
	virtual find_fft_peak_ret_t find_fft_peak(const TWave1D<Tbigger,double>& tw_T);
	
	virtual void setones() = 0;
	
	virtual void definewindow()
	{
		win = Array1D<Tbigger>(window_size);
		for(int i=0; i<win.dim(); i++)
		{
//			win[i] = (Tbigger)(2 * MATH_ONE * sin(i*PI/win.dim())*sin(i*PI/win.dim()));
			double a0 = 0.54;
			double a1 = 0.46;
//			double a0 = 0.50;
//			double a1 = 0.50;
			double a2 = 0;
			double a3 = 0;
			win[i] = (Tbigger)(MATH_ONE*a0 - 
				MATH_ONE*a1*cos(2*PI*i/(win.dim())) +
				MATH_ONE*a2*cos(4*PI*i/(win.dim())) -
				MATH_ONE*a3*cos(6*PI*i/(win.dim()))
				);
			
			// use if you need a box window for testing
			//				win[i] = factor;
		}
	}
	
	// variables that control our behavior
	double JITTER_THRESHOLD;	// allowable frequency jitter from one point to the next.
	int window_size;			// averaging length for one data point
	Tbigger MATH_ONE;		// value of 1.0 in fixed point math. automatically set to be slightly higher
								// resolution than the incoming data
	int JITTER_AVERAGING_LENGTH;	// number of continuous points that must be jitter free,
									// also used as the averaging length for applying the
									// phase-coherence correction to the frequency
	int CORRECTION_LEVEL;	// how much of our machinery for dropping jittery points and phase correction
							// should we use? (0=none, but do report phase jitter, 1=drop points where we
							// can't determine phase coherence, 2=correct the frequency to maintain phase
							// coherence, smoothing by JITTER_AVERAGING_LENGTH)
	
	int instance_state;
};

template <class T, class Tbigger>
const double FourierDemodulator<T,Tbigger>::PI = 3.14159265358979323846264338327950288419716939937510;

template <class T, class Tbigger>
TWave1D<Tbigger,double> FourierDemodulator<T,Tbigger>::applywindow(const TWave1D<T,double>& data)
{
	TWave1D<Tbigger,double> wout(data.dim());
	wout.dt(data.dt());
	wout.startt(data.startt());
	
	Tbigger sum = 0;
	for(int i=0; i<wout.dim(); i++)
	{
		wout[i] = win[i]*data[i];
		sum += wout[i];
	}
	
	// subtract common mode
	sum /= data.dim();
	for(int i=0; i<wout.dim(); i++)
	{
		wout[i] -= sum;
	}
	
	return wout;
}


template <class T, class Tbigger>
find_fft_peak_ret_t FourierDemodulator<T,Tbigger>::find_fft_peak(const TWave1D<Tbigger,double>& tw_T)
{
	TWave1D<double,double> tw(tw_T.dim());
	tw.dt(tw_T.dt());
	tw.startt(tw_T.startt());
	
	for(int i=0; i<tw.dim(); i++)
	{
		tw[i] = tw_T[i];
	}
	
	FWave1D<double,double> fw = tw;
//cerr << "first few lines of fw: " << fw[0] << "\t" << fw[1] << "\t" << fw[2] << "\t" << fw[3] << "\t" << fw[4] << "\t" << fw[5] << endl;
	
	double trial;
	double max = 0;
	double maxf = 0;
	for(int i=MAX(fw.dim()/50,6); i<fw.dim()-1; i+=2)	// don't start at exactly zero, because window might give a false peak
	{
		trial = fw[i]*fw[i]+fw[i+1]*fw[i+1];
		if(trial > max)
		{
			max = trial;
			maxf = fw.f(i);
		}
	}
//cerr << "find_fft_peak: " << maxf << " : " << max << endl;
	find_fft_peak_ret_t r = {maxf, max};
	return r;
}


// brent's method for finding the peak of a fourier transform. last three params are just pass-through
template <class T, class Tbigger>
brent_four_ret_t FourierDemodulator<T,Tbigger>::brent_four(double ax, double bx, double cx, double tol, const TWave1D<Tbigger,double>& data)
{
	const int ITMAX = 100;
	const double CGOLD = 0.3819660;
	const double ZEPS = 1.0e-10;
//cout << "Starting Brent with ax=" << ax << ", bx=" << bx << ", cx=" << cx << "\n";
	
    int iter;
    double a,b,d,etemp, p,q,r,tol1,tol2,u,v,w,x,xm;
	fint_ret_t fu,fv,fw,fx;
    d = 0;	// otherwise it starts uninitialized. no idea if that's possibly bad, but i'm sick of the compiler warning
	double e=0.0;	// This will be the distance moved on the step before last.
    
    a=(ax < cx ? ax : cx);	// a and b must be in ascending order, but input abscissas need not be.
    b=(ax > cx ? ax : cx);
    x=w=v=bx;	// Initializations...
    
		fw=fv=fx=four(data,x);

//cout << data.startt()+data.dim()/2*data.dt() << "\t" << x << "\t" << (double)fw.mag << endl;
	    
    
    
    for (iter=1;iter<=ITMAX;iter++) {	// main program loop
        xm=0.5*(a+b);
        tol2=2.0*(tol1=tol*abs<double>(x)+ZEPS);
        if (abs<double>(x-xm) <= (tol2-0.5*(b-a))) {	// test for done here.
			//cout << "Apparently, Brent thinks we're done already. tol1=" << tol1 << ", tol2=" << tol2 << "\n";            
			// *xmin=x; // from original
			// return fx;	// from original
            brent_four_ret_t ret = {x, fx.a, fx.b};
//cout << endl;
			return ret;
        }
        if (abs<double>(e) > tol1) {	// construct a trial parabolic fit.
            r=(x-w)*(fv.mag-fx.mag);
            q=(x-v)*(fw.mag-fx.mag);
            p=(x-v)*q-(x-w)*r;
            q=2.0*(q-r);
            if (q > 0.0) p = -p;
            q=abs<double>(q);
            etemp=e;
            e=d;
            if (abs<double>(p) >= abs<double>(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
                d=CGOLD*(e=(x >= xm ? a-x : b-x));
            /* The above conditions determine the acceptability of the parabolic fit. Here we take the 
			 golden section into the larger of the two segments. */
            else {
                d=p/q;	// Take the parabolic step.
                u=x+d;
                if (u-a < tol2 || b-u < tol2)
                    d=SIGN<double>(tol1,xm-x);
            }
        } else {
            d = CGOLD*(e=(x >= xm ? a-x : b-x));
        }
        u=(abs<double>(d) >= tol1 ? x+d : x+ SIGN<double>(tol1,d));
        // This is the one function evaluation per iteration
		fu = four(data,u);

//		cout << data.startt()+data.dim()/2*data.dt() << "\t" << u << "\t" << (double)fu.mag << endl;
        
        if (fu.mag > fx.mag) {	// Now decide what to do with our function evaluation.
            if (u >= x) a=x; else b=x;
            SHFT<double>(v,w,x,u);
            SHFT<fint_ret_t>(fv,fw,fx,fu);
        } else {
            if (u < x) a=u; else b=u;
            if (fu.mag > fw.mag || w == x) {
                v=w;
                w=u;
                fv=fw;
                fw=fu;
            } else if (fu.mag > fv.mag || v== x || v == w) {
                v=u;
                fv=fu;
            }
        }
    }	// done with housekeeping. back for another iteration.
    throw "Too many iterations in brent.";
	brent_four_ret_t ret = {0.0, 0.0, 0.0};
	return ret;
}

template <class T, class Tbigger>
tfapj_demodulation_t FourierDemodulator<T,Tbigger>::demod()
{
	if(instance_state == FOURIER_DEMOD_UNINITIALIZED)
	{
		throw FOURIER_DEMOD_UNINITIALIZED_ERROR;
	}
	// since this is virtual, we can't call it in the constructor
	setones();
	definewindow();

	freqs = TWave1D<double,double>(wave.dim()*2/window_size-1);
	amps = TWave1D<double,double>(freqs.dim());
	phases = TWave1D<double,double>(freqs.dim());
	jitters = TWave1D<double,double>(freqs.dim());
	jitters_smth = TWave1D<double,double>(freqs.dim());
	jitter_free_points = TWave1D<int,double>(freqs.dim());
	
	freqs.dt(window_size*wave.dt()/2);
	freqs.startt(wave.startt()+window_size*wave.dt());
	
	amps.dt(freqs.dt());
	phases.dt(freqs.dt());
	jitters.dt(freqs.dt());
	jitters_smth.dt(freqs.dt());
	jitter_free_points.dt(freqs.dt());
	
	amps.startt(freqs.startt());
	phases.dt(freqs.dt());
	jitters.startt(freqs.startt());
	jitters_smth.startt(freqs.startt());
	jitter_free_points.startt(freqs.startt());
	
//cerr << "Starting demodulation" << endl;
	#ifdef _OPENMP
		#pragma omp parallel for
	#endif
	for(int j=0; j<freqs.dim(); j++)
	{
//#pragma omp critical
//{ cerr << "j=" << j << endl; }

		if(cancel_demod)
		{
		    continue;
		}
		TWave1D<Tbigger,double> windowed_wave = applywindow(wave.subwave(j*window_size/2,window_size));
		find_fft_peak_ret_t fft_peak = find_fft_peak(windowed_wave);

		double fmid = fft_peak.freq;
		double fmin = fmid - 1.0/wave.dt()/window_size;
		double fmax = fmid + 1.0/wave.dt()/window_size;
		//cerr << "(" << j << "/" << outwave.dim() << ") find_fft_peak: " << fft_peak.freq << " : " << fft_peak.magsq << endl;
		if(fmin <= 0)
		{
			fmin = 0;
			fmid = (fmin+fmax)/2;
		}
		if(fmax >= 0.5/wave.dt())
		{
			fmax = 0.5/wave.dt();
			fmid = (fmin+fmax)/2;
		}
		
		brent_four_ret_t r = brent_four(fmin,fmid,fmax, 1e-12, windowed_wave);
		freqs[j] = r.freq;
		amps[j] = sqrt(r.a*r.a+r.b*r.b)/MATH_ONE;
		phases[j] = (r.a == 0 ? SIGN<double>(PI,r.b) : atan2(r.b,r.a));
		
		percent_done = 95.0*j/freqs.dim();
	}
	//cerr << "Multithreading complete... postprocessing" << endl;

	// Postprocessing options:
	// "correction level"
	//	0) always calculate jitter
	//	1) can also throw out jittery points
	//	2) frequency correction to maintain phase-lock
	//
	//
	
	// strategy here is to try and determine if this point is phase coherent with its neighbors. if it's not,
	// we have no use for it.
	int num_low_jitter_points_in_row = 1;	// first point always counts
	for(int j=0; j<freqs.dim()-1; j++)
	{
		if(j>0)
		{
			// first find out a prediction of our phase shift from the frequencies and window interval
			const double wholel = (freqs[j]+freqs[j-1])/2*(freqs.t(j)-freqs.t(j-1));
			const double fpartl = (wholel - (int)wholel)*2*PI;
			const double wholer = (freqs[j+1]+freqs[j])/2*(freqs.t(j+1)-freqs.t(j));
			const double fpartr = (wholer - (int)wholer)*2*PI;
			
			// now see what the phase information itself tells us. use the interval [0:2pi]
			double pdiffl = phases[j-1]-phases[j];
			double pdiffr = phases[j]-phases[j+1];
			if(pdiffl<0)
				pdiffl+=2*PI;
			if(pdiffr<0)
				pdiffr+=2*PI;
			//			cout << "\t" << pdiff;
			//			cout << "\t" << fpart;
			
			// this section effectively rounds to the nearest oscillation. if that's a bad
			// assumption, we've clearly lost phase coherence anyway
			double jitterl = fpartl - pdiffl;
			if(jitterl < -PI)
				jitterl += 2*PI;
			if(jitterl >= PI)
				jitterl -= 2*PI;
			double jitterr = fpartr - pdiffr;
			if(jitterr < -PI)
				jitterr += 2*PI;
			if(jitterr >= PI)
				jitterr -= 2*PI;

			// this used to convert the phase difference to a frequency correction
			double jitter = (jitterl+jitterr)/2;
			

			jitters[j] = jitter;
			jitters_smth[j] = jitter;
			// check for jitter and frequencies consistent with phase coherence
			if(abs<double>(jitterr - jitterl) < (PI/2) &&
			   abs<double>(freqs[j]-freqs[j-1]) < JITTER_THRESHOLD &&
			   abs<double>(freqs[j+1]-freqs[j]) < JITTER_THRESHOLD)
			{
				num_low_jitter_points_in_row++;
			}
			else
			{
				num_low_jitter_points_in_row = 0;
			}
			jitter_free_points[j] = num_low_jitter_points_in_row;
		}
		else	// first/last points
		{
			jitter_free_points[0] = 1;	// end points are always jitter-free
			jitter_free_points[jitter_free_points.dim()-1] = 1;
		}
		// from 95% to 97%
		percent_done = 95.0 + 2.0*j/freqs.dim();
	}

	if(cancel_demod)
	{
	    tfapj_demodulation_t null_demod;
	    return null_demod;
	}

	// now iterate over all the points again, throw out the bad ones, and smooth the phase error (jitter) wave
	// before we apply the correction to freq

	SavitzkyGolay<double> sg(2,JITTER_AVERAGING_LENGTH);
	int good_point_counter = 0;
	for(int j=0; j<jitters_smth.dim(); j++)
	{
		if(jitter_free_points[j] == 1)
		{
			// is the minimally last jitter-free point of this set within the wave?
			if(j+JITTER_AVERAGING_LENGTH <= jitters.dim())
			{
				if(jitter_free_points[j+JITTER_AVERAGING_LENGTH-1] == JITTER_AVERAGING_LENGTH)
				{
					int i;
					for(i=j; jitter_free_points[i] >= 1 && i<jitters.dim(); i++)
					{
						// for each jitter free point of the set
						good_point_counter++;
					}
					i--;	// we overshot by one. put it back.
					// do the filtering on this coherent set. = operator is a reference,
					// so we'll be operating on this portion of jitters[]
					TWave1D<double,double> coherent_set = jitters_smth.subwave(j,i-j);
					
					// smooth the jitter data using a Savitzky-Golay filter
					TWave1D<double,double> temp = sg.filter(coherent_set);
					
					// copy back over from temp to coherent_set (and therefore jitters_smth)
					for(int k=0; k<temp.dim(); k++)
						coherent_set[k] = temp[k];
					
					// advance j to the next point
					j = i+1;
				}
				else
					jitter_free_points[j] = 0;
			}
			else
				jitter_free_points[j] = 0;
				
		}
		else
			jitter_free_points[j] = 0;
	
		// from 97% to 98%
		percent_done = 97.0 + 2.0*j/freqs.dim();
	}
	
	// if correction level is zero, we want all the data, not just the good stuff
	if(correctionlevel() < 1)
		good_point_counter = freqs.dim();
		
	// last loop. make the output wave with jitter corrections to the frequency and return it
	tfapj_demodulation_t r;
	r.time = TWave1D<double,double>(good_point_counter);
	r.freq = TWave1D<double,double>(good_point_counter);
	r.amp = TWave1D<double,double>(good_point_counter);
	r.phase = TWave1D<double,double>(good_point_counter);
	r.jitter = TWave1D<double,double>(good_point_counter);
	int j = 0;	// iterator along freqs/amps/etc...
	for(int i = 0; i<good_point_counter; i++)
	{
		// skip all the bad ones (correction level >= 1)
		if(correctionlevel() >= 1)
		{
			while(jitter_free_points[j] == 0)
				j++;
		}

		r.time[i] = freqs.t(j);
		r.amp[i] = amps[j];
		r.phase[i] = phases[j];

		if(correctionlevel() >= 2)	// correction level 2 says we should phase-cohere our frequencies
			r.freq[i] = freqs[j] - jitters_smth[j]/window_size/wave.dt()/PI;
		else
			r.freq[i] = freqs[j];
		 
		
		// report jitter as a frequency
		if(correctionlevel() >= 2)
			r.jitter[i] = (jitters[j] - jitters_smth[j])/window_size/wave.dt()/PI;
		else	
			r.jitter[i] = jitters[j]/window_size/wave.dt()/PI;
		
		j++;
		
		// from 98% to 99%
		percent_done = 98.0 + 1.0*j/freqs.dim();
	}
	
	r.time.copyTime(freqs);
	r.freq.copyTime(freqs);
	r.amp.copyTime(freqs);
	r.phase.copyTime(freqs);
	
	// claim that we're done
	percent_done = 100.0;
	return r;
}


template <class T, class Tbigger>
class FixptFourierDemodulator : public FourierDemodulator<T, Tbigger>
{
	public:
	FixptFourierDemodulator(const TWave1D<T,double>& wave_, const int window_size_) : FourierDemodulator<T,Tbigger>(wave_,window_size_) { }
	FixptFourierDemodulator(const FixptFourierDemodulator<T,Tbigger>& old) : FourierDemodulator<T,Tbigger>(old) { }
	FixptFourierDemodulator() : FourierDemodulator<T,Tbigger>() { }

	protected:
	fint_ret_t four(const TWave1D<Tbigger,double>& data, double freq);	
	
	void setones()
	{
		// trig bits allows division by 2^n to be done as the faster bitwise operations
		TRIG_BITS = sizeof(Tbigger)*8/2-2;
		TRIG_ONE = (Tbigger)1<<TRIG_BITS;
			
		FourierDemodulator<T,Tbigger>::MATH_ONE = (Tbigger) ((int)1<<(sizeof(T)*8));
	}
	
	int TRIG_BITS;	// TRIG_BITS and TRIG_ONE automatically set. size of 1.0 in fixed point math for trig
	Tbigger TRIG_ONE;

	
};

/** Compute a Fourier integral on integer input data in fixed point as fast as possible */
template <class T,class Tbigger>
fint_ret_t FixptFourierDemodulator<T,Tbigger>::four(const TWave1D<Tbigger,double>& data, double freq)
{
	double delta = 2*FourierDemodulator<T,Tbigger>::PI*freq*data.dt();
	
	Tbigger alpha = (Tbigger)((2.0*sin(delta/2.0)*sin(delta/2.0))*TRIG_ONE);
	Tbigger beta = (Tbigger)(sin(delta)*TRIG_ONE);
	
	Tbigger costh = TRIG_ONE;
	Tbigger sinth = 0;
	Tbigger temp;
	Tbigger a = 0;
	Tbigger b = 0;
	
	//	cout << "alpha = " << alpha << "(" << (2.0*sin(delta/2.0)*sin(delta/2.0)) << ")" << endl;
	//	cout << "beta = " << beta << "(" << sin(delta) << ")" << endl;
	
	// here's the most important loop
	for(int i=0; i < data.dim(); i++)
	{
		// fast trig using recursion
		//		temp = costh - (alpha*costh + beta*sinth)/TRIG_ONE;
		//		sinth = sinth - (alpha*sinth - beta*costh)/TRIG_ONE;
		//cout << "costh = {" << costh - (alpha*costh + beta*sinth)/TRIG_ONE << ", ";
		// 
		temp = costh - ((alpha*costh + beta*sinth)>>TRIG_BITS);
		sinth = sinth - ((alpha*sinth - beta*costh)>>TRIG_BITS);
		
		costh = temp;
		//cout << costh << "}" << endl;
		
		b += (data[i] * (sinth));
		a += (data[i] * (temp));
	}
	
	fint_ret_t r;
	r.a = (double)a/TRIG_ONE/data.dim()*2;
	r.b = (double)b/TRIG_ONE/data.dim()*2;
	r.mag = sqrt(r.a*r.a+r.b*r.b);
	return r;
}
	

template <class T, class Tbigger>
class FloatFourierDemodulator : public FourierDemodulator<T, Tbigger>
{
public:
	FloatFourierDemodulator(const TWave1D<T,double>& wave_, const int window_size_) : FourierDemodulator<T,Tbigger>(wave_,window_size_) { }
	FloatFourierDemodulator(const FloatFourierDemodulator<T,Tbigger>& old) : FourierDemodulator<T,Tbigger>(old) { }
	FloatFourierDemodulator() : FourierDemodulator<T,Tbigger>() { }

protected:
	fint_ret_t four(const TWave1D<Tbigger,double>& data, double freq);

	void setones()
	{
		FourierDemodulator<T,Tbigger>::MATH_ONE = 1.0;
	}
	
};

template <class T, class Tbigger>
fint_ret_t FloatFourierDemodulator<T,Tbigger>::four(const TWave1D<Tbigger,double>& data, double freq)
{
	double omega = 2*FourierDemodulator<T,Tbigger>::PI*freq*data.dt();
	
	double a=0;
	double b=0;
		
	for(int i=0; i<data.dim(); i++)
	{
		a += data[i]*cos(omega*i);
		b += data[i]*sin(omega*i);
	}
		
	fint_ret_t r;
	r.a = a/0.5/data.dim();
	r.b = b/0.5/data.dim();
	
	r.mag = sqrt(r.a*r.a+r.b*r.b);
	
	return r;
}
	


#endif
