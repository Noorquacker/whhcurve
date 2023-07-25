#ifndef _PFLT_TDO_H
#define _PFLT_TDO_H
/*
 *  pflt_tdo.h
 *  
 *
 *  Created by Andy on 8/26/08.
 *
 */

//#include "ffat.h"
#include "tnt/tnt.h"
#include "savitzkygolay.h"

//Ffat loadFfat(const char* fieldfile, const char* freqfile, const char* shotname, bool subtract_pre_shot=true, double freqoffset=0, int freqsign=1);
//Array1D<double> integratePU(Array1D<double> pu);
//void splitPulse(const Ffat& shot, Ffat& up, Ffat& down);

class PfltField
{
	// by default, everything is public. if you need a piece of data, be careful with it.
	public:
		PfltField(const PfltField& old) { (*this) = old; }
		const PfltField& operator =(const PfltField& old)
		{
			time = old.time;
			field = old.field;
			pu = old.pu;
			pearson = old.pearson;
			ch3 = old.ch3;
			ch4 = old.ch4;
		
			return (*this);
		}
		PfltField(const char* filename, const char* cols) {loadField(filename,cols);}
		
		// cols: t = time, f = Field, p = Pickup coil, x = Pearson current transformer, 3 = Ch. 3, 4 = Ch. 4, ... more?
		void loadField(const char* filename, const char* cols);
		void integratePU(double dt, int pre_trigger_points, double PU_AREA, double DIV_BY_FACTOR);
		void integratePU(double PU_AREA, double DIV_BY_FACTOR);

		void savgol(SavitzkyGolay<double>& sg);

		Array1D<double> time;
		Array1D<double> field;
		Array1D<double> pu;
		Array1D<double> pearson;
		Array1D<double> ch3;
		Array1D<double> ch4;
		
		
};

class PfltFreq
{
	public:
		PfltFreq(const PfltFreq& old) { (*this) = old; }
		const PfltFreq& operator =(const PfltFreq& old)
		{
			time = old.time;
			freq = old.freq;
			amp = old.amp;
			phase = old.phase;
			jitter = old.jitter;
			
			return (*this);
		}
		PfltFreq(const char* filename, const char* cols) {loadFreq(filename,cols);}

		// cols: t = time, f = freq, a = amp, j = jitter
		void loadFreq(const char* filename, const char* cols);
		void savgol(SavitzkyGolay<double>& sg);

		Array1D<double> time;
		Array1D<double> freq;
		Array1D<double> amp;
		Array1D<double> phase;
		Array1D<double> jitter;
	
};

#endif

