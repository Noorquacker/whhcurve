#ifndef _FFAT_H
#define _FFAT_H

#include <string>
#include <iostream>
#include <fstream>

#include "twave1d.h"
#include "spline.h"
#include "savitzkygolay.h"
#include "pflt_tdo.h"

using namespace std;

#define FFAT_MODE_NONE	0
#define FFAT_MODE_TIME	1
#define FFAT_MODE_FIELD	2
#define FFAT_MODE_INVFIELD	3
#define FFAT_MODE_FFTMAGSQ	4

#define FFAT_DIRECTION_NONE	0
#define FFAT_DIRECTION_UP	1
#define FFAT_DIRECTION_DOWN	2

class Ffat
{
	public:
	// constructors
	Ffat() { FFAT_MODE = FFAT_MODE_NONE; }
	Ffat(const char* n) { name(n); FFAT_MODE = FFAT_MODE_NONE; }
	Ffat(const char* n,
		 const Array1D<double>& fi,
		 const Array1D<double>& fr,
		 const Array1D<double>& am,
		 const Array1D<double>& ti,
		 const Array1D<double>& ph,
		 const Array1D<double>& ji,
		 const Array1D<double>& pu_,
		 const Array1D<double>& pe,
		 const Array1D<double>& c3,
		 const Array1D<double>& c4)
	{ name(n); populate(fi,fr,am,ti,ph,ji,pu_,pe,c3,c4); }
	Ffat(const Ffat& old) { ref(old); }
	Ffat(const char* n, PfltField pfltfi, PfltFreq pfltfr, double smoothing);
	
	const Ffat& operator = (const Ffat& old) { ref(old); return (*this); }
	
	// methods
	Ffat copy();
	void ref(const Ffat& old);

	void splitPulse(Ffat& up, Ffat& down);
	void tspline(double tstart, double dt, double smoothing);
	void bspline(double bstart, double db, double smoothing);
	void invbspline(double invbmin, double invbmax, double dinvb, double smoothing);
	void fftmagsq();
	void tsort();
	void rtsort();
	void bsort();
	void rbsort();
	void ensurebspace(double space);
	void ensuretspace(double space);
	void reverse();
	void subtractPreShot();

	void savgol(SavitzkyGolay<double>& sg);
	void savgolx(SavitzkyGolay<double>& sg);

	void populate(const Array1D<double>& fi,
				  const Array1D<double>& fr,
				  const Array1D<double>& am,
				  const Array1D<double>& ti,
				  const Array1D<double>& ph,
				  const Array1D<double>& ji,
				  const Array1D<double>& pu_,
				  const Array1D<double>& pe,
				  const Array1D<double>& c3,
				  const Array1D<double>& c4)
		{
			field.ref(fi); freq.ref(fr); amp.ref(am); time.ref(ti); phase.ref(ph); jitter.ref(ji); pu.ref(pu_); pearson.ref(pe); ch3.ref(c3); ch4.ref(c4);
			isValid(); FFAT_MODE = FFAT_MODE_NONE;  FFAT_DIRECTION = FFAT_DIRECTION_NONE; }

	const char* name() const { return name_.c_str(); }
	void name(const char* snew) { name_.assign(snew); }
	void name(const string& snew) { name_.assign(snew); }
	
	Array1D<double> initialValues() const { return initial_values; }
	Array1D<double> initialValues(Array1D<double> iv_) { return initial_values = iv_; }
	
	Array1D<double> initialStdDev() const { return initial_stddev; }
	Array1D<double> initialStdDev(Array1D<double> isd_) { return initial_stddev = isd_; }
	

	protected:
	void isValid()
		{	if(field.dim() != freq.dim() || freq.dim() != amp.dim() || amp.dim() != time.dim()
			   || time.dim() != jitter.dim() || jitter.dim() != phase.dim() || jitter.dim() != pearson.dim() || pearson.dim() != pu.dim() || pearson.dim() != ch3.dim() || ch3.dim() != ch4.dim() )
				throw "Waves in Ffat must be same size.";	}
	int maxIndex(const Array1D<double> &arr)
	{
		double maxy = arr[0];
		int maxx = 0;
		
		for(int i=0; i<arr.dim(); i++)
		{
			if(arr[i] > maxy)
			{
				maxy = arr[i];
				maxx = i;        
			}
		}
		return maxx;
	}

	// data members
	char FFAT_MODE;
	char FFAT_DIRECTION;
	
	Array1D<double> initial_values;
	Array1D<double> initial_stddev;
	
	// names
	string name_;

	public:
	// data arrays are public. be careful with them
	TReal1D field;
	TReal1D freq;
	TReal1D amp;
	TReal1D time;
	TReal1D phase;
	TReal1D jitter;
	TReal1D pu;
	TReal1D pearson;
	TReal1D ch3;
	TReal1D ch4;
	
};

ostream& operator << (ostream& out, const Ffat& ffat);
istream& operator >> (ifstream& in, Ffat& ffat);
#endif	// _FFAT_H
