#include <strstream>
#include <fstream>
#include <stdlib.h>

#include "ffat.h"
#include "filters.h"
#include "textwavedata.h"
#include "smoothsort.h"

using namespace std;

const double PI = 3.141592654;

/** load a ffat from loaded field and frequency data */
Ffat::Ffat(const char* n, PfltField pfltfi, PfltFreq pfltfr, double smoothing = 0.0)
{
	name(n);
	
	if(pfltfi.pu.dim()==0)
		throw "PU array must exist when splining field and frequency into a Ffat.";
	if(pfltfr.time.dim()==0)
		throw "Frequency data must exist when splining into a Ffat.";
	
	if(pfltfi.time.dim()!=pfltfi.pu.dim())
		pfltfi.time = Array1D<double>(pfltfi.pu.dim(),0.0);
	if(pfltfi.field.dim()!=pfltfi.pu.dim())
		pfltfi.field = Array1D<double>(pfltfi.pu.dim(),0.0);
	if(pfltfi.pearson.dim()!=pfltfi.pu.dim())
		pfltfi.pearson = Array1D<double>(pfltfi.pu.dim(),0.0);
	if(pfltfi.ch3.dim()!=pfltfi.pu.dim())
		pfltfi.ch3 = Array1D<double>(pfltfi.pu.dim(),0.0);
	if(pfltfi.ch4.dim()!=pfltfi.pu.dim())
		pfltfi.ch4 = Array1D<double>(pfltfi.pu.dim(),0.0);
	
	if(pfltfr.freq.dim()!=pfltfr.time.dim())
		pfltfr.freq = Array1D<double>(pfltfr.time.dim(),0.0);
	if(pfltfr.amp.dim()!=pfltfr.time.dim())
		pfltfr.amp = Array1D<double>(pfltfr.time.dim(),0.0);
	if(pfltfr.phase.dim()!=pfltfr.time.dim())
		pfltfr.phase = Array1D<double>(pfltfr.time.dim(),0.0);
	if(pfltfr.jitter.dim()!=pfltfr.time.dim())
		pfltfr.jitter = Array1D<double>(pfltfr.time.dim(),0.0);

	
	// before we do anything else, gather the pre-shot data
	Array1D<double> averages(9);
	Array1D<double> stddevs(9);
	
	// field first:
	double avfield = 0; double stfield = 0;
	double avpu = 0; double stpu = 0;
	double avpearson = 0; double stpearson = 0;
	double avch3 = 0; double stch3 = 0;
	double avch4 = 0; double stch4 = 0;
	int i = 0;
	
	// average. use a short window of time before the trigger
	for(int j=0; pfltfi.time[j] < (pfltfi.time[0]*0.05); j++)
	{
		if(pfltfi.time[j] > pfltfi.time[0]*0.1)
		{
			i++;
			avfield += pfltfi.field[j];
			avpu += pfltfi.pu[j];
			avpearson += pfltfi.pearson[j];
			avch3 += pfltfi.ch3[j];
			avch4 += pfltfi.ch4[j];
		}
	}
	averages[0] = avfield /= i;
	averages[5] = avpu /= i;
	averages[6] = avpearson /= i;
	averages[7] = avch3 /= i;
	averages[8] = avch4 /= i;

	// standard deviation. uses same data as above
	for(int j=0; pfltfi.time[j] < (pfltfi.time[0]*0.05); j++)
	{
		if(pfltfi.time[j] > pfltfi.time[0]*0.1)
		{
			stfield += (pfltfi.field[j] - avfield)*(pfltfi.field[j] - avfield);
			stpu += (pfltfi.pu[j] - avpu)*(pfltfi.pu[j] - avpu);
			stpearson += (pfltfi.pearson[j] - avpearson)*(pfltfi.pearson[j] - avpearson);
			stch3 += (pfltfi.ch3[j] - avch3)*(pfltfi.ch3[j] - avch3);
			stch4 += (pfltfi.ch4[j] - avch4)*(pfltfi.ch4[j] - avch4);
		}
	}
	stddevs[0] = stfield = sqrt(stfield/i);
	stddevs[5] = stpu = sqrt(stpu/i);
	stddevs[6] = stpearson = sqrt(stpearson/i);
	stddevs[7] = stch3 = sqrt(stch3/i);
	stddevs[8] = stch4 = sqrt(stch4/i);
	
	// now frequency stuff
	double avfreq = 0; double stfreq = 0;
	double avamp = 0; double stamp = 0;
	double avphase = 0; double stphase = 0;
	double avjitter = 0; double stjitter = 0;
	// average
	i = 0;
	for(int j=0; pfltfr.time[j] < (pfltfr.time[0]*0.05); j++)
	{
		if(pfltfr.time[j] > pfltfr.time[0]*0.1)
		{
			i++;
			avfreq += pfltfr.freq[j];
			avamp += pfltfr.amp[j];
			avphase += pfltfr.phase[j];
			avjitter += pfltfr.jitter[j];
		}
	}
	averages[1] = avfreq /= i;
	averages[2] = avamp /= i;
	averages[3] = avphase /= i;
	averages[4] = avjitter /= i;
	// standard deviation
	for(int j=0; pfltfr.time[j]<(pfltfr.time[0]*0.05); j++)
	{
		if(pfltfr.time[j] > pfltfr.time[0]*0.1)
		{
			stfreq += (pfltfr.freq[j] - avfreq)*(pfltfr.freq[j] - avfreq);
			stamp += (pfltfr.amp[j] - avamp)*(pfltfr.amp[j] - avamp);
			stphase += (pfltfr.phase[j] - avphase)*(pfltfr.phase[j] - avphase);
			stjitter += (pfltfr.jitter[j] - avjitter)*(pfltfr.jitter[j] - avjitter);
		}
	}
	stddevs[1] = stfreq = sqrt(stfreq/i);
	stddevs[2] = stamp = sqrt(stamp/i);
	stddevs[3] = stphase = sqrt(stphase/i);
	stddevs[4] = stjitter = sqrt(stjitter/i);

	initialValues(averages);
	initialStdDev(stddevs);

	// use whatever interval has the slowest acquisition rate
	double dt = MAX<double>(pfltfi.time[1]-pfltfi.time[0], pfltfr.time[1]-pfltfr.time[0]);
	double tstart = MAX<double>(pfltfi.time[0], pfltfr.time[0]);
	double maxt = MIN<double>(pfltfi.time[pfltfi.time.dim()-1], pfltfr.time[pfltfr.time.dim()-1]);
	
	CubicSpline<double> csfreq(pfltfr.time, (pfltfr.freq.dim()==pfltfr.time.dim() ? pfltfr.freq : Array1D<double>(pfltfr.time.dim(),0.0)), stfreq, smoothing);
	CubicSpline<double> csamp(pfltfr.time, (pfltfr.amp.dim()==pfltfr.time.dim() ? pfltfr.amp : Array1D<double>(pfltfr.time.dim(),0.0)), stamp, smoothing);
	CubicSpline<double> csphase(pfltfr.time, (pfltfr.phase.dim()==pfltfr.time.dim() ? pfltfr.phase : Array1D<double>(pfltfr.time.dim(),0.0)), stphase, smoothing);
	CubicSpline<double> csjitter(pfltfr.time, (pfltfr.jitter.dim()==pfltfr.time.dim() ? pfltfr.jitter : Array1D<double>(pfltfr.time.dim(),0.0)), stjitter, smoothing);
	CubicSpline<double> csfield(pfltfi.time, (pfltfi.field.dim()==pfltfi.time.dim() ? pfltfi.field : Array1D<double>(pfltfi.time.dim(),0.0)), stfield, smoothing);
	CubicSpline<double> cspu(pfltfi.time, (pfltfi.pu.dim()==pfltfi.time.dim() ? pfltfi.pu : Array1D<double>(pfltfi.time.dim(),0.0)), stpu, smoothing);
	CubicSpline<double> cspearson(pfltfi.time, (pfltfi.pearson.dim()==pfltfi.time.dim() ? pfltfi.pearson : Array1D<double>(pfltfi.time.dim(),0.0)), stpearson, smoothing);
	CubicSpline<double> csch3(pfltfi.time, (pfltfi.ch3.dim()==pfltfi.time.dim() ? pfltfi.ch3 : Array1D<double>(pfltfi.time.dim(),0.0)), stch3, smoothing);
	CubicSpline<double> csch4(pfltfi.time, (pfltfi.ch4.dim()==pfltfi.time.dim() ? pfltfi.ch4 : Array1D<double>(pfltfi.time.dim(),0.0)), stch4, smoothing);
	
    this->freq = csfreq.toWaveDT(dt, tstart, maxt);
    this->amp = csamp.toWaveDT(dt, tstart, maxt);
    this->field = csfield.toWaveDT(dt, tstart, maxt);
    this->phase = csphase.toWaveDT(dt, tstart, maxt);
    this->jitter = csjitter.toWaveDT(dt, tstart, maxt);
    this->pu = cspu.toWaveDT(dt, tstart, maxt);
    this->pearson = cspearson.toWaveDT(dt, tstart, maxt);
    this->ch3 = csch3.toWaveDT(dt, tstart, maxt);
    this->ch4 = csch4.toWaveDT(dt, tstart, maxt);
    this->time = (this->freq).t();

//cerr << "first time point: " << this->time[0] << endl;
	FFAT_MODE = FFAT_MODE_TIME;
}

/** copy constructor. uses references. if you want a seperate copy, use copy() */
void Ffat::ref(const Ffat& old)
{
	populate(old.field,old.freq,old.amp,old.time,old.phase,old.jitter,old.pu,old.pearson,old.ch3,old.ch4);
	isValid();
	name(old.name());
	freq.dt(old.freq.dt());
	freq.startt(old.freq.startt());
	amp.copyTime(freq);
	time.copyTime(freq);
	field.copyTime(freq);
	phase.copyTime(freq);
	jitter.copyTime(freq);
	pu.copyTime(freq);
	pearson.copyTime(freq);
	ch3.copyTime(freq);
	ch4.copyTime(freq);
	
	initialValues(old.initialValues());
	initialStdDev(old.initialStdDev());
//cerr << "in Ffat::ref(), DB = " << old.freq.dt() << endl;
}

/** return a seperate copy of this Ffat. */
Ffat Ffat::copy()
{
	Ffat a(name(),field.copy(),freq.copy(),amp.copy(),time.copy(),phase.copy(),jitter.copy(),pu.copy(),pearson.copy(),ch3.copy(),ch4.copy());

	return a;
}

/** splitPulse() splits the pulse into up and down sweeps. Does not sort the resulting Ffats. Pass the
 original Ffat in shot, pass empty ones in up and down. splitPulse() will populate the up/down parameters.
 */
void Ffat::splitPulse(Ffat& up, Ffat& down)
{
	int maxfieldindex = maxIndex(field);

	up.name(name_ + string("u"));
	up.populate(
		field.subarray(0,maxfieldindex).copy(),
		freq.subarray(0,maxfieldindex).copy(),
		amp.subarray(0,maxfieldindex).copy(),
		time.subarray(0,maxfieldindex).copy(),
		phase.subarray(0,maxfieldindex).copy(),
		jitter.subarray(0,maxfieldindex).copy(),
		pu.subarray(0,maxfieldindex).copy(),
		pearson.subarray(0,maxfieldindex).copy(),
		ch3.subarray(0,maxfieldindex).copy(),
		ch4.subarray(0,maxfieldindex).copy());

	down.name(name_ + string("d"));
	down.populate(
		field.subarray(maxfieldindex,field.dim()-1).copy(),
		freq.subarray(maxfieldindex,field.dim()-1).copy(),
		amp.subarray(maxfieldindex,field.dim()-1).copy(),
		time.subarray(maxfieldindex,field.dim()-1).copy(),
		phase.subarray(maxfieldindex,field.dim()-1).copy(),
		jitter.subarray(maxfieldindex,field.dim()-1).copy(),
		pu.subarray(maxfieldindex,field.dim()-1).copy(),
		pearson.subarray(maxfieldindex,field.dim()-1).copy(),
		ch3.subarray(maxfieldindex,field.dim()-1).copy(),
		ch4.subarray(maxfieldindex,field.dim()-1).copy());

	down.reverse();
	
	up.initialValues(initialValues());
	up.initialStdDev(initialStdDev());
	down.initialValues(initialValues());
	down.initialStdDev(initialStdDev());
}

/** spline data in constant time. smoothing is for follows CubicSpline::constrainN() */
void Ffat::tspline(double tstart, double dt, double smoothing = 0.0)
{
	// sort first, or cubic splines might fail. (built in to cubic spline)
	//	tsort();
	
	CubicSpline<double> csfield, csfreq, csamp, csphase, csjitter, cspu, cspearson, csch3, csch4;
	if(!smoothing || initial_stddev.dim()==0)
	{
		csfield.create(this->time, this->field);
		csfreq.create(this->time, this->freq);
		csamp.create(this->time, this->amp);
		csphase.create(this->time, this->phase);
		csjitter.create(this->time, this->jitter);
		cspu.create(this->time, this->pu);
		cspearson.create(this->time, this->pearson);
		csch3.create(this->time, this->ch3);
		csch4.create(this->time, this->ch4);
    }
	else
	{
		csfield.create(this->time, this->field);
		csfreq.create(this->time, this->freq, this->initial_stddev[1], smoothing);
		csamp.create(this->time, this->amp, this->initial_stddev[2], smoothing);
		csphase.create(this->time, this->phase, this->initial_stddev[3], smoothing);
		csjitter.create(this->time, this->jitter, this->initial_stddev[4], smoothing);
		cspu.create(this->time, this->pu, this->initial_stddev[5], smoothing);
		cspearson.create(this->time, this->pearson, this->initial_stddev[6], smoothing);
		csch3.create(this->time, this->ch3, this->initial_stddev[7], smoothing);
		csch4.create(this->time, this->ch4, this->initial_stddev[8], smoothing);
	}
	
    double mint = csfreq.minx();
    double maxt = csfreq.maxx();
	//cerr << "min: " << mint << ", max: " << maxt << endl;
	
	// find start and end points, aligned on db spacing from zero
	mint = MAX<double>(tstart,mint);
	int tstart_int = (int)(mint/dt);
	if(mint > 0 && mint != (tstart_int*dt))
		tstart = (tstart_int+1)*dt;
	else
		tstart = tstart_int*dt;
	
	//	int i;
	//    for(i = 0 ; i*dt < tstart || i*dt < mint; i++);
	//    tstart = i*dt;
	//cerr << "mint established, ";
    
    this->freq = csfreq.toWaveDT(dt, tstart, maxt);
    this->amp = csamp.toWaveDT(dt, tstart, maxt);
    this->field = csfield.toWaveDT(dt, tstart, maxt);
    this->phase = csphase.toWaveDT(dt, tstart, maxt);
    this->jitter = csjitter.toWaveDT(dt, tstart, maxt);
    this->pu = cspu.toWaveDT(dt, tstart, maxt);
    this->pearson = cspearson.toWaveDT(dt, tstart, maxt);
    this->ch3 = csch3.toWaveDT(dt, tstart, maxt);
    this->ch4 = csch4.toWaveDT(dt, tstart, maxt);
    this->time = (this->freq).t();
	
	FFAT_MODE = FFAT_MODE_TIME;
}

/** spline data in constant field. smoothing is for follows CubicSpline::constrainN() */
void Ffat::bspline(double bstart, double db, double smoothing = 0.0)
{
	// sort first, or cubic splines might fail. (built in to cubic spline)
	//	bsort();
	CubicSpline<double> csfreq, csamp, cstime, csphase, csjitter, cspu, cspearson, csch3, csch4;
	if(!smoothing || initial_stddev.dim()==0)
	{
		csfreq.create(this->field, this->freq);
		csamp.create(this->field, this->amp);
		cstime.create(this->field, this->time);
		csphase.create(this->field, this->phase);
		csjitter.create(this->field, this->jitter);
		cspu.create(this->field, this->pu);
		cspearson.create(this->field, this->pearson);
		csch3.create(this->field, this->ch3);
		csch4.create(this->field, this->ch4);		
	}
	else
	{
		csfreq.create(this->field, this->freq, this->initial_stddev[1], smoothing);
		csamp.create(this->field, this->amp, this->initial_stddev[2], smoothing);
		cstime.create(this->field, this->time);
		csphase.create(this->field, this->phase, this->initial_stddev[3], smoothing);
		csjitter.create(this->field, this->jitter, this->initial_stddev[4], smoothing);
		cspu.create(this->field, this->pu, this->initial_stddev[5], smoothing);
		cspearson.create(this->field, this->pearson, this->initial_stddev[6], smoothing);
		csch3.create(this->field, this->ch3, this->initial_stddev[7], smoothing);
		csch4.create(this->field, this->ch4, this->initial_stddev[8], smoothing);
		
	}
    
    double minB = csfreq.minx();
    double maxB = csfreq.maxx();
	
	cerr << "minB=" << minB << " maxB=" << maxB << "..." << endl;

	// find start and end points, aligned on db spacing from zero
	minB = MAX<double>(bstart,minB);
	int bstart_int = (int)(minB/db);
	if(minB > 0 && minB != (bstart_int*db))
		bstart = (bstart_int+1)*db;
	else
		bstart = bstart_int*db;
	
	cerr << "bstart=" << bstart << " db=" << db << "..." << endl;


	//	int i;
	//    for(i = 0 ; i*db < bstart || i*db < minB; i++);
	//    bstart = i*db;
	
    this->freq = csfreq.toWaveDT(db, bstart, maxB);
    this->amp = csamp.toWaveDT(db, bstart, maxB);
    this->time = cstime.toWaveDT(db, bstart, maxB);
    this->phase = csphase.toWaveDT(db, bstart, maxB);
    this->jitter = csjitter.toWaveDT(db, bstart, maxB);
    this->pu = cspu.toWaveDT(db, bstart, maxB);
    this->pearson = cspearson.toWaveDT(db, bstart, maxB);
    this->ch3 = csch3.toWaveDT(db, bstart, maxB);
    this->ch4 = csch4.toWaveDT(db, bstart, maxB);
    this->field = (this->freq).t();

	cerr << "done." << endl;

	FFAT_MODE = FFAT_MODE_FIELD;
}

/** spline data in constant inverse field for SdH calculation. smoothing follows CubicSpline::constrainN() */
void Ffat::invbspline(double invbmin, double invbmax, double dinvb, double smoothing = 0.0)
{
	// first invert the field
	for(int i=0; i<=this->field.dim(); i++)
	{
		this->field[i] = 1.0/this->field[i];
	}

	// sort first, or cubic splines might fail. (built in to cubic spline)
	//	bsort();
	CubicSpline<double> csfreq, csamp, cstime, csphase, csjitter, cspu, cspearson, csch3, csch4;
	
	if(!smoothing || initial_stddev.dim()==0)
	{
		csfreq.create(this->field, this->freq);
		csamp.create(this->field, this->amp);
		cstime.create(this->field, this->time);
		csphase.create(this->field, this->phase);
		csjitter.create(this->field, this->jitter);
		cspu.create(this->field, this->pu);
		cspearson.create(this->field, this->pearson);
		csch3.create(this->field, this->ch3);
		csch4.create(this->field, this->ch4);		
	}
	else
	{
		csfreq.create(this->field, this->freq, this->initial_stddev[1], smoothing);
		csamp.create(this->field, this->amp, this->initial_stddev[2], smoothing);
		cstime.create(this->field, this->time);
		csphase.create(this->field, this->phase, this->initial_stddev[3], smoothing);
		csjitter.create(this->field, this->jitter, this->initial_stddev[4], smoothing);
		cspu.create(this->field, this->pu, this->initial_stddev[5], smoothing);
		cspearson.create(this->field, this->pearson, this->initial_stddev[6], smoothing);
		csch3.create(this->field, this->ch3, this->initial_stddev[7], smoothing);
		csch4.create(this->field, this->ch4, this->initial_stddev[8], smoothing);
		
	}
    
    double minB = csfreq.minx();
    double maxB = MIN<double>(csfreq.maxx(),invbmax);	// you want to pass invbmax, so you don't look at 0.1 T data
	double bstart;
	
	// find start and end points, aligned on db spacing from zero
	minB = MAX<double>(invbmin,minB);
	int bstart_int = (int)(minB/dinvb);
	if(minB > 0 && minB != (bstart_int*dinvb))
		bstart = (bstart_int+1)*dinvb;
	else
		bstart = bstart_int*dinvb;
	
	//	int i;
	//    for(i = 0 ; i*db < bstart || i*db < minB; i++);
	//    bstart = i*db;
	
    this->freq = csfreq.toWaveDT(dinvb, bstart, maxB);
    this->amp = csamp.toWaveDT(dinvb, bstart, maxB);
    this->time = cstime.toWaveDT(dinvb, bstart, maxB);
    this->phase = csphase.toWaveDT(dinvb, bstart, maxB);
    this->jitter = csjitter.toWaveDT(dinvb, bstart, maxB);
    this->pu = cspu.toWaveDT(dinvb, bstart, maxB);
    this->pearson = cspearson.toWaveDT(dinvb, bstart, maxB);
    this->ch3 = csch3.toWaveDT(dinvb, bstart, maxB);
    this->ch4 = csch4.toWaveDT(dinvb, bstart, maxB);
    this->field = (this->freq).t();
	
	FFAT_MODE = FFAT_MODE_INVFIELD;
}

void Ffat::fftmagsq()
{
	if(FFAT_MODE == FFAT_MODE_FFTMAGSQ)
	{
		throw "Already in FFT-MAGNITUDE-SQURE MODE in fftmagsq()";
	}

//	cerr << "first field point in subtraction: " << (this->field)[0] << endl;

	FReal1D tempf;
	tempf.isPeriodic(1);

	if(FFAT_MODE == FFAT_MODE_TIME)
	{
		// do stuff in the time domain
		tempf.startt(this->time[0]);
		tempf.dt(((this->time[(this->time.dim())-1])-(this->time[0]))/(this->time.dim()-1));

		tempf = this->field; this->field = tempf.magsqt();
		tempf = this->time; this->time = tempf.magsqtf();
		
	}
	if(FFAT_MODE == FFAT_MODE_FIELD || FFAT_MODE == FFAT_MODE_INVFIELD )
	{
		// do stuff in the field or 1/field domain
		tempf = this->time;
		tempf.startt(this->field[0]);
		tempf.dt(((this->field[(this->field.dim())-1])-(this->field[0]))/(this->field.dim()-1));
		this->time = tempf.magsqt();

		tempf = this->field;
		tempf.startt(this->field[0]);
		tempf.dt(((this->field[(this->field.dim())-1])-(this->field[0]))/(this->field.dim()-1));
		this->field = tempf.magsqtf();
//cerr << "Domain parameters of tempf: dt:" << tempf.dt() << ", startt: " << tempf.startt() << endl;

	}
	
	// now do everything else that doesn't depend on domain

	tempf = this->freq; this->freq = tempf.magsqt();
	tempf = this->amp; this->amp = tempf.magsqt();
	tempf = this->phase; this->phase = tempf.magsqt();
	tempf = this->jitter; this->jitter = tempf.magsqt();
	tempf = this->pu; this->pu = tempf.magsqt();
	tempf = this->pearson; this->pearson = tempf.magsqt();
	tempf = this->ch3; this->ch3 = tempf.magsqt();
	tempf = this->ch4; this->ch4 = tempf.magsqt();

	FFAT_MODE = FFAT_MODE_FFTMAGSQ;
}

/** ensurebspace()
	To avoid misbehaving cubic splines when splining in field on the downsweep, throw out (perhaps
	later we should change this to average) data points that are too close together in field.
*/
void Ffat::ensurebspace(double space)
{
	bsort();

	Array1D<double> fi(field.dim()), fr(field.dim()), am(field.dim()), ti(field.dim()),
		 ph(field.dim()), ji(field.dim()), pu_(field.dim()), pe(field.dim()), c4(field.dim()), c3(field.dim());
	
	fi[0] = field[0]; fr[0] = freq[0]; am[0] = amp[0]; ti[0] = time[0];
	ph[0] = phase[0]; ji[0] = jitter[0]; pu_[0] = pu[0]; pe[0] = pearson[0]; c3[0] = ch3[0]; c4[0] = ch4[0];
	
	int j = 0;
	for(int i=1; i<field.dim(); i++)
	{
		if(field[i] > (fi[j] + space))
		{
			// point is good. add it to the list
			j++;
			fi[j] = field[i];
			fr[j] = freq[i];
			am[j] = amp[i];
			ti[j] = time[i];
			ph[j] = phase[i];
			ji[j] = jitter[i];
			pu_[j] = pu[i];
			pe[j] = pearson[i];
			c3[j] = ch3[i];
			c4[j] = ch4[i];			
		}
		// do nothing with this point
	}
	populate(fi.subarray(0,j),fr.subarray(0,j),am.subarray(0,j),ti.subarray(0,j),ph.subarray(0,j),ji.subarray(0,j),pu_.subarray(0,j),pe.subarray(0,j),c3.subarray(0,j),c4.subarray(0,j));
}
/** ensuretspace()
        To avoid misbehaving cubic splines when splining in field on the downsweep, throw out (perhaps
        later we should change this to average) data points that are too close together in field.
*/
void Ffat::ensuretspace(double space)
{
  tsort();

	Array1D<double> fi(field.dim()), fr(field.dim()), am(field.dim()), ti(field.dim()),
	ph(field.dim()), ji(field.dim()), pu_(field.dim()), pe(field.dim()), c4(field.dim()), c3(field.dim());
	
	fi[0] = field[0]; fr[0] = freq[0]; am[0] = amp[0]; ti[0] = time[0];
	ph[0] = phase[0]; ji[0] = jitter[0]; pu_[0] = pu[0]; pe[0] = pearson[0]; c3[0] = ch3[0]; c4[0] = ch4[0];
	
	int j = 0;
	for(int i=1; i<field.dim(); i++)
    {
      if(time[i] > (ti[j] + space))
	{
	  // point is good. add it to the list
	  j++;
		fi[j] = field[i];
		fr[j] = freq[i];
		am[j] = amp[i];
		ti[j] = time[i];
		ph[j] = phase[i];
		ji[j] = jitter[i];
		pu_[j] = pu[i];
		pe[j] = pearson[i];
		c3[j] = ch3[i];
		c4[j] = ch4[i];			
	}
      // do nothing with this point
    }
	populate(fi.subarray(0,j),fr.subarray(0,j),am.subarray(0,j),ti.subarray(0,j),ph.subarray(0,j),ji.subarray(0,j),pu_.subarray(0,j),pe.subarray(0,j),c3.subarray(0,j),c4.subarray(0,j));
}

/** reverse()
	Significant performance gains can be found by starting arrays in the proper order before ensurebspace().
	Otherwise, sort() must go through them and perform a very expensive reversal operation. Do it here in
	linear time. Update: sorting functions take care of direction if you set it properly. Reverse function
	included in case you need it for whatever reasonn.
*/
void Ffat::reverse()
{
	int i = 0;
	int j = this->freq.dim() - 1;
	while(i < j)
	{
		// swap the ith and jth entries
		double temp;
		temp = this->field[i]; this->field[i] = this->field[j]; this->field[j] = temp;
		temp = this->freq[i]; this->freq[i] = this->freq[j]; this->freq[j] = temp;
		temp = this->amp[i]; this->amp[i] = this->amp[j]; this->amp[j] = temp;
		temp = this->time[i]; this->time[i] = this->time[j]; this->time[j] = temp;
		temp = this->phase[i]; this->phase[i] = this->phase[j]; this->phase[j] = temp;
		temp = this->jitter[i]; this->jitter[i] = this->jitter[j]; this->jitter[j] = temp;
		temp = this->pu[i]; this->pu[i] = this->pu[j]; this->pu[j] = temp;
		temp = this->pearson[i]; this->pearson[i] = this->pearson[j]; this->pearson[j] = temp;
		temp = this->ch3[i]; this->ch3[i] = this->ch3[j]; this->ch3[j] = temp;
		temp = this->ch4[i]; this->ch4[i] = this->ch4[j]; this->ch4[j] = temp;
		
		// increment i, decrement j
		i++; j--;
	}
	// now make sure the timing information is reasonable
	double dt, startt;
	dt = this->freq.dt();
	startt = this->freq.startt();
	
	startt = startt + dt*this->freq.dim();
	dt = -dt;
	this->field.dt(dt); this->field.startt(startt);
	this->freq.dt(dt); this->freq.startt(startt);
	this->amp.dt(dt); this->amp.startt(startt);
	this->time.dt(dt); this->time.startt(startt);
	this->phase.dt(dt); this->phase.startt(startt);
	this->jitter.dt(dt); this->jitter.startt(startt);
	this->pu.dt(dt); this->pu.startt(startt);
	this->pearson.dt(dt); this->pearson.startt(startt);
	this->ch3.dt(dt); this->ch3.startt(startt);
	this->ch4.dt(dt); this->ch4.startt(startt);
}

void Ffat::subtractPreShot()
{
	// it's just a subtraction for frequency and phase, but a normalization amplitude
	// no transformation is performed on the jitter.
	for(int i=0; i<freq.dim(); i++)
	{
		freq[i] -= initial_values[1];
		amp[i] = amp[i] / initial_values[2];
		phase[i] -= initial_values[3]; if(phase[i]< -PI) phase[i] += PI; if(phase[i] > PI) phase[i] -= PI;
//		jitter[i] -= initial_values[4];
		pu[i] -= initial_values[5];
		pearson[i] -= initial_values[6];
		ch3[i] -= initial_values[7];
		ch4[i] -= initial_values[8];
	}
}


/** do a savitzky-golay filtering on !(time and field) */
void Ffat::savgol(SavitzkyGolay<double>& sg)
{
	this->freq = sg.filter(this->freq);
	this->amp = sg.filter(this->amp);
	this->phase = sg.filter(this->phase);
	this->jitter = sg.filter(this->jitter);
	this->pu = sg.filter(this->pu);
	this->pearson = sg.filter(this->pearson);
	this->ch3 = sg.filter(this->ch3);
	this->ch4 = sg.filter(this->ch4);
}
/** do a savitzky-golay filtering on either time or field, depending on mode. */
void Ffat::savgolx(SavitzkyGolay<double>& sg)
{
	if(FFAT_MODE == FFAT_MODE_TIME)
		this->field = sg.filter(this->field);
	if(FFAT_MODE == FFAT_MODE_FIELD)
		this->time = sg.filter(this->time);
}

/**	When splining, we must make sure the arrays are sorted. There's a performance penalty here, but let's
	do it ourselves to make sure everything is ready to go. This is particularly important (and costly)
	when playing with downsweeps. 
*/
void Ffat::tsort()
{
	if(time[0] < time[time.dim()-1])
	{
		SmoothSort<double> s(time);
		s.reorder(field);
		s.reorder(freq);
		s.reorder(amp);
		s.reorder(phase);
		s.reorder(jitter);
		s.reorder(pu);
		s.reorder(pearson);
		s.reorder(ch3);
		s.reorder(ch4);
	}
	else
	{
		SmoothSort<double,true> s(time);
		s.reverse();
		s.reorder(field);
		s.reorder(freq);
		s.reorder(amp);
		s.reorder(phase);
		s.reorder(jitter);
		s.reorder(pu);
		s.reorder(pearson);
		s.reorder(ch3);
		s.reorder(ch4);
	}
}
void Ffat::rtsort()
{
	if(time[0] < time[time.dim()-1])
	{
		SmoothSort<double> s(time);
		s.reverse();
		s.reorder(field);
		s.reorder(freq);
		s.reorder(amp);
		s.reorder(phase);
		s.reorder(jitter);
		s.reorder(pu);
		s.reorder(pearson);
		s.reorder(ch3);
		s.reorder(ch4);
	}
	else
	{
		SmoothSort<double,true> s(time);
		s.reorder(field);
		s.reorder(freq);
		s.reorder(amp);
		s.reorder(phase);
		s.reorder(jitter);
		s.reorder(pu);
		s.reorder(pearson);
		s.reorder(ch3);
		s.reorder(ch4);
	}
}
void Ffat::bsort()
{
	if(field[0] < field[field.dim()-1])
	{
		SmoothSort<double> s(field);
		s.reorder(time);
		s.reorder(freq);
		s.reorder(amp);
		s.reorder(phase);
		s.reorder(jitter);
		s.reorder(pu);
		s.reorder(pearson);
		s.reorder(ch3);
		s.reorder(ch4);
	}
	else
	{
		SmoothSort<double,true> s(field);
		s.reverse();
		s.reorder(time);
		s.reorder(freq);
		s.reorder(amp);
		s.reorder(phase);
		s.reorder(jitter);
		s.reorder(pu);
		s.reorder(pearson);
		s.reorder(ch3);
		s.reorder(ch4);
	}
}
void Ffat::rbsort()
{
	if(field[0] < field[field.dim()-1])
	{
		SmoothSort<double> s(field);
		s.reverse();
		s.reorder(time);
		s.reorder(freq);
		s.reorder(amp);
		s.reorder(phase);
		s.reorder(jitter);
		s.reorder(pu);
		s.reorder(pearson);
		s.reorder(ch3);
		s.reorder(ch4);
	}
	else
	{
		SmoothSort<double,true> s(field);
		s.reorder(time);
		s.reorder(freq);
		s.reorder(amp);
		s.reorder(phase);
		s.reorder(jitter);
		s.reorder(pu);
		s.reorder(pearson);
		s.reorder(ch3);
		s.reorder(ch4);
	}
}

/** output inserter */
ostream& operator << (ostream& out, const Ffat& ffat)
{
    out.flags(out.flags() | ios::showpoint | ios::scientific);
    out.precision(7);

    // spit out some column headers
	out << "#s" << ffat.name() << "time\t" << "s" << ffat.name() << "field\t" << "s" << ffat.name() << "freq\t" << "s" << ffat.name() << "amp\t" << "s" << ffat.name() << "phase\t" << "s" << ffat.name() << "jitter\t" << "s" << ffat.name() << "pickup\t" << "s" << ffat.name() << "pearson\t" << "s" << ffat.name() << "ch3\t" << "s" << ffat.name() << "ch4\t" << endl;

    for(int i=0; i<ffat.time.dim(); i++)
    {
        out << ffat.time[i] << "\t" << ffat.field[i] << "\t" << ffat.freq[i] << "\t" << ffat.amp[i] << "\t" << ffat.phase[i] << "\t" << ffat.jitter[i] << "\t" << ffat.pu[i] << "\t" << ffat.pearson[i] << "\t" << ffat.ch3[i] << "\t" << ffat.ch4[i] << endl;
    }
	
	return out;
}

/** input extractor */
istream& operator >> (ifstream& in, Ffat& ffat)
{
	char pound;
	in >> pound;
	if(pound != '#')
	{
		cerr << "Error reading Ffat: first character must be a pound sign (#)." << endl;
		string s;
		getline(in,s);
		cerr << "First line: " << s << endl;
		return in;
	}
	
	string s;
	getline(in,s);
	ostrstream os;
	os << atoi((s.c_str()+2));
	ffat.name(os.str());
	
	in.seekg(0);

    TextWaveData<double> fd(in);
    fd.read(10); // read 10 columns

	ffat.populate(fd[1].copy(),fd[2].copy(),fd[3].copy(),fd[0].copy(),fd[4].copy(),fd[5].copy(),fd[6].copy(),fd[7].copy(),fd[8].copy(),fd[9].copy());

/*cerr << "Loaded Ffat." << endl << "Field: " << ffat.field[0] << endl << "Freq: " << ffat.freq[0] << endl
	<< "Amp: " << ffat.amp[0] << endl << "Time: " << ffat.time[0] << endl;*/

	return in;
}
