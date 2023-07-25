/*
 *  pflt_tdo.cpp
 *  
 *
 *  Created by Andy on 8/26/08.
 */

#include <cstring>

#include "pflt_tdo.h"
#include "textwavedata.h"


int maxIndex(const Array1D<double> &arr);

/** loadFfat() gets field and frequency data from text files and returns them in a Ffat for our use.
	NOTE ON OPTIONAL PARAMETERS:
		subtract_pre_shot will adjust the frequency so pre-shot data is zero.
		offset is added on to the frequency.
		freqsign is an integer multiplier for frequency (before adding offset). use when adding on a mixing frequency
		 to determine raw frequency.
	For reference on these parameters, examine the line marked FREQUENCY ADJUSTMENT.
*/
/*Ffat loadFfat(const char* fieldfile, const char* freqfile, const char* shotname, bool subtract_pre_shot, double freqoffset, int freqsign)
{
	// time-domain cubic spline constants
	double CS_TSTART = 0.4e-3;
	double CS_TEND = 790e-3;	// no longer used? now end is calculated from shortest input wave.
	double CS_DT = 10e-6;

cerr << "Loading field " << fieldfile << endl;

	// istreams to read data
    ifstream fin(fieldfile);
cerr << "and freq " << freqfile << endl;
    ifstream tin(freqfile);

	// load field first
    TextWaveData<double> fd(fin);
cerr << "fd.read(1);" << endl;
    fd.read(1); // read 1 column
cerr << "close" << endl;
    fin.close();
cerr << "copy" << endl;
    Array1D<double> field(fd[0].copy());
cerr << "integratePU()" << endl;
cerr << fd[0][0] << endl; cerr.flush();
    Array1D<double> fieldtime = integratePU(field);
cerr << "cs" << endl;
    CubicSpline<double> fieldcs(fieldtime,field);

cerr << "Field done... ";

    // now do the frequency/amp
    TextWaveData<double> td(tin);
    td.read(3); // read 3 columns
    tin.close();

cerr << "Freq loaded... ";

    // load the columns into arrays
    Array1D<double> freq(td[0].copy());
    Array1D<double> famp(td[1].copy());
    Array1D<double> ftime(td[2].copy());

//cerr << endl << "freq: " << freq.dim() << endl;
	// optional frequency adjustments
	if(subtract_pre_shot || freqoffset != 0.0 || freqsign != 1)
	{
		double initialfreq = 0.0;
		if(subtract_pre_shot)
		{
			int i;
			for(i=0; ftime[i] < -7e-4; i++)
			{
				initialfreq += freq[i];
			}
			initialfreq = initialfreq / i;
		}
		for(int i=0; i < freq.dim(); i++)
		{
			freq[i] = (freq[i] - initialfreq) * freqsign + freqoffset;	// FREQUENCY ADJUSTMENT
		}
	}
	
    CubicSpline<double> freqcs(ftime,freq);
    CubicSpline<double> ampcs(ftime,famp);
//cerr << "d2/d2x points found...";

//cerr << "freq: " << freqcs.maxx() << ", field: " << fieldcs.maxx() << endl;
	double maxtime = (freqcs.maxx() < fieldcs.maxx() ? freqcs.maxx() : fieldcs.maxx());
//cerr << "maxtime: " << maxtime << endl;
	TReal1D frequency = freqcs.toWaveDT(CS_DT, CS_TSTART, maxtime);
	
//cerr << "splined... ";
	
	Ffat ffat(
		fieldcs.toWaveDT(CS_DT, CS_TSTART, maxtime),
		frequency,
		ampcs.toWaveDT(CS_DT, CS_TSTART, maxtime),
		frequency.t(),
		shotname);

//cerr << "Ffat created." << endl;
	return ffat;
}

/**	integratePU() -- Integrates pickup coil data into field data.
	memory and code-saving trickery here: this integrates the pickup
	signal in place, replacing it with the field. probably best to
	grab it from the file and call it field, remembering that you have
	to integrate it to actually GET the field. returns time, which
	you need. */
/*Array1D<double> integratePU(Array1D<double> pu)
{


cerr << "pu: "; cerr.flush(); cerr << pu[0] << endl;
	// field acquisition constants
	int PRE_SHOT_POINTS = 1000;
	double NIDAQ_DT = 5.0e-6;
	double PU_AREA = 1.0573e-2;
	double DIVIDE_BY_10_FACTOR = 11.04;	// 11.04 w/ box, 1.0 without box

    Array1D<double> time(pu.dim());
	cerr << "foo" << endl;

    // first step is to average some pre-shot data
    double v_preshot = 0;
    for(int i=0; i<PRE_SHOT_POINTS*95/100; i++)
    {
        v_preshot += pu[i];
    }
    v_preshot /= (PRE_SHOT_POINTS*95/100);
	
    // and integrate, correcting for pre-shot offset
    double curr_field = 0;
    for(int i=0; i<pu.dim(); i++)
    {
        curr_field += (pu[i]-v_preshot)/PU_AREA*DIVIDE_BY_10_FACTOR*NIDAQ_DT;
        pu[i] = curr_field;
        time[i] = (i-PRE_SHOT_POINTS)*NIDAQ_DT;
    }    

    return time;
}

/** Returns the index of the maximum element in an Array1D<T>. Used by splitPulse() */
/*int maxIndex(const Array1D<double> &arr)
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

/** splitPulse() splits the pulse into up and down sweeps. Does not sort the resulting Ffats. Pass the
	original Ffat in shot, pass empty ones in up and down. splitPulse() will populate the up/down parameters.
*/
/*void splitPulse(const Ffat& shot, Ffat& up, Ffat& down)
{
    int maxfieldindex = maxIndex(shot.field());

	Array1D<double> currfield(shot.field());
	Array1D<double> currfreq(shot.freq());
	Array1D<double> curramp(shot.amp());
	Array1D<double> currtime(shot.time());

	Array1D<double> timeup(maxfieldindex+1);
    Array1D<double> frequp(maxfieldindex+1);
    Array1D<double> ampup(maxfieldindex+1);
    Array1D<double> fieldup(maxfieldindex+1);
	Array1D<double> timedown(currtime.dim()-maxfieldindex);
    Array1D<double> freqdown(currfreq.dim()-maxfieldindex);
    Array1D<double> ampdown(curramp.dim()-maxfieldindex);
    Array1D<double> fielddown(currfield.dim()-maxfieldindex);
    
    // now copy the data into all six new arrays
    // upsweep goes up
    for(int i=0; i <= maxfieldindex; i++)
    {
		timeup[i] = currtime[i];
        frequp[i] = currfreq[i];    
        ampup[i] = curramp[i];    
        fieldup[i] = currfield[i];
    }
    
	// --> skipping this step, as we do it each time we spline any Ffat
	// sort by field
    // piksr4(fieldup,frequp,ampup,timeup);
    
    // downsweep gets copied in reverse
    for(int i=0; i < currfreq.dim()-maxfieldindex; i++)
    {
		timedown[i] = currtime[currfreq.dim()-i-1];
        freqdown[i] = currfreq[currfreq.dim()-i-1];    
        ampdown[i] = curramp[currfreq.dim()-i-1];    
        fielddown[i] = currfield[currfreq.dim()-i-1];    
    }
	// --> skipping this step, as we do it each time we spline any Ffat
    // piksr4(fielddown,freqdown,ampdown,timedown);
	
	string name(shot.basename());
	
	up.populate(fieldup,frequp,ampup,timeup);
	up.basename(name + "u"); up.makenames();
	
	down.populate(fielddown,freqdown,ampdown,timedown);
	down.basename(name + "d"); down.makenames();
}
*/
void PfltField::loadField(const char* filename, const char* cols)
{
	ifstream fin(filename);
	
	TextWaveData<double> wd(fin);
	wd.read(strlen(cols));
	
	for(int i=0; i<(int)strlen(cols); i++)
	{
		switch(tolower(cols[i]))
		{
			case 't':	// time
				time = wd[i].copy();
				break;
			case 'f':
				field = wd[i].copy();
				break;
			case 'p':
			case '1':
				pu = wd[i].copy();
				break;
			case 'x':
			case '2':
				pearson = wd[i].copy();
				break;
			case '3':
				ch3 = wd[i].copy();
				break;
			case '4':
				ch4 = wd[i].copy();
				break;
		}
	}
}

void PfltField::integratePU(double dt, int pre_trigger_points, double PU_AREA, double DIV_BY_FACTOR = 1.0)
{	
	if(pu.dim() == 0)
		throw "No pickup data.";
		
//cerr << "pu size " << pu.dim() << endl;
		
	time = Array1D<double>(pu.dim());
	field = Array1D<double>(pu.dim());
	
    // first step is to average some pre-shot data
    double v_preshot = 0;
    for(int i=0; i<pre_trigger_points*19/20; i++)
    {
        v_preshot += pu[i];
    }
    v_preshot /= (pre_trigger_points*19/20);
	
    // and integrate, correcting for pre-shot offset
    double curr_field = 0;
    for(int i=0; i<pu.dim(); i++)
    {
        curr_field += (pu[i]-v_preshot)/PU_AREA*DIV_BY_FACTOR*dt;
        field[i] = curr_field;
        time[i] = (i-pre_trigger_points)*dt;
    }
}

void PfltField::integratePU(double PU_AREA, double DIV_BY_FACTOR = 1.0)
{	
	if(pu.dim() == 0)
		throw "No pickup data.";
	if(time.dim() == 0)
		throw "integratePU called in this form needs time data. Use the long form if you do not have it.";
	if(time[0] >= 0)
		throw "integratePU needs pre-trigger points denoted by time[0] < 0.0.";
		
	int pre_trigger_points;
	for(pre_trigger_points=0; pre_trigger_points < time.dim() && time[pre_trigger_points] < 0.0; pre_trigger_points--);
	
    // first step is to average some pre-shot data
    double v_preshot = 0;
    for(int i=0; i<pre_trigger_points*95/100; i++)
    {
        v_preshot += pu[i];
    }
    v_preshot /= (pre_trigger_points*95/100);
	
    // and integrate, correcting for pre-shot offset
    double curr_field = (pu[0]-v_preshot)/PU_AREA*DIV_BY_FACTOR*0.5*(time[1]-time[0]);
	field[0] = curr_field;
    for(int i=1; i<pu.dim()-1; i++)
    {
        curr_field += (pu[i]-v_preshot)/PU_AREA*DIV_BY_FACTOR*(0.5*(time[i]-time[i-1])+0.5*(time[i+1]-time[i]));
        field[i] = curr_field;
    }
	curr_field += (pu[pu.dim()-1]-v_preshot)/PU_AREA*DIV_BY_FACTOR*(0.5*(time[pu.dim()-1]-time[pu.dim()-2]));
	field[pu.dim()-1] = curr_field;
}

void PfltField::savgol(SavitzkyGolay<double>& sg)
{
	if(field.dim()) this->field = sg.filter(this->field);
	if(pu.dim()) this->pu = sg.filter(this->pu);
	if(pearson.dim()) this->pearson = sg.filter(this->pearson);
	if(ch3.dim()) this->ch3 = sg.filter(this->ch3);
	if(ch4.dim()) this->ch4 = sg.filter(this->ch4);
}



void PfltFreq::loadFreq(const char* filename, const char* cols)
{
	ifstream fin(filename);
	
	TextWaveData<double> wd(fin);
	wd.read(strlen(cols));
	
	for(int i=0; i<(int)strlen(cols); i++)
	{
		switch(tolower(cols[i]))
		{
			case 't':	// time
				time = wd[i].copy();
				break;
			case 'f':	// freq
				freq = wd[i].copy();
				break;
			case 'a':	// amplitude
				amp = wd[i].copy();
				break;
			case 'p':
				phase = wd[i].copy();
				break;
			case 'j':	// jitter
				jitter = wd[i].copy();
				break;
		}
	}
}

void PfltFreq::savgol(SavitzkyGolay<double>& sg)
{
	if(freq.dim()) this->freq = sg.filter(this->freq);
	if(amp.dim()) this->amp = sg.filter(this->amp);
	if(jitter.dim()) this->jitter = sg.filter(this->jitter);
	// phase filtering doesn't really make sense, but do it anyway by default, assuming the user
	// knows what they are doing. maybe they're using using the demodulator as a fancy lock-in
	if(phase.dim()) this->phase = sg.filter(this->phase);
}	
		
