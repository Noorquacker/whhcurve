#ifndef _FILTERS_H
#define _FILTERS_H

#include "math.h"

#include "twave1d.h"
#include "fwave1d.h"

/** Filters that use FFTs. Be very careful using on
	waves that do not start and end at the same value
*/
template <class T>
TWave1D<T,T> filterLP(TWave1D<T,T> inwave, double fCutoff, T attPerOct)
{
    FReal1D fcs;
    fcs.isPeriodic(false);
    fcs = inwave;	// do FFT
    
    for(unsigned i=0; i<fcs.dim(); i++)
    {
        if(fcs.f(i) > 0)
			fcs[i] /= (1.0 + pow(fcs.f(i) / fCutoff, (attPerOct / 6.0)));

// old one doesn't work
//            fcs[i] *= exp( -1.0/fCutoff*log(2.0)/attPerOct*(fcs.f(i)-fCutoff) );
    }
    
    TReal1D yfilt = fcs;	// reverse FFT
    return yfilt;
}
/** Filters that use FFTs. Be very careful using on
	waves that do not start and end at the same value
*/
template <class T>
TWave1D<T,T> filterHP(TWave1D<T,T> inwave, double fCutoff, T attPerOct)
{
    FReal1D fcs;
    fcs.isPeriodic(false);
    fcs = inwave;	// do FFT
    
    for(unsigned i=0; i<fcs.dim(); i++)
    {
        if(fcs.f(i) <= 0)
            fcs[i] = 0.0;
		else
			fcs[i] /= (1.0 + pow(fCutoff / fcs.f(i), (attPerOct / 6.0)));

//            fcs[i] *= exp( -1.0/fCutoff*log(2.0)/attPerOct*(fCutoff-fcs.f(i)) );
    }
    
    TReal1D yfilt = fcs;	// reverse FFT
    return yfilt;
}
/** Quick and dirty slope formula. */
template <class T>
T qdSlope(const T x1, const T y1, const T x2, const T y2)
{
    return (y2-y1)/(x2-x1);
}

/** Quick and dirty (qd) derivatives. At the endpoints, takes a two-point slope. For all other points
	(a,b,c), it takes m_ab/4 + m_bc/4 + m_ac/2. Superceded by SavitkyGolay class.
*/
template <class T>
TWave1D<T,T> qdDeriv(const TWave1D<T,T> wave)
{
    TWave1D<T,T> deriv(wave.dim());
    deriv.copyTime(wave);
    
    for(int i=0; i<wave.dim(); i++)
    {
        if(i==0)
        {
            deriv[0] = qdSlope<double>(wave.t(0),wave[0],wave.t(1),wave[1]);
            continue;
        }
        if(i==(wave.dim()-1))
        {
            deriv[i] = qdSlope<double>(wave.t(i-1),wave[i-1],wave.t(i),wave[i]);
            continue;
        }
        deriv[i] = 
            qdSlope<double>(wave.t(i-1),wave[i-1],wave.t(i),wave[i])/4 +
            qdSlope<double>(wave.t(i),wave[i],wave.t(i+1),wave[i+1])/4 +
            qdSlope<double>(wave.t(i-1),wave[i-1],wave.t(i+1),wave[i+1])/2;
    }
    return deriv;
}

#endif
