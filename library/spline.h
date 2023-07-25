#ifndef _SPLINE_H
#define _SPLINE_H

#include <iostream>
using namespace std;

#include "tnt/tnt.h"

using namespace TNT;

#include "tntnr.h"
#include "twave1d.h"
#include "smoothsort.h"

extern const char* SPLINE_ARRAYS_MUST_BE_SAME_SIZE;
extern const char* BAD_XA_TO_SPLINT;
extern const char* SPLINT_OUT_OF_BOUNDS_LOW;
extern const char* SPLINT_OUT_OF_BOUNDS_HIGH;
extern const char* NO_CONVERGENCE_IN_ADAPTIVE_SMOOTH_SPLINE;

/** <h1>Cubic Spline Interpolation</h1>
<p>Does cubic spline interpolation on arrays of the same length. If function
is multi-valued, only the first value encountered is used.</p>

<pre>try
    {
    string sShotNumber = argv[1];
    string sFilename = string("DataAnalysis/20071009_probe_vibration_analysis/C3_") + sShotNumber + string(".txt");
    
    TextWaveData<double> wd;
    wd.openr(sFilename.c_str());
    wd.read();

    CubicSpline<double> csf(wd[2],wd[0]);
    
    TReal1D tdofreq = csf.toWave(65536);
    }//try
    catch(const char* err)
    {
        cerr << err << endl;
    }
    
    // tdofreq is now a TReal1D with the frequency given at regular intervals
    // that may be plotted against a similarly described magnetic field
}</pre>
*/
template <class T>
class CubicSpline
{
    public:
		CubicSpline(const Array1D<T> &x_, const Array1D<T> &y_) { create(); create(x_,y_); }
		CubicSpline(const Array1D<T> &x_, const Array1D<T> &y_, const Array1D<T> &sigma_) { create(); create(x_,y_,sigma_); }
		CubicSpline(const Array1D<T> &x_, const Array1D<T> &y_, const Array1D<T> &sigma_, double smoothing) { create(); create(x_,y_,sigma_,smoothing); }
		CubicSpline(const Array1D<T> &x_, const Array1D<T> &y_, double sigma_, double smoothing) { create(); create(x_,y_,sigma_,smoothing); }
        CubicSpline() { create(); }
        
		void create(const Array1D<T> &x_, const Array1D<T> &y_) { create(x_,y_,Array1D<T>()); }
		void create(const Array1D<T> &x_, const Array1D<T> &y_, const Array1D<T> &sigma_);
		void create(const Array1D<T> &x_, const Array1D<T> &y_, const Array1D<T> &sigma_, double smoothing) { constrainS(smoothing); create(x_,y_,sigma_); }
		void create(const Array1D<T> &x_, const Array1D<T> &y_, double sigma_, double smoothing) { constrainS(smoothing); create(x_,y_,Array1D<T>(x_.dim(),sigma_)); }

        T splint(const T x);
        /** See splint() */
        T operator[](const T x) { return splint(x); }
		
		inline const Array1D<T>& getX() const { return x; }
		inline const Array1D<T>& getY() const { return y; }
		inline const Array1D<T>& getd2YdX2() const { return y2; }
        
        /** Returns the minimum x-value */
        T minx() { return x[0]; }
        /** Returns the maximum x-value */
        T maxx() { return x[x.dim()-1]; }
        
        /** toWave is an alias for toWaveN */
        TWave1D<T,T> toWave(int nPoints) { return toWaveN(nPoints); }
        /** toWave is an alias for toWaveN */
        TWave1D<T,T> toWave(int nPoints, T startX, T endX) { return toWaveN(nPoints,startX,endX); }
        /** See toWaveN(int nPoints, T startX, T endX). Uses minx() and maxx() */
        TWave1D<T,T> toWaveN(int nPoints) { return toWave(nPoints,minx(),maxx()); }
        TWave1D<T,T> toWaveN(int nPoints, T startX, T endX);
        /** See toWaveDT(T deltaX, T startX, T endX). Uses minx() and maxx() */
        TWave1D<T,T> toWaveDT(T deltaX) { return toWave(deltaX,minx(),maxx()); }
        TWave1D<T,T> toWaveDT(T deltaX, T startX, T endX);
        
        /**
        By default, endpoint 2nd derivatives are set to zero. Alternatively, they
        may be set with dydxInitial() and dydxFinal():
        <ul><li>Create an empty (uninitialized) CubicSpline.</li>
        <li>Set dydxInitial() and dydxFinal() as desired.</li>
        <li>Call create(x,y) to do the interpolation.</li></ul>
        
        Note that dydx calls have no effect after passing x and y arrays to
        CubicSpline. Make these calls first, then pass arrays to
        create().
        */
        void dydxInitial(T q) { yp1 = q; yp1_natural = false; }
        /** See dydxInitial(). */
        void dydxFinal(T q) { ypn = q; ypn_natural = false; }
        
        /** Set optimization method.
        <ul><li><b>false</b> <i>(default)</i> optimizes spline for sequential access.</li>
        <li><b>true</b> does no optimization and searches anew for each call.</li></ul>
        */
        bool optimizeRandom(bool q) { return bisect_all = q; }
        bool optimizeRandom() const { return bisect_all; }
        /** Is source data already sorted? <i>(default false)</i> */
        bool sourcePresorted(bool q) { return presorted = q; }
        bool sourcePresorted() const { return presorted; }
		
		/** constrain smoothing amount by some mean deviation */
		T constrainS(const T s_) { return constrain = s_; }
		T constrainS() { return constrain; }
		
		/** constrain smoothing amount statistically by number of points. requires
			standard deviation to be set */
		T constrainN(const T n_) { return n_ > 1 ? constrain = 1.0-sqrt(2/(n_+1)) : 0.0 ; }

    protected:
        Array1D<T> x;
        Array1D<T> y;
		Array1D<T> dy;
        Array1D<T> y2;
        int n; // = x.dim();

		void create();
        
        // controls behavior at endpoints
        T yp1, ypn;
        bool yp1_natural, ypn_natural;

        // used to speed up search process when interpolating
        int klo;
        int khi;
        bool first_splint;
        bool bisect_all;
        bool presorted;
		T constrain;

        /** Makes the 2nd derivative array */
        void make_d2ydx();
		void make_d2ydx_smooth(T s);
		
        /** Sorts data; necessary for spline. Uses smoothsort on the assumption
			that the data is already somewhat ordered. */
        void sort();
};

/** Default/empty spline creation. */
template <class T>
void CubicSpline<T>::create()
{
    first_splint = true;
    optimizeRandom(true);    // default
    sourcePresorted(false);	// default. override for optimization.
    klo = khi = 0;
    yp1_natural = ypn_natural = true;
	constrain = 0.0;
}

/** Creates the spline from the arrays. Usually called from constructior, except in special cases. */
template <class T>
void CubicSpline<T>::create(const Array1D<T> &x_, const Array1D<T> &y_, const Array1D<T> &sigma_)
{
    if(x.dim() != y.dim())
        throw SPLINE_ARRAYS_MUST_BE_SAME_SIZE;

	// don't want to set defaults again
    // create();
    
    // if the caller insists the arrays are ready to go, make copies of them for us. otherwise,
    // sort them and remove duplicates ourselves. sort() copies its finished products into x and y.
	x = x_.copy();
	y = y_.copy();
	dy = sigma_.copy();
	
    if(!sourcePresorted())
        sort();
    
    y2 = Array1D<T>(x.dim());

    n = x.dim();
    make_d2ydx();

    klo = 0;
    khi = x.dim()-1;
}

/**
Creates a TWave1D from the spline with nPpoints points that starts at startX and ends at endX.
*/
template <class T>
TWave1D<T,T> CubicSpline<T>::toWaveN(int nPoints, T startX, T endX)
{
    TWave1D<T,T> a(nPoints);
    a.dt((endX-startX)/(nPoints-1));
    a.startt(startX);

    for(int i=0; i < a.dim()-1; i++)
    {
        a[i] = (*this)[i*a.dt()+a.startt()];
    }
    a[a.dim()-1] = (*this)[endX];	// do the last one with an exact value to avoid possible overrun

    return a;
}

/**
Creates a TWave1D from the spline. Use delta as the interval between points,
and go from startX to endX.
*/
template <class T>
TWave1D<T,T> CubicSpline<T>::toWaveDT(T deltaX, T startX, T endX)
{
	// this function used to call toWaveN, but there was a bad accumulation error,
	// so it now has its own definition

    // it appears that in an especially painful case, a fp precision
    // error can cause startX+(nPoints-1)*deltaX to be a bit larger
    // than endX, which makes splint() crash.	
    int nPoints = (int)((endX-startX)/(deltaX)) + 1;
	if(startX+(nPoints-1)*deltaX > endX)
		nPoints--;

    TWave1D<T,T> a(nPoints);
    a.dt(deltaX);
    a.startt(startX);
	
    for(int i=0; i < a.dim()-1; i++)
    {
        a[i] = (*this)[i*a.dt()+a.startt()];
    }
    a[a.dim()-1] = (*this)[endX];	// do the last one with an exact value to avoid possible overrun
	
    return a;
}


// Cubic splines from Numerical Recipes in C:
// Given arrays x[1..n] and y[1..n] containing a tabulated function, i.e., y_i = fs_i), with s_1 < x_2 < .. < x_N and given values yp1 and ypn for the first derivative of the interpolating function at points 1 and n, respectively, this routine returns an array y2[1..n] that contains the second derivatives of the interpolating function at the tabulated points x_i. If yp1 and/or ypn are equal to 1e10^30 or larger, the routine is singaled to set the corresponding boundary condition for a natrual spline, with zero second derivative on that boundary.
template <class T>
void CubicSpline<T>::make_d2ydx()
{
	if(constrain > 0.0)
	{
		make_d2ydx_smooth(constrain);
		return;
	}

    int i,k;
    T p,qn,sig,un;

    Array1D<T> u(n-1);

    if(yp1_natural)	// The lower boundary condition is set either to be "natural" or else to have a specified first derivative
        y2[0]=u[0]=0.0;
    else {
        y2[0] = -0.5;
        u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
    }
    for(i=1;i<n-1;i++) {	// This is thte decomposition of hte tridiagonal algorithm. y2 and u are used for temporary storage of hte decomposed factors.
        sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
        p=sig*y2[i-1]+2.0;
        y2[i]=(sig-1.0)/p;
        u[i]=(y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
        u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
    }
    if(ypn_natural)	// The upper boundary condition is set either to be "natural" or else to have a specified first derivative
        qn=un=0.0;
    else {
        qn=0.5;
        un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
    }
    y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);
    for(k=n-2;k>=0;k--)	// this is the backsubstitution loop of the tridiagonal algorithm
        y2[k]=y2[k]*y2[k+1]+u[k];
}

/**
Spline interpolation. Return the interpolated y-value corresponding
to the given x-value q.<br/ ><br/ >

TO-DO/BUGS:
<ul><li>check for klo/khi out of bounds before using as indices</li></ul>
*/
template <class T>
T CubicSpline<T>::splint(const T q)
{
	int k;
	T h,b,a;

	if(q < x[0])
    {
//cerr << "q=" << q << ", x[0]=" << x[0] << endl;
	    throw SPLINT_OUT_OF_BOUNDS_LOW;
	}
	if(q > x[x.dim()-1])
    {
//cerr << "q=" << q << ", x[x.dim()-1]=" << x[x.dim()-1] << endl;
//cerr << "difference: " << (q-x[x.dim()-1]) << endl;
        throw SPLINT_OUT_OF_BOUNDS_HIGH;
    }

    // force a full bisection the first run
    if(first_splint || bisect_all)
    {
        klo = 0;
        khi = n-1;
        first_splint = false;
    }

    // first bracket the q
    if(x[khi] < q)
    {
        // move right
        int kstep = 1;
        while(x[klo] > q || x[khi] < q)
        {
            klo = khi;
            khi += kstep;
            kstep += kstep;
            if(khi >= n)
            {
                khi = n-1;
                break;
            }
            
        }
    }
    if(x[klo] > q)
    {
        // move left
        int kstep = -1;
        while(x[klo] > q || x[khi] < q)
        {
            khi = klo;
            klo += kstep;
            kstep += kstep;
            if(klo <= 0)
            {
                klo = 0;
                break;
            }            
        }
    }

    // then bisect
	while (khi-klo > 1)
	{
		k=(khi+klo) >> 1;
		if (x[k] > q)
		    khi=k;
		else
		    klo=k;
	}
	
	h=x[khi]-x[klo];
	if (h == 0.0)
	{
	  //throw x[khi];
            throw BAD_XA_TO_SPLINT;
        }
	a=(x[khi]-q)/h;
	b=(q-x[klo])/h;

    return a*y[klo]+b*y[khi]+((a*a*a-a)*y2[klo]+(b*b*b-b)*y2[khi])*(h*h)/6.0;
}


/**
Follows ALGOL listing in C.H. Reinsch. Numerische Mathematik 10, 177--183 (1967).
Go read the original and drool over a code listing that is typeset quite nicely, despite
some spacing issues and occasional inconsistency in italics.

I've tweaked the meaning of s in this implementation. The N dependence is introduced
automatically, and the input is also squared so it has units of Y deviation. As defined,
it may be used in place of the standard deviation at each point. For example, s = 1.0
and dy[0..N] = 0.001 produces identical results to s = 0.001 and dy[0..N] = 1.0 (default
behavior if you don't pass a dy/sigma).
*/
template <class T>
void CubicSpline<T>::make_d2ydx_smooth(T s)
{
	int MAX_ITER = 100;	// allow lots of iterations
	int n1 = 0;
	int n2 = x.dim()-1;
	s *= x.dim() * s;
	
	if(dy.dim() != x.dim())
		dy = Array1D<T>(x.dim(),1.0);

	// begin smooth
	int i, m1, m2;
	Array1D<T> a(x.dim());
	Array1D<T> b(x.dim());
	Array1D<T> c(x.dim());
	Array1D<T> d(x.dim());
	
	T e, f, f2, g=0, h, p;	// Andy: compiler complains about g being used uninitialized, which is hogwash, but we'll placate it.

	// these arrays are a bit longer. to convince TNT to allow the necessary -1 index,
	// we allocate them long, then take a subarray() to shift the data_ pointer. This
	// will not work if bounds checking is enabled, but that should not be the case
	// when running optimized code.
	Array1D<T> r(x.dim()+2); r = r.subarray(1,x.dim());
	Array1D<T> r1(x.dim()+2); r1 = r1.subarray(1,x.dim());
	Array1D<T> r2(x.dim()+2); r2 = r2.subarray(1,x.dim());
	Array1D<T> t(x.dim()+2); t = t.subarray(1,x.dim());
	Array1D<T> t1(x.dim()+2); t1 = t1.subarray(1,x.dim());
	Array1D<T> u(x.dim()+2); u = u.subarray(1,x.dim());
	Array1D<T> v(x.dim()+2); v = v.subarray(1,x.dim());

	m1 = n1 - 1;
	m2 = n2 + 1;
	r[m1] = r[n1] = r[n2] = r[n2] = r2[m2] = u[m1] = u[n1] = u[n2] = u[m2] = p = 0.0;
	m1 = n1+1;
	m2 = n2 - 1;

	h = x[m1]-x[n1];
	f = (y[m1] - y[n1])/h;

	for(i=m1; i<=m2; i++)
	{
		g = h;
		h = x[i+1] - x[i];
		e = f;
		f = (y[i+1]-y[i])/h;
		a[i] = f-e;
		t[i]=2*(g+h)/3;
		t1[i] = h/3;
		r2[i] = dy[i-1]/g;
		r[i] = dy[i+1]/h;
		r1[i]= -dy[i]/g - dy[i]/h;
	}
	
	for(i=m1; i<=m2; i++)
	{
		b[i] = r[i]*r[i] + r1[i]*r1[i] + r2[i]*r2[i];
		c[i] = r[i]*r1[i+1] + r1[i]*r2[i+1];
		d[i] = r[i]*r2[i+2];
	}
	
	f2 = -s;
	
	int iter;
	for(iter=0; iter<MAX_ITER; iter++)
	{
		// next iteration:
		for(i=m1; i<=m2; i++)
		{
			r1[i-1] = f*r[i-1];
			r2[i-2] = g*r[i-2];
			r[i] = 1.0/(p*b[i] + t[i] - f*r1[i-1] - g* r2[i-2]);
			u[i] = a[i] - r1[i-1]*u[i-1] - r2[i-2]*u[i-2];
			f = p*c[i] + t1[i] - h*r1[i-1];
			g = h;
			h = d[i] * p;
		}
	
		for(i=m2; i >= m1; i--)
		{
			u[i] = r[i]*u[i] - r1[i]*u[i+1] - r2[i]*u[i+2];
		}
		e = h = 0;
	
		for(i=n1; i<=m2; i++)
		{
			g = h;
			h = (u[i+1] - u[i])/(x[i+1] - x[i]);
			v[i] = (h-g) * dy[i] * dy[i];
			e = e + v[i] * (h-g);
		}
	
		g = v[n2] = -h*dy[n2]*dy[n2];
		e = e - g*h;
		g = f2;
		f2 = e*p*p;
		if(f2 >= s || f2 <= g )	// THEN GOTO fin
			break;
	
			f = 0;
		h = (v[m1] - v[n1])/(x[m1] - x[n1]);
	
		for(i=m1; i<=m2; i++)
		{
			g = h;
			h = (v[i+1] - v[i]) / (x[i+1] - x[i]);
			g = h - g - r1[i-1]*r[i-1] - r2[i-2]*r[i-2];
			f = f + g*r[i]*g;
			r[i] = g;
		}
		h = e - p*f;
	
		if(h <= 0) // THEN GOTO fin
			break;

		p = p + (s-f2) / ((sqrt(s/e) + p) * h);	// GOTO next iteration
	
		// comment. Use negative branch of square root, if the sequence of abscissae x[i]
		// is strictly decreasing;
	} // while(1)
	if(iter==MAX_ITER)
	{
		cerr << "Max iterations reached in Reinch algorithm.";
		constrain = 0.0;
		make_d2ydx();
	}
	else // finish up
	{

		// fin:
		for(i=n1; i<=n2; i++)
		{
			a[i] = y[i] - p*v[i];
			c[i] = u[i];
		}
		for(i=n1; i<=m2; i++)
		{
			h = x[i+1] - x[i];
			d[i] = (c[i+1] - c[i])/(3*h);
			b[i] = (a[i+1] - a[i])/h - (h*d[i] + c[i])*h;
		}
	
		// end smooth
	
		// Output of Reinsh's routine produce coefficient arrays a, b, c, d. We wish to keep only
		// the function value
		for(i=1; i<x.dim()-1; i++)
		{
			y[i] = a[i];
			y2[i] = 2*c[i];
		}
	}
}


/**
When a wave comes in, sort() makes sure the x-values are always increasing. 
<ol><li>Look at first and last points. If x_first &gt; x_last, copy the reverse wave.</li>
<li>Sort the wave using smoothsort.</li>
<li>Remove duplicate entries.</li>
<li>Put the finished product in x and y. (Update 12/8/09: happens in-place now)</li></ol>
*/
template <class T>
void CubicSpline<T>::sort()
{
	if(x[1] > x[0])
	{
		SmoothSort<T> s(x);
		s.reorder(y);
		if(dy.dim()==x.dim())
			s.reorder(dy);
	}
	else
	{
		SmoothSort<T,true> s(x);
		s.reverse();
		s.reorder(y);
		if(dy.dim()==x.dim())
			s.reorder(dy);
	}
	
    // now march down the array, removing pairs with duplicate x's
    // current behavior keeps the first of each x it finds and obliterates the rest. there
    // are certainly other (often more appropriate) ways to handle this situation.
    int i=1;
    int dupcount = 0;
    T lastval = x[0];
    while((i+dupcount)<x.dim())
    {
        if(x[i+dupcount]==lastval)
        {
            dupcount++;
//cerr << "At i=" << i << ", dupcount=" << dupcount << ", lastval=" << lastval << endl;
        }
        else
        {
            if(dupcount)
            {
                x[i] = x[i+dupcount];
                y[i] = y[i+dupcount];
//cerr << "Copying from source[" << i+dupcount << "] to dest[" << i << "]" << endl;
            }
            lastval = x[i];
            i++;
        }
    }
    
    // put the arrays in their proper locations, masking off just the relevant parts.
    x = x.subarray(0,x.dim()-dupcount-1);
    y = y.subarray(0,y.dim()-dupcount-1);
}

template <class T>
Array1D<T> smoothSpline(const Array1D<T>& x, const Array1D<T>& y, T sigma, T alpha)
{
cerr << "smoothSpline: dim(" << x.dim() << "," << y.dim() << "), sigma=" << sigma << ", alpha=" << alpha << endl;
	CubicSpline<T> cs(x,y,sigma,alpha);
	return cs.getY();
}

template <class T>
T arrayMean(const Array1D<T>& arr, int lb, int ub)
{
	if(arr.dim()==0) return 0;
	T accum = 0;
	for(int i=MAX<T>(0,lb); i < MIN<T>(arr.dim(),ub); i++)
	{
		accum += arr[i];
	}
	return accum/arr.dim();
}

template <class T>
inline T arrayMean(const Array1D<T>& arr)
{
	return arrayMean(arr,0,arr.dim());
}

template<class T>
T stdDev(const Array1D<T>& arr, const T avg, int lb, int ub)
{
	if(arr.dim()==0) return 0;
	
	T accum = 0;
	for(int i=MAX<T>(0,lb); i < MIN<T>(arr.dim(),ub); i++)
	{
		T q = arr[i] - avg;
		accum += q*q;
	}
	return sqrt(accum/arr.dim());	
}

template<class T>
T stdDev(const Array1D<T>& arr, int lb, int ub)
{
	if(arr.dim()==0) return 0;
	T avg = arrayMean<T>(arr,lb,ub);
	
	return stdDev<T>(arr,avg,lb,ub);
}

template<class T>
inline T stdDev(const Array1D<T>& arr, const T avg)
{
	return stdDev(arr,avg,0,arr.dim());
}

template<class T>
inline T stdDev(const Array1D<T>& arr)
{
	return stdDev(arr,0,arr.dim());
}

template<class T>
T stdErr(const Array1D<T>& a, const Array1D<T>& b, int lb, int ub)
{
	if(a.dim()==0 || a.dim() != b.dim()) return 0;
	
	T accum = 0;
	for(int i=MAX<T>(0,lb); i < MIN<T>(a.dim(),ub); i++)
	{
		T q = a[i] - b[i];
		accum += q*q;
	}
	return sqrt(accum/a.dim());	
}

template<class T>
inline T stdErr(const Array1D<T>& a, const Array1D<T>& b)
{
	return stdErr(a,b,0,a.dim());
}

template <class T>
Array1D<T> adaptiveSmoothSpline(const Array1D<T>& x, const Array1D<T>& y, T sigma, T& alpha)
{
//cerr << "Starting adaptiveSmoothSpline with sigma = " << sigma << endl;
	const int MAXITER = 50;
	const T MINALPHA = 0.01;
	const T MAXALPHA = 1e8;
	const int LB = x.dim()/10;
	const int UB = x.dim() - LB;
	const T EPSILON = 0.1; // fractional tolerance
	const T logsigma = log10(sigma);
	T err;
	bool isBracketed;
	
	T p0 = 0;	// starting points in log10 space
	T p1 = 4;
	T p2;
	T yp2;
	Array1D<T> wavep2;
	
	Array1D<T> wavep0 = smoothSpline(x,y,sigma,pow(10,p0));
	for(int j=0; j<y.dim(); j++) { wavep0[j] -= y[j]; }
	{ SmoothSort<T> s(wavep0); }
	T yp0 = log10(stdDev<T>(wavep0,0.0,LB,UB)) - logsigma;

	Array1D<T> wavep1 = smoothSpline(x,y,sigma,pow(10,p1));
	for(int j=0; j<y.dim(); j++) { wavep1[j] -= y[j]; }
	{ SmoothSort<T> s(wavep1); }
	T yp1 = log10(stdDev<T>(wavep1,0.0,LB,UB)) - logsigma;
	
//cerr << " p0 = " << p0 << ", yp0 = " << yp0 << ", 10^p[" << 0 << "] = " << pow(10,p0) << ", 10^yp[" << 0 << "] = " << pow(10,yp0) << endl;
//cerr << " p1 = " << p1 << ", yp1 = " << yp1 << ", 10^p[" << 1 << "] = " << pow(10,p1) << ", 10^yp[" << 1 << "] = " << pow(10,yp1) << endl;

	int i=0;
	for( ; i<MAXITER; i++)
	{
		isBracketed = (yp1*yp0 < 0);	// if we have bracketed the root, bisection can stabilize convergence
		
//cerr << "At beginning of iteration, yp1 = " << yp1 << endl;
//cerr << "p2 = " << p1 << " - (" << yp1 << "*(" << p1 << " - " << p0 << ")) / (" << yp1 << " - " << yp0 << ")" << endl;
		p2 = p1 - (yp1*(p1-p0)) / (yp1 - yp0);
		
		// bounds check
		if(p2 > log10(MAXALPHA))
		   p2 = log10(MAXALPHA);
		if(p2 < log10(MINALPHA))
		   p2 = log10(MINALPHA);

		
		wavep2 = smoothSpline(x,y,sigma,pow(10,p2));
		for(int j=0; j<y.dim(); j++) { wavep2[j] -= y[j]; }
		{ SmoothSort<T> s(wavep2); }
		yp2 = log10(stdDev<T>(wavep2,0.0,LB,UB)) - logsigma;
//cerr << "Smoothing with (" << sigma << ", " << pow(10,p2) << ") => " << pow(10,yp2+logsigma) << endl;
				   
		if(isBracketed && ((p2 > p1 && p2 > p0) || (p2 < p1 && p2 < p0) || (p1*p2>=0)) || (p2 > log10(MAXALPHA)) || (p2 < log10(MINALPHA)))	// always accept if we haven't bracketed the root
		{
//cerr << "Trial " << i << " failed. Bisecting." << endl;
			p2 = (p0+p1)/2;
			wavep2 = smoothSpline(x,y,sigma,pow(10,p2));
			for(int j=0; j<y.dim(); j++) { wavep2[j] -= y[j]; }
			{ SmoothSort<T> s(wavep2); }
			yp2 = log10(stdDev<T>(wavep2,0.0,LB,UB)) - logsigma;			
//cerr << "Smoothing with (" << sigma << ", " << pow(10,p2) << ") => " << pow(10,yp2+logsigma) << endl;
			
			if(yp2*yp0<0)
			{
				p1 = p2;
				wavep1 = wavep2;
				yp1 = yp2;
			}
			else
			{
				p0 = p2;
				wavep0 = wavep2;
				yp0 = yp2;
			}

		}
		else
		{
			p0 = p1;
			wavep0 = wavep1;
			yp0 = yp1; 
			
			p1 = p2;
			wavep1 = wavep2;
			yp1 = yp2;
		}
		
		err = abs<T>((pow(10,yp2+logsigma)-sigma)/sigma);		

//cerr << "End of iteration: p1 = " << p1 << ", yp1 = " << yp1 << ", 10^p[" << i+2 << "] = " << pow(10,p1) << ", 10^(yp[" << i+2 << "]+logsigma) = " << pow(10,yp1+logsigma) << ", err = " << err << endl;
//cerr << "    (p0 = " << p0 << ", yp0 = " << yp0 << ")" << endl;
		if(err < EPSILON) break;
		
//cerr << "At end of iteration, yp1 = " << yp1 << endl;
	}
	if(i==MAXITER)
//		throw NO_CONVERGENCE_IN_ADAPTIVE_SMOOTH_SPLINE;
		return y;	// just return the old one if we can't deal with it
		
//cerr << "Adaptive Smoothing convergence: alpha = " << pow(10,p2) <<  ", sigma = " << pow(10,yp2+logsigma)  << endl;

	alpha = pow(10,p2);
	return smoothSpline(x,y,sigma,pow(10,p2));
}

template <class T>
Array1D<T> adaptiveSmoothSpline(const Array1D<T>& x, const Array1D<T>& y, T sigma)
{
	T alpha = 0;
	return adaptiveSmoothSpline(x,y,sigma,alpha);
}


#endif
