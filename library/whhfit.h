#ifndef _WHH_FIT_H 
#define _WHH_FIT_H

#include <cstdlib>
#include <math.h>
#include "tnt/tnt.h"
#include "mrqmin.h"

using namespace TNT;

// from the theory of Werthamer, Helfand, and Hohenberg (1966)

template <class T>
class WHHSolver
{
    public:
	// input seed values to constructor
	WHHSolver(T t_, T h_, T alpha_, T lambda_so_)
	{
	    t = t_;
	    h = h_;
	    alpha = alpha_;
	    lambda_so = lambda_so_;

	    ITERMAX = 1000;
	    good_values = false;
	}
	WHHSolver() : WHHSolver(0,0,0,0) { }

	// solve the whh. see secant_solve or newton_solve for param values	
	void solve(int param, bool reset = false)
	{
	    if(reset) good_values = false;

	    if(good_values)
		newton_solve(param,0.001,0.001);
	    else
		secant_solve(param,0.0,10.0,1e-8); 
	}

	// solve by the Secant method starting from a window of pmin
	// and pmax
	//
	// param indicates which of the 4 WHH parameters should be
	// adjusted
	void secant_solve(int param, T pmin, T pmax, T tol)
	{
//cerr << "secant_solve" << endl;
	    T x, xn, y, yn, xnn, ynn;
	    x = pmin;
	    xn = pmax;
//cerr << "getting initial values" << endl;
	    switch(param)
	    {
		case 0:	// adjust t
		    y = WHHEval(x, h, alpha, lambda_so);
		    yn = WHHEval(xn, h, alpha, lambda_so);
		    t = x;
		    break;
		case 1:	// adjust h
		    y = WHHEval(t, x, alpha, lambda_so);
		    yn = WHHEval(t, xn, alpha, lambda_so);
		    h = x;
		    break;
		case 2:	// adjust alpha
		    y = WHHEval(t, h, x, lambda_so);
		    yn = WHHEval(t, h, xn, lambda_so);
		    alpha = x;
		    break;
		case 3:	// adjust lambda_so
		    y = WHHEval(t, h, alpha, x);
		    yn = WHHEval(t, h, alpha, xn);
		    lambda_so = x;
		    break;
		default:
		    throw "WHHSolve: param incorrect in secant_solve";
	    }
	    
	    // loop until (x-xn) is small (or we hit ITERMAX)
	    int i;
	    for(i=2; (i<ITERMAX) && (abs(x-xn)>tol); i++)
	    {
//cerr << "iteration " << i << ", x=" << x << ", xn=" << xn << endl;
		    
		    // first priority is to take points of opposite sign (bracketing points)
		    // second priority is to take the closest point
		    // last priority is to take the newest point
		    //
		    // so case-by-case: both yn and ynn are opposite sign to y
		    if(y*yn<0 && y*ynn<0)
		    {
			// take the closer one
			if(abs(ynn)<abs(yn))
			{
			    xn = xnn;
			    yn = ynn;
			} // otherwise, leave it alone, and _nn will drop off
		    }
		    else if(y*ynn<0)
		    {
			xn = xnn; yn = ynn;
		    }

		    T temp = x - y*(x - xn)/(y - yn);
		    xnn = xn; ynn = yn;
		    xn = x; yn = y;
		    x = temp;

		    switch(param)
		    {
			case 0:	// adjust t
			    y = WHHEval(x, h, alpha, lambda_so);
			    t = x;
			    break;
			case 1:	// adjust h
			    y = WHHEval(t, x, alpha, lambda_so);
			    h = x;
			    break;
			case 2:	// adjust alpha
			    y = WHHEval(t, h, x, lambda_so);
			    alpha = x;
			    break;
			case 3:	// adjust lambda_so
			    y = WHHEval(t, h, alpha, x);
			    lambda_so = x;
			    break;
			default:
			    throw "WHHSolve: param incorrect in secant_solve";
		    }
	    } // back for another iteration
//cerr << "done. x=" << x << ", xn=" << xn << endl;
	    if(i>=ITERMAX)
		throw "No convergence in WHHSolve::secant_solve.";
	    good_values = true;
	}

	// solve by Newton-Raphson iteration starting from the
	// point defined by this object
	//
	// param indicates which of the 4 WHH parameters should be
	// adjusted
	void newton_solve(int param, T dx, T tol)
	{
//cerr << "newton_solve" << endl;
	    T x, xn, y, yn, lastval;
	    int i;
	    for(i=0;(i<ITERMAX/2) && (abs(x-lastval)>tol);i++) 
	    {
//cerr << "iteration " << i << endl;
		    switch(param)
		    {
			case 0: // t 
			    x = t + dx/2.0;
			    xn = t - dx/2.0;
			    y = WHHEval(x, h, alpha, lambda_so);
			    yn = WHHEval(xn, h, alpha, lambda_so);
			    lastval = t = x - y*(x - xn)/(y - yn);
			    break;
			case 1: // t 
			    x = h + dx/2.0;
			    xn = h - dx/2.0;
			    y = WHHEval(t, x, alpha, lambda_so);
			    yn = WHHEval(t, xn, alpha, lambda_so);
			    lastval = h = x - y*(x - xn)/(y - yn);
			    break;
			case 2: // t 
			    x = alpha + dx/2.0;
			    xn = alpha - dx/2.0;
			    y = WHHEval(t, h, x, lambda_so);
			    yn = WHHEval(t, h, xn, lambda_so);
			    lastval = alpha = x - y*(x - xn)/(y - yn);
			    break;
			case 3: // t 
			    x = lambda_so + dx/2.0;
			    xn = lambda_so - dx/2.0;
			    y = WHHEval(t, h, alpha, x);
			    yn = WHHEval(t, h, alpha, xn);
			    lastval = lambda_so = x - y*(x - xn)/(y - yn);
			    break;
			default:
			    throw "WHHSolve: param incorrect in secant_solve";
		    }
	    } 
	    if(i>=ITERMAX/2)
		throw "No convergence in WHHSolve::newton_solve.";
	    good_values = true;
	}

	static T WHHEval(T t_, T h_, T alpha_, T lambda_so_, int NMAX=10000)
	{
//cerr << "WHHEval: " << t_ << "\t" << h_ << "\t" << alpha_ << "\t" << lambda_so_ << "\t" << NMAX;
	    T accum = 0.0;
	    for(int n=-NMAX; n<=NMAX; n++)
	    {	
		accum += 1.0/abs(2*n+1) - 1.0/(abs(2*n+1) + h_/t_ + alpha_*alpha_*h_*h_/t_/t_/(abs(2*n+1) + (h_ + lambda_so_) / t_));
	    }
//cerr << "\tReturning: " << (accum - log(1.0/t_)) << endl;
	    return accum - log(1.0/t_);
	}

    // protection removed to allow lightweight tweaking
    //protected:
	T t;
	T h;
	T alpha;
	T lambda_so;

    protected:
	int ITERMAX;
	bool good_values;
};



template <class T>
class WHHFit : public LevenbergMarquardtFit<T>
{
    public:
    
    WHHFit(unsigned xndatamax, unsigned xmamax)
        : LevenbergMarquardtFit<T>::LevenbergMarquardtFit(xndatamax,xmamax) { }
    
    T start(const Array1D<T>& qx, const Array1D<T>& qy, const Array1D<T>& qsig, const Array1D<T>& guessa)
    {
        return LevenbergMarquardtFit<T>::start(qx,qy,qsig,guessa, WHHFit<T>::whhfunc);
    }
    const Array1D<T> fit(const Array1D<T>& qx, const Array1D<T>& qy, const Array1D<T>& qsig, const Array1D<T>& guessa)
    {
        return LevenbergMarquardtFit<T>::fit(qx,qy,qsig,guessa, WHHFit<T>::whhfunc);
    }
//PUT PROTECTION BACK    
//    protected:
    
    static T whhfunc(T x, const Array1D<T> a, Array1D<T> dyda);
};

/** Here's the meat of the WHHFit class. */
template <class T>
    T WHHFit<T>::whhfunc(T x, const Array1D<T> a, Array1D<T> dyda)
    {
//	// commented out from gaussianfunc
        T fac,ex,arg,y;
//        
//        arg=(x-a[1])/a[2];
//        ex=exp(-arg*arg);
//        fac=a[0]*ex*2.0*arg;
//        y = a[0]*ex;
//        
//        dyda[0]=ex;
//        dyda[1]=fac/a[2];
//        dyda[2]=fac*arg/a[2];
	// a[0] is t
	// a[1] is h
	// a[2] is alpha
	// a[3] is lambda_so 
	const int NRANGE = 100;
	

	return y;
    }

#endif // _WHH_FIT_H
