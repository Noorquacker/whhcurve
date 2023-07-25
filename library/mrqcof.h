#ifndef _MRQCOF_H
#define _MRQCOF_H

#include "tnt/tnt.h"
#include "tntnr.h"

using namespace TNT;

extern const char* MRQCOF_DIM_MISMATCH;
extern const char* FGAUSS_DIM_MISMATCH;


template <class T>
void mrqcof(
    Array1D<T> x,
    Array1D<T> y, 
    Array1D<T> sig,
    Array1D<T> a,
    Array1D<int> ia,
    Array2D<T> cofalpha,
    Array1D<T> cofbeta,
    T& chisq,
    T (*funcs)(T, const Array1D<T>, Array1D<T>))
{
	int i, j, k, l, m, mfit=0;
	T ymod,wt,sig2i,dy;
	
	if(x.dim() != y.dim() || y.dim() != sig.dim()
	    || a.dim() != ia.dim())
	        throw MRQCOF_DIM_MISMATCH;

	int ndata = x.dim();
	int ma = a.dim();

	Array1D<T> dyda(ma);
	
    for(j=0; j<ma; j++)
        if(ia[j]) mfit++;
    for(j=0; j<mfit; j++)
    {
        for(k=0;k<=j;k++)
            cofalpha[j][k] = 0.0;
		cofbeta[j] = 0.0;
	}
	chisq = 0.0;
	for(i=0; i<ndata; i++)
	{
        ymod = (*funcs)(x[i],a,dyda);
        sig2i = 1.0/(sig[i]*sig[i]);
        dy = y[i] - ymod;
        for(j=-1, l=0; l<ma; l++)
        {
            if(ia[l])
            {
                wt=dyda[l]*sig2i;
                for(j++,k=-1,m=0; m<=l; m++)
                    if(ia[m]) cofalpha[j][++k] += wt * dyda[m];
                cofbeta[j] += dy * wt;
            }
		}
		chisq += dy*dy*sig2i;
	}
	for(j=1; j<mfit; j++)
		for(k=0; k<j-1; k++)
		    cofalpha[k][j] = cofalpha[j][k];
}

template <class T>
T fgauss(T x, const Array1D<T>& a, Array1D<T>& dyda)
{
    int na;
    #ifdef TNT_NR_DIM_CHECK
    if(a.dim() != dyda.dim())
        throw FGAUSS_DIM_MISMATCH;
    #endif
    
    int i;
    T fac, ex, arg, y = 0.0;
    
    for(i=0; i<na-1; i+=3)
    {
        arg = (x-a[i+1])/a[i+2];
        ex=exp(-arg*arg);
        fac=a[i]*ex*2.0*arg;
        y += a[i]*ex;
        dyda[i]=ex;
        dyda[i+1]=fac/a[i+2];
        dyda[i+2]=fac*arg/a[i+2];
    }
    return y;
}

#endif
// _MRQCOF_H
