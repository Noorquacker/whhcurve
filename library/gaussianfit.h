#ifndef _GAUSSIAN_FIT_H
#define _GAUSSIAN_FIT_H

#include <math.h>
#include "mrqmin.h"

using namespace TNT;


template <class T>
class GaussianFit : public LevenbergMarquardtFit<T>
{
    public:
    
    GaussianFit(unsigned xndatamax, unsigned xmamax)
        : LevenbergMarquardtFit<T>::LevenbergMarquardtFit(xndatamax,xmamax) { }
    
    T start(const Array1D<T>& qx, const Array1D<T>& qy, const Array1D<T>& qsig, const Array1D<T>& guessa)
    {
        return LevenbergMarquardtFit<T>::start(qx,qy,qsig,guessa, GaussianFit<T>::gaussianfunc);
    }
    const Array1D<T> fit(const Array1D<T>& qx, const Array1D<T>& qy, const Array1D<T>& qsig, const Array1D<T>& guessa)
    {
        return LevenbergMarquardtFit<T>::fit(qx,qy,qsig,guessa, GaussianFit<T>::gaussianfunc);
    }
//PUT PROTECTION BACK    
//    protected:
    
    static T gaussianfunc(T x, const Array1D<T> a, Array1D<T> dyda);
    
};

/** Here's the meat of the GaussianFit class. */
template <class T>
    T GaussianFit<T>::gaussianfunc(T x, const Array1D<T> a, Array1D<T> dyda)
    {
        T fac,ex,arg,y;
        
        arg=(x-a[1])/a[2];
        ex=exp(-arg*arg);
        fac=a[0]*ex*2.0*arg;
        y = a[0]*ex;
        
        dyda[0]=ex;
        dyda[1]=fac/a[2];
        dyda[2]=fac*arg/a[2];
        
        return y;
    }

#endif
