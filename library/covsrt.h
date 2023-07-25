#ifndef _COVSRT_H
#define _COVSRT_H

#include "tnt/tnt.h"
#include "tntnr.h"

using namespace TNT;


extern const char* COVSRT_DIM_MISMATCH;

// covar is the covariance matrix we would like to expand
// ia[j] is 0 for parameters held fixed and nonzero for fitted ones
// mfit is the number of fitted parameters (number of nonzeros in ia)
template <class T>
void covsrt(Array2D<T> covar, Array1D<int> ia, int mfit)
{
    int i, j, k;
    int ma = covar.dim1();
    
    #ifdef TNT_NR_DIM_CHECK    
    if(covar.dim1() != covar.dim2() || covar.dim2() != ia.dim())
        throw COVSRT_DIM_MISMATCH;
    #endif
    
    for(i=mfit-1; i<ma; i++)
        for(j=0; j<=i; j++) covar[i][j] = covar[j][i] = 0.0;

    k = mfit-1;
    for(j=ma-1; j>=0; j--)
    {
        if(ia[j])
        {
            for(i=0; i<ma; i++) SWAP<T>(covar[i][k],covar[i][j]);
            for(i=0; i<ma; i++) SWAP<T>(covar[k][i],covar[j][i]);
            k--;
        }        
    }
}

#endif
// _COVSRT_H
