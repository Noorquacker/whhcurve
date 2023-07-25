#ifndef _GAUSSJ_H
#define _GAUSSJ_H

#include <math.h>
#include "tnt/tnt.h"
#include "tntnr.h"

using namespace TNT;


extern const char* GAUSSJ_DIM_MISMATCH;
extern const char* GAUSSJ_SINGULAR_1;
extern const char* GAUSSJ_SINGULAR_2;


template <class T>
void gaussj(Array2D<T> a, Array2D<T> b)
{
    int n, m;
    
    #ifdef TNT_NR_DIM_CHECK    
    if(a.dim1() != a.dim2() || a.dim2() != b.dim1())
        throw GAUSSJ_DIM_MISMATCH;
    #endif
    
    n = a.dim1();
    m = b.dim2();

    Array1D<int> indxc(n);
    Array1D<int> indxr(n);
    Array1D<int> ipiv(n);

    int i, icol, irow, j, k, l, ll;
    T big, dum, pivinv;

    for(j=0; j<n; j++)
        ipiv[j]=0;
        
    for(i=0; i<n; i++)
    {
        big=0.0;
        for(j=0; j<n; j++)
        {
            if(ipiv[j] != 1)
            {
                for (k=0; k<n; k++)
                {
                    if(ipiv[k] == 0)
                    {
                        if(fabs(a[j][k]) >= big)
                        {
                            big=fabs(a[j][k]);
                            irow=j;
                            icol=k;
                        }
                    }
                    else if(ipiv[k] > 1)
                        throw GAUSSJ_SINGULAR_1;
                }
            }
        }
        ++(ipiv[icol]);
        if(irow != icol)
        {
            for(l=0;l<n;l++) SWAP<T>(a[irow][l],a[icol][l]);
            for(l=0;l<m;l++) SWAP<T>(b[irow][l],b[icol][l]);
        }

        indxr[i]=irow;
        indxc[i]=icol;
            
        if(a[icol][icol] == 0.0)
            throw GAUSSJ_SINGULAR_2;

        pivinv=1.0/a[icol][icol];
        a[icol][icol]=1.0;

        for(l=0;l<n;l++) a[icol][l] *= pivinv;
        for(l=0;l<m;l++) b[icol][l] *= pivinv;
        for(ll=0;ll<n;ll++)
        {
            if (ll != icol)
            {
                dum=a[ll][icol];
                a[ll][icol]=0.0;
                for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
                for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
            }
        }
    }
    for(l=n-1;l>=0;l--)
    {
        if(indxr[l] != indxc[l])
            for(k=0;k<n;k++)
                SWAP<T>(a[k][indxr[l]],a[k][indxc[l]]);
    }
}

#endif
// _GAUSSJ_H
