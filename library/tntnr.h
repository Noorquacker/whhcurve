#ifndef _TNT_NR_H
#define _TNT_NR_H

#include "tnt/tnt.h"
#include <algorithm>

using namespace TNT;

// Bounds checking on dimensions as they are passed between functions
//#define TNT_NR_DIM_CHECK

// Same thing for TNT array indexing
//#define TNT_BOUNDS_CHECK

//#ifndef __OBJC__
//template <class T>
//inline T MIN(T a, T b) { return a > b ? b : a; }
//
//template <class T>
//inline T MAX(T a, T b) { return a > b ? a : b; }
//#endif

#define MAX std::max
#define MIN std::min

template <class T>
inline void SWAP(T& a, T& b)
{
    T temp = a;
    a = b;
    b = temp;
}

template <class T>
inline void SHFT(T& a, T& b, T& c, const T& d)
{
    a=b;
    b=c;
    c=d;
}

template <class T>
inline T SIGN(const T& a, const T& b)
{
     return b >= 0.0
        ? (a>=0 ? a : -a)
        : (a>0 ? -a : a);
}

template <class T>
inline T abs(const T& a)
{
    return (a>=0.0 ? a : -a);
}

/*
// these are already declared by TNT
template <class T>
ostream& operator << (ostream& out, Array1D<T> arr)
{
    for(unsigned i=0; i < arr.dim(); i++)
    {
        cout << arr[i] << " ";
    }
    cout << endl;
}
template <class T>
ostream& operator << (ostream& out, Array2D<T> arr)
{
    for(unsigned i=0; i < arr.dim1(); i++)
    {
        for(unsigned j=0; i < arr.dim2(); j++)
        {
            cout << arr[i][j] << " ";
        }
        cout << endl;
    }
}
*/
#endif
