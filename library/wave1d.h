#ifndef _WAVE1D_H
#define _WAVE1D_H

#include "tnt/tnt.h"

using namespace TNT;

// defined in cpp file
extern const char* DIMENSION_DISAGREEMENT;
extern const char* CANNOT_CONSTRUCT_FROM_WAVE1D;
extern const char* FOURIER_MUST_BE_POWER_OF_TWO;
extern const char* INCOMPATIBLE_TWAVE_FWAVE;
int nextPowerOfTwo(int n);

template <class T, class TimeType>
class Wave1D;

template <class T, class TimeType>
class TWave1D;

template <class T, class TimeType>
class FWave1D;

template <class T>
void realft(T data[],int n, int isign);

template <class T>
void four1(T data[], int nn, int isign);

/**
<h1>Example of building functions that operate on waves</h1>
<pre>REAL_FUNCTION(sin)
REAL_FUNCTION(cos)

void doSomething()
{
    double PI = 3.141592654;
    TReal1D x(65536);
    x.mint(-0.5);
    x.maxt(0.5);
    x = x.t();

    // VERY important to express data types correctly here.
    // 2 is not the same as 2.0, as they have different types
    TReal1D sin2x = sin(2.0*x);
    TReal1D cos2x = cos(2.0*x);
}</pre>

*/
template <class T, class TimeType>
class Wave1D : public Array1D<T>
{
    public:
    
    Wave1D() : Array1D<T>::Array1D() {initTime();}
    Wave1D(int n) : Array1D<T>::Array1D(n) {initTime();}
    Wave1D(int n, const T &a) : Array1D<T>::Array1D(n,a) {initTime();}
    Wave1D(int n,  T *a) : Array1D<T>::Array1D(n,a) {initTime();}
    Wave1D(const Array1D<T> &A) : Array1D<T>::Array1D(A) {initTime();}
    Wave1D(const Wave1D &A) : Array1D<T>::Array1D(A) {copyTime(A);}

    Wave1D<T,TimeType> & ref(const Wave1D<T,TimeType> &A) {Array1D<T>::ref(A); copyTime(A); return (*this);}
    Wave1D<T,TimeType> copy() const;
    Wave1D<T,TimeType> & inject(const Wave1D<T,TimeType> &A) {Array1D<T>::inject(A); copyTime(A); return (*this);}

    /** Low-level access function.
    Use with care. Sometimes you just need to be able to address a block from a pointer.
    */
    T* data() { return &((*this)[0]); }

    /** Get time interval between wave points */
    inline TimeType dt() const { return dt_; }
    /** Set time interval between wave points */
    inline TimeType dt(const TimeType &dt_new) { return dt_ = dt_new; }
    /** Get time corrosponding to array index 0 */
    inline TimeType startt() const { return startt_; }
    /** Set time corrosponding to array index 0 */
    inline TimeType startt(const TimeType &startt_new) { return startt_ = startt_new; }
    
    /** Copy time information from wave A to this object */
    virtual void copyTime(const Wave1D<T,TimeType> &A) {dt(A.dt()); startt(A.startt());}
    /** Check to see if wave A has the same time characteristics as this. Arithmetic may only be performed on two arrays that match according to this function. */
    virtual bool isCompatible(const Wave1D<T,TimeType> &A) const {return (dt()==A.dt() && startt()==A.startt() && Array1D<T>::dim()==A.dim()); }

    protected:
    
    /** Default time information initialization */
    virtual inline void initTime() { dt_=1; startt_=0; }

    TimeType dt_;
    TimeType startt_;
};

/** Returns a duplicate copy of this object. See also TNT::Array1D<T>::copy() */
template <class T, class TimeType>
Wave1D<T,TimeType> Wave1D<T,TimeType>::copy() const
{
    Wave1D<T,TimeType> A(Array1D<T>::copy());
    A.copyTime((*this));
    
    return A;
}

#define REAL_FUNCTION(fn) inline TReal1D fn(const TReal1D &x) { return byMember<double,double,double,fn>(x); } inline FReal1D fn(const FReal1D &x) { return byMember<double,double,double,fn>(x); }

typedef TWave1D<double,double> TReal1D;
typedef FWave1D<double,double> FReal1D;

#endif
