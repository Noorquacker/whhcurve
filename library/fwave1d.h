/** \file
	Waves in Fourier space
*/
#ifndef _FWAVE1D_H
#define _FWAVE1D_H

#include "tnt/tnt.h"
#include "wave1d.h"
#include "tntnr.h"

using namespace TNT;

/** Waves in Fourier space.
*/
template <class T, class TimeType>
class FWave1D : public Wave1D<T,TimeType>
{
    public:
    
    FWave1D() : Wave1D<T,TimeType>::Wave1D() { initTime(0); }
    FWave1D(int n) : Wave1D<T,TimeType>::Wave1D(n) { initTime(this->dim()); }
    FWave1D(int n, const T &a) : Wave1D<T,TimeType>::Wave1D(n,a) { initTime(this->dim()); }
    FWave1D(int n,  T *a) : Wave1D<T,TimeType>::Wave1D(n,a) { initTime(this->dim()); }
    FWave1D(const Array1D<T> &A) : Wave1D<T,TimeType>::Wave1D(A) { initTime(this->dim()); }
    FWave1D(const Wave1D<T,TimeType> &A) { throw CANNOT_CONSTRUCT_FROM_WAVE1D; }
    FWave1D(const TWave1D<T,TimeType> &A) { (*this) = A; }
    FWave1D(const FWave1D &B) { (*this) = B; }

    inline TimeType minf() const { return 0; }
    inline TimeType maxf() const { return 1/(2*(*this).dt()); }
    inline TimeType maxf(const TimeType &max) { (*this).dt(1.0/(2.0*max)); return max; }
    inline int originalLength() const { return original_length_; }
    inline int originalLength(int ol) { return original_length_ = ol; }
    inline bool isPeriodic() const { return is_periodic_; }
    inline bool isPeriodic(const bool ip) { return is_periodic_ = ip; }
    
    FWave1D<TimeType,TimeType> f() const;
    TimeType f(int i) const { if(i==0) return minf(); if(i==1) return maxf(); return minf() + (i)*(maxf()-minf())/(*this).dim(); }
	TWave1D<T,TimeType> magsqt() const;
	TWave1D<T,TimeType> magsqtf() const;
	
    // duplicate the assignment functions with proper return types
    FWave1D<T,TimeType> & ref(const FWave1D<T,TimeType> &A) { Wave1D<T,TimeType>::ref(A); return (*this); }
    FWave1D<T,TimeType> copy() const;
    FWave1D<T,TimeType> & inject(const FWave1D<T,TimeType> &A) {Wave1D<T,TimeType>::ref(A); return (*this);}

    virtual FWave1D<T,TimeType>& operator = (const FWave1D<T,TimeType>& A);
    virtual FWave1D<T,TimeType>& operator = (const TWave1D<T,TimeType>& y);

    virtual void copyTime(const FWave1D<T,TimeType> &A) { Wave1D<T,TimeType>::copyTime(A); originalLength(A.originalLength()); isPeriodic(A.isPeriodic()); }

    protected:
    virtual void initTime(int ol);
 
    int original_length_;
    bool is_periodic_;
};

/** get a wave whose values are the frequency values at each point
*/
template <class T, class TimeType>
FWave1D<TimeType,TimeType> FWave1D<T,TimeType>::f() const
{
    FWave1D<TimeType,TimeType> freqs((*this).dim());
    freqs.copyTime(*this);

    freqs[0] = minf();
    freqs[1] = maxf();
    for(int i=2;i<(*this).dim();i+=2)
        freqs[i] = freqs[i+1] = minf() + (i)*(maxf()-minf())/(*this).dim();
    return freqs;

}

/**	Aside from the types, this looks just like Wave1D::copy(). Unfortunately, we can't declare
	a return type virtual.
*/
template <class T, class TimeType>
FWave1D<T,TimeType> FWave1D<T,TimeType>::copy() const
{
    FWave1D<T,TimeType> A(Array1D<T>::copy());
    A.copyTime((*this));
    
    return A;
}


//#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

/** From Numerical Recipes in C
	Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as 1;
	or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform, if isign
	is input as -1. data is a complex array of length nn or, equivalently, a real array
	of length 2*nn. nn MUST be an integer power of 2 (this is not checked for!). */
template <class T>
void four1(T data[], int nn, int isign)
{
    const T PI = 3.141592654;

    int n,mmax,m,j,istep,i;
    T wtemp,wr,wpr,wpi,wi,theta;
    T tempr,tempi;
    
    n=nn << 1;
    j=1;
    for (i=1;i<n;i+=2) {	// this is the bit-reversal section of the routine.
        if (j > i) {
            SWAP<T>(data[j],data[i]);	// exchange the two complex numbers
            SWAP<T>(data[j+1],data[i+1]);
        }
        m=n >> 1;
        while (m >= 2 && j > m) {
            j -= m;
            m >>= 1;
       }
        j += m;
    }
    
    // Here begins the Danielson-Lanczos section of the routine.
    mmax=2;
    while(n > mmax) {
        istep=mmax << 1;
        theta=isign*(2*PI/mmax);	// initialize the trigonometric recurrance
        wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp;
        wpi=sin(theta);
        wr=1.0;
        wi=0.0;
        for (m=1;m<mmax;m+=2) {	// here are the two nested loops
            for (i=m; i<=n;i+=istep) {
                j=i+mmax;	// This is the Danielson-Lanczos formula
                tempr=wr*data[j]-wi*data[j+1];
                tempi=wr*data[j+1]+wi*data[j];data[j]=data[i]-tempr;
                data[j+1]=data[i+1]-tempi;
                data[i] += tempr;
                data[i+1] += tempi;
            }
            wr=(wtemp=wr)*wpr-wi*wpi+wr;	// trigonometric recurrance
            wi=wi*wpr+wtemp*wpi+wi;
        }
        mmax=istep;
    }
}

/** From Numerical Recipes in C
	Calculates the Fourier transform of a set of n real-valued data points. Replaces this data 
	(which is stored in array data[1..n]) by the positive frequency half of its complex Fourier 
	transform. The real-valued first and last components of hte complex transform are returned as 
	elements data[1] and data[2], respectively. n must be a power of 2. This routine also calculates 
	the inverse transform of a complex data array if it is the transform of real data. (Result in 
	this case must be multiplied by 2/n.)
*/
template <class T>
void realft(T data[],int n, int isign)
{
    const T PI = 3.141592654;

    int i,i1,i2,i3,i4,np3;
    T c1=0.5,c2,h1r,h1i,h2r,h2i;
    T wr,wi,wpr,wpi,wtemp,theta;	// double precision for hte trigonometric recurrences
    theta=PI/(T) (n>>1);	// initialize the recurrence;
    if (isign == 1) {
        c2 = -0.5;
        four1<T>(data,n>>1,1);	// the forward transform is here.
    } else {
        c2=0.5;	// otherwise set up for an inverse transform.
        theta = -theta;
    }
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0+wpr;
    wi=wpi;
    np3=n+3;
    for (i=2;i<=(n>>2);i++) {	// case i=1 done separately below.
        i4=1+(i3=np3-(i2=1+(i1=i+i-1)));
        h1r=c1*(data[i1]+data[i3]);	// the two separate transforms are separated our of data.
        h1i=c1*(data[i2]-data[i4]);
        h2r = -c2*(data[i2]+data[i4]);
        h2i=c2*(data[i1]-data[i3]);
        data[i1]=h1r+wr*h2r-wi*h2i;	// Here they are recombined to form the true transform 
        data[i2]=h1i+wr*h2i+wi*h2r;	// of the original real data.
        data[i3]=h1r-wr*h2r+wi*h2i;
        data[i4] = -h1i+wr*h2i+wi*h2r;
        wr=(wtemp=wr)*wpr-wi*wpi+wr;	// the recurrence
        wi=wi*wpr+wtemp*wpi+wi;
    }
    if (isign == 1) {
        data[1] = (h1r=data[1])+data[2];	// Squeeze the first and last data together to get 
        data[2] = h1r-data[2];			// them all within the original array.
    } else {
        data[1]=c1*((h1r=data[1])+data[2]);
        data[2]=c1*(h1r-data[2]);
        four1<T>(data,n>>1,-1);		// This is the inverse transform for the case isign=-1.
    }
}



template<class T, class TimeType>
FWave1D<T,TimeType>& FWave1D<T,TimeType>::operator = (const FWave1D<T,TimeType>& A)
{
    (*this).copyTime(A);
    ref(A);
    return (*this);
}

/** Despite the innocuous appearance, this assignment operator does
	a Fourier transform
*/
template <class T,class TimeType>
FWave1D<T,TimeType>& FWave1D<T,TimeType>::operator = (const TWave1D<T,TimeType>& y)
{
    if(isPeriodic())
    {
        (*this) = FWave1D<T,TimeType>(nextPowerOfTwo(y.dim()));

        // copy time across different data types
        (*this).dt(y.dt());
        (*this).startt(y.startt());
        (*this).originalLength(y.dim());
        (*this).isPeriodic(true);
        
        // first load everything into the fourier array
        int i=0;
        for(; i<y.dim(); i++)
            ((Array1D<T>) (*this))[i] = ((const Array1D<T>) y)[i];
        for(; i<(*this).dim(); i++)
            ((Array1D<T>) (*this))[i] = 0.0;
    }
    else	// if we're doing an operation on a non-periodic function, double the length and
                // play the various games to make it appear periodic about the boundaries
    {
        (*this) = FWave1D<T,TimeType>(nextPowerOfTwo(y.dim()*2));

        // copy time across different data types
        (*this).dt(y.dt());
        (*this).startt(y.startt());
        (*this).originalLength(y.dim());
        (*this).isPeriodic(false);
        
        // first load everything into the fourier array
        for(int i=0; i<y.dim(); i++)
            ((Array1D<T>) (*this))[i] = y[i];
        for(int i=y.dim(); i<(y.dim()+y.dim()/2); i++)
            ((Array1D<T>) (*this))[i] = y[y.dim()-1]*2-y[y.dim()*2-i-2];
        for(int i=(y.dim()+y.dim()/2); i<((*this).dim()-y.dim()/2); i++)
            ((Array1D<T>) (*this))[i] = y[y.dim()/2];
        for(int i=((*this).dim()-y.dim()/2); i<(*this).dim(); i++)
            ((Array1D<T>) (*this))[i] = y[0]*2-y[(*this).dim()-i];
    }
    // do the transform
    realft<T>( (T*)((Array1D<T>) (*this))-1, (*this).dim(), 1);
    
    return (*this);
}

template <class T, class TimeType>
void FWave1D<T,TimeType>::initTime(int ol)
{
    // first do a check to make sure we've got a properly sized array
    if(this->dim() != nextPowerOfTwo(this->dim()))
        throw FOURIER_MUST_BE_POWER_OF_TWO;

    // if they make this mistake, just quietly correct it.
    if(ol > this->dim())
        ol = this->dim();
    
    original_length_=ol;
    is_periodic_ = true;
}

typedef FWave1D<double,double> FReal1D;

template <class T, class U, class TimeType, U (*func)(const T)>
FWave1D<U,TimeType> byMember(const FWave1D<T,TimeType> &x)
{
    FWave1D<U,TimeType> q(x.dim());
    q.copyTime(x);
    for(int i=0;i<x.dim();i++)
    {
        q[i] = (*func)(x[i]);    
    }
    return q;
}

template <class T, class TimeType>
TWave1D<T,TimeType> FWave1D<T,TimeType>::magsqt() const
{
	TWave1D<T,TimeType> twave(this->dim()/2+1);
	twave[0] = (*this)[0];
	twave[twave.dim()-1] = (*this)[1];
	for(int i=1; i<(twave.dim()-1); i++)
	{
		twave[i] = ((*this)[2*i])*((*this)[2*i]) + ((*this)[2*i+1])*((*this)[2*i+1]);
	}
	return twave;
}

template <class T, class TimeType>
TWave1D<T,TimeType> FWave1D<T,TimeType>::magsqtf() const
{
	TWave1D<T,TimeType> twave(this->dim()/2+1);
	FWave1D<T,TimeType> freqs = this->f();
	twave[0] = freqs[0];
	twave[twave.dim()-1] = freqs[1];
	for(int i=1; i<(twave.dim()-1); i++)
	{
		twave[i] = freqs[2*i];
	}
	return twave;
}


#endif

