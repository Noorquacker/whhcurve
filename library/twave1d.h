#ifndef _TWAVE1D_H
#define _TWAVE1D_H

#include "tnt/tnt.h"
#include "wave1d.h"

using namespace TNT;

template <class T, class TimeType>
class TWave1D : public Wave1D<T, TimeType>
{
public:
    TWave1D() : Wave1D<T, TimeType>::Wave1D() {}
    TWave1D(int n) : Wave1D<T, TimeType>::Wave1D(n) {}
    TWave1D(int n, const T &a) : Wave1D<T, TimeType>::Wave1D(n, a) {}
    TWave1D(int n, T *a) : Wave1D<T, TimeType>::Wave1D(n, a) {}
    TWave1D(const Array1D<T> &A) : Wave1D<T, TimeType>::Wave1D(A) {}
    TWave1D(const Wave1D<T, TimeType> &A) { throw CANNOT_CONSTRUCT_FROM_WAVE1D; }
    TWave1D(const TWave1D<T, TimeType> &A) { (*this) = A; }
    TWave1D(const FWave1D<T, TimeType> &B) { (*this) = B; }

    inline TimeType mint() const { return (*this).startt(); }
    inline TimeType maxt() const { return (*this).startt() + (*this).dim() * (*this).dt(); }
    inline TimeType mint(const TimeType &min);
    inline TimeType maxt(const TimeType &max);

    TWave1D<TimeType, TimeType> t() const;
    TimeType t(int i) const { return (*this).startt() + i * (*this).dt(); }
    int index(const TimeType t) const { return (int)((t - (this->startt())) / (this->dt()) + 0.5); }

    // duplicate the assignment functions with proper return types
    TWave1D<T, TimeType> &ref(const TWave1D<T, TimeType> &A)
    {
        Wave1D<T, TimeType>::ref(A);
        return (*this);
    }
    TWave1D<T, TimeType> copy() const;
    TWave1D<T, TimeType> &inject(const TWave1D<T, TimeType> &A)
    {
        Wave1D<T, TimeType>::ref(A);
        return (*this);
    }
    TWave1D<T, TimeType> subwave(const int x0, const int length) const;

    TWave1D<T, TimeType> &operator=(const TWave1D<T, TimeType> &A);
    TWave1D<T, TimeType> &operator=(const FWave1D<T, TimeType> &f);
};

template <class T, class TimeType>
TimeType TWave1D<T, TimeType>::mint(const TimeType &min)
{
    TimeType max = maxt();
    Wave1D<T, TimeType>::startt(min);
    Wave1D<T, TimeType>::dt((max - min) / Array1D<T>::dim());
    return min;
}

template <class T, class TimeType>
TimeType TWave1D<T, TimeType>::maxt(const TimeType &max)
{
    Wave1D<T, TimeType>::dt((max - mint()) / Array1D<T>::dim());
    return max;
}

/** get a wave whose values are the time values at each point, matching this one */
template <class T, class TimeType>
TWave1D<TimeType, TimeType> TWave1D<T, TimeType>::t() const
{
    TWave1D<TimeType, TimeType> times(Array1D<T>::dim());
    times.copyTime(*this);
    for (int i = 0; i < (*this).dim(); i++)
        times[i] = (*this).startt() + i * (*this).dt();
    return times;
}

template <class T, class TimeType>
TWave1D<T, TimeType> TWave1D<T, TimeType>::copy() const
{
    TWave1D<T, TimeType> A(Array1D<T>::copy());
    A.copyTime((*this));

    return A;
}

template <class T, class TimeType>
TWave1D<T, TimeType> TWave1D<T, TimeType>::subwave(const int x0, const int length) const
{
    // the type cast ((Array1D<T>) *this) must happen b/c Array1D::subarray() isn't defined
    // as const (although it is, sort of)
    TWave1D<T, TimeType> a(((Array1D<T>)*this).subarray(x0, length + x0 - 1));
    a.dt(this->dt());
    a.startt(this->startt() + (this->dt()) * x0);
    return a;
}

template <class T, class TimeType>
TWave1D<T, TimeType> &TWave1D<T, TimeType>::operator=(const TWave1D<T, TimeType> &A)
{
    (*this).copyTime(A);
    ref(A);
    return (*this);
}

/** convert an FWave1D to our TWave1D by performing an inverse Fourier transform */
template <class T, class TimeType>
TWave1D<T, TimeType> &TWave1D<T, TimeType>::operator=(const FWave1D<T, TimeType> &f)
{
    (*this) = TWave1D<T, TimeType>(f.originalLength());
    FWave1D<T, TimeType> tf = f.copy(); // work from a copy

    // copy time across different data types
    (*this).dt(tf.dt());
    (*this).startt(tf.startt());

    // do the inverse transform
    realft<T>((T *)((Array1D<T>)tf) - 1, tf.dim(), -1);

    // move the data over into the wave array
    for (int i = 0; i < (*this).dim(); i++)
        ((Array1D<T>)(*this))[i] = (2.0 / tf.dim()) * ((Array1D<T>)tf)[i];

    return (*this);
}

typedef TWave1D<double, double> TReal1D;

template <class T, class U, class TimeType, U (*func)(const T)>
TWave1D<U, TimeType> byMember(const TWave1D<T, TimeType> &x)
{
    TWave1D<U, TimeType> q(x.dim());
    q.copyTime(x);
    for (int i = 0; i < x.dim(); i++)
    {
        q[i] = (*func)(x[i]);
    }
    return q;
}

#endif
