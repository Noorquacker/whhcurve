#ifndef _PEAKDATA_H
#define _PEAKDATA_H

#include "../tntnr.h"

template <class T>
class PeakData
{
    public:
    PeakData() { init(0,0,0); }
    PeakData(const T& t, const T& a, int d) { init(t,a,d); }
    PeakData(const PeakData<T>& other) { init(other.time(),other.amp(),other.direction()); }

    const PeakData<T>& operator=(const PeakData<T>& other) { init(other.time(),other.amp(),other.direction()); return (*this); }

    /**
        calculates frequency between this peak and peak q
    */
    T freq(const PeakData& q)
    {
        // peaks are both in same direction (full wavelength)
        if(this->data_direction * q.data_direction > 0)
            return abs<T>(1.0/(this->data_time - q.data_time));
        // peaks are in opposite directions (half wavelength)
        else
            return abs<T>(0.5/(this->data_time - q.data_time));            
    }
    
    T time(const T& t) { return data_time = t; }
    T amp(const T& a) { return data_amp = a; }
    int direction(const int d) { return data_direction = d; }
    T time() const { return data_time; }
    T amp() const { return data_amp; }
    int direction() const { return data_direction; }
    
    protected:
    T data_time;
    T data_amp;
    int data_direction;
    
    void init(const T& t, const T& a, int d)
    {
        data_time = t;
        data_amp = a;
        data_direction = d;
    }         
        
};

#endif // _PEAKDATA_H
