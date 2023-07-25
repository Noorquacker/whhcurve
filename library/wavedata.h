/** \file
	TO-DO/BUGS LIST:
	in open(fn), test for a bad file opening and chuck an error out
	at the end of readText, iterate over that getline() until we make sure we 
	 get the whole thing out
	needs a copy constructor
	need to subclass this to read wave names
	needs assignment/constructor from an array
	needs ability to WRITE files
	perhaps need to break this apart into a general superclass, binary, and text
	bug: need bounds checking in getDataColumn
*/
#ifndef _WAVEDATA_H
#define _WAVEDATA_H

using namespace std;

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "tnt/tnt.h"
using namespace TNT;

extern const char* NO_INPUT_STREAM;
extern const char* NO_OUTPUT_STREAM;


template <class T>
class WaveData
{
    public:
        WaveData() { create(); }
        WaveData(const char* fn) { create(); open(fn); }
        WaveData(istream &xpin) { create(); open(xpin); }
        WaveData(int cols, int rows);
        virtual ~WaveData();
        
        void open(const char* fn) { openr(fn); }    // for compatibility
        void openr(const char* fn);
        void openw(const char* fn);
        void open(istream &xpin) { pin = &xpin; should_delete_pin = false; }
        void open(ostream &xpout) { pout = &xpout; should_delete_pout = false; }
        void close();   // closes ALL streams

        // Subclasses should overload this function to read their data
        //virtual void read(int maxlength = -1);
        
        // Note that operator[] is used to get a column. We are a stream
        // reader, not an array.
        // 
        // Also notice that we return NON-const arrays, and Array1D uses
        // copy-by-reference. You may modify them in-place, and save to
        // another file.

		// 12/1/09: this business of returning a reference may be a problem. taking
		// it out for a bit
//        Array2D<T>& getData2D() { return data2d; }
        Array2D<T> getData2D() { return data2d; }
        // note the difference between the above direct (by-reference) function
        // and the below by-value ones. if you want to assign an Array to WaveData,
        // you must do so with getData (as of right now)
        Array1D<T> getData() { return data1d; }
        Array1D<T> getDataColumn(unsigned c) { return Array1D<T>(data2d.dim2(),data2d[c]); }
        Array1D<T> operator[](unsigned c) { return getDataColumn(c); }

    protected:
        void create();
    
        istream* pin;
        ostream* pout;
        bool should_delete_pin;
        bool should_delete_pout;

        // first (or single) column data1d pointer
        Array1D<T> data1d;
        
        // multi-column data1d pointer array
        Array2D<T> data2d;
        
        double data_dt;
};

template <class T>
WaveData<T>::WaveData(int cols, int rows)
{
    create();
    (this->data2d) = Array2D<T>(cols,rows);
    (this->data1d) = Array1D<T>(rows,(this->data2d)[0]);

}

template <class T>
void WaveData<T>::create()
{
    pin = 0;
    pout = 0;
}

template <class T>
WaveData<T>::~WaveData()
{
    if(should_delete_pin)
        delete pin;
    if(should_delete_pout)
        delete pout;
}

template <class T>
void WaveData<T>::openr(const char* fn)
{
    pin = new ifstream(fn);
    should_delete_pin = true;
}

template <class T>
void WaveData<T>::openw(const char* fn)
{
    pout = new ofstream(fn);
    should_delete_pout = true;
}

template <class T>
void WaveData<T>::close()
{
    if(should_delete_pin)
    {
        ((ifstream*)pin)->close();
        delete pin;
    }
    if(should_delete_pout)
    {
        ((ofstream*)pout)->close();
        delete pout;
    }
    pin = pout = 0;
    should_delete_pin = should_delete_pout = false;
}

#endif // _WAVEDATA_H
