/** \file
	Load binary data from a file
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

#ifndef _BINARYWAVEDATA_H
#define _BINARYWAVEDATA_H

using namespace std;

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

#include "wavedata.h"
#include "wave1d.h"
#include "endian.h"

#include "tnt/tnt.h"
using namespace TNT;

//extern const char* TEXT_DATA_DOES_NOT_MATCH_DESCRIPTION;	// unused for the moment
extern const char* TIME_OUT_OF_RANGE;
extern const char* INDEX_OUT_OF_RANGE;
extern const char* DT_MUST_BE_POSITIVE;

/** used to copy strings around, but keep unused memory tidy and null. */
void strcpy_pad_null(char* dest, const char* src, unsigned long destsize);


/**	Binary load loads all data1d into memory verbatim, then places (T*)data1d at
	the head of the array
*/
template <class T>
class BinaryWaveData : public WaveData<T>
{
    public:
        BinaryWaveData() : WaveData<T>::WaveData() { create(); }
        BinaryWaveData(const char* fn) : WaveData<T>::WaveData(fn) { create(); }
        BinaryWaveData(istream &xpin) : WaveData<T>::WaveData(xpin) { create(); }
        BinaryWaveData(int pts) : WaveData<T>::WaveData(1,pts) { create(); }
        
        void read(int maxlength = -1);
        void write();

        int indexAtTime(double t);
        double timeAtIndex(int i);

        double dt(double zdt);
        double dt() { return data_dt; }
        
        /** Assume that the file was generated according to "little endian"
        	byte ordering (probably a safe bet). If false, no checking or
        	processing is done. Default true, set in constructor.
        */
        bool ASSUME_LOFIRST_BYTEORDER;
        /**	For a little more control over the process, subclasses may find
        	it useful to set this if they need the bytes of input types
        	reversed for whatever reason.
        */
        bool REVERSE_BYTES_ANYWAY;
        
    protected:
        void create();

        double data_dt;
};

template <class T>
void BinaryWaveData<T>::create()
{
    ASSUME_LOFIRST_BYTEORDER = true;
    REVERSE_BYTES_ANYWAY = false;
    
    // default time spacing
    data_dt = 1;
}

template <class T>
void BinaryWaveData<T>::read(int maxlength)
{
    int start = (*(this->pin)).tellg();
    int length;
//	cout << "binarywavedata start: " << start << endl;

    // get length of file:
    (*(this->pin)).seekg(0, ios::end);
    length = (*(this->pin)).tellg();    // this line seperated from the next b/c ISO C++
//	cout << "binarywavedata end: " << length << endl;
    length -= start;            // complains about an operator-() type conversion
    (*(this->pin)).seekg(start, ios::beg);

    if(maxlength != -1 && maxlength < length)
        length = maxlength;
//cout << "binarywavedata length: " << length << endl;
    (this->data2d) = Array2D<T>(1,length);
    (this->data1d) = Array1D<T>(length,(this->data2d)[0]);
    
    // read data1d as a block:
    (*(this->pin)).read((char*) &((this->data1d)[0]), length*sizeof(T));
	
    // take care of byte-ordering if this needs to be done
    if((ASSUME_LOFIRST_BYTEORDER && machine_endian() == MACHINE_BIG_ENDIAN)
        || REVERSE_BYTES_ANYWAY)
    {
        for(int i=0; i<(this->data1d.dim)(); i++)
        {
            (this->data1d)[i] = reverse_bytes<T>((this->data1d)[i]);
        }
    }

}

template <class T>
void BinaryWaveData<T>::write()
{
    if(this->pout==0)
        throw NO_OUTPUT_STREAM;    

    Wave1D<T,T> w = this->data1d;
    
    this->pout->write((const char*) w.data(), w.dim());
}

template <class T>
int BinaryWaveData<T>::indexAtTime(double t)
{
    int result = (int)(t/data_dt+0.5);
    if(result > ((this->data1d).dim()-1) || result < 0)
        throw TIME_OUT_OF_RANGE;
    return result;
}

template <class T>
double BinaryWaveData<T>::timeAtIndex(int i)
{
    if(i > ((this->data1d).dim()-1))
        throw INDEX_OUT_OF_RANGE;

    return i*data_dt;
}

template <class T>
double BinaryWaveData<T>::dt(double zdt)
{
    if(zdt <= 0)
        throw DT_MUST_BE_POSITIVE;
        
    return data_dt = zdt;
}

#endif // _WAVEDATA_H
