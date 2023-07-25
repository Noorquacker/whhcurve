/** \file
	TO-DO LIST:
	Index of the trigger
	Convert an integer value to a voltage
	Store/spit out the contents of any text block present
*/
#ifndef _LECROY_WAVE_DATA_H
#define _LECROY_WAVE_DATA_H

#include <cstring>

#include "lecroyheader.h"
#include "binarywavedata.h"
#include "endian.h"

using namespace lecroyheader;

extern const char* NOT_A_LECROY_TRACE;

/** used to copy strings around, but keep unused memory tidy and null. */
void strcpy_pad_null(char* dest, const char* src, int destsize);

/** Subclass of BinaryWaveData to read a Lecroy file */
template <class T>
class LecroyWaveData : public BinaryWaveData<T>
{
    public:
        LecroyWaveData() : BinaryWaveData<T>::BinaryWaveData() { }
        LecroyWaveData(const char* fn) : BinaryWaveData<T>::BinaryWaveData(fn) { }
        LecroyWaveData(istream &xpin) : BinaryWaveData<T>::BinaryWaveData(xpin) { }
        LecroyWaveData(unsigned long pts) : BinaryWaveData<T>::BinaryWaveData(pts) { }
        
        void read() { readHeader(); BinaryWaveData<T>::read(); }
        void readHeader();
        void newHeader();
        void write();
        
        lecroyheader::WAVEDESC_BLOCK header;
    
    protected:
        void headerSwitchEndian();
    
};

template <class T>
void LecroyWaveData<T>::readHeader()
{
    char buffer[128];

//    cout << "I8: " << sizeof(I8) << endl;
//    cout << "I16: " << sizeof(I16) << endl;
//    cout << "I32: " << sizeof(I32) << endl;
//    cout << "I64: " << sizeof(I64) << endl;
//    cout << "F32: " << sizeof(F32) << endl;
//    cout << "F64: " << sizeof(F64) << endl;
	
    (*WaveData<T>::pin).seekg(0,ios::beg);
    (*WaveData<T>::pin).read(buffer,127);
    buffer[127] = '\0';
    
    char* pstr = strstr(buffer,"WAVEDESC");
    if(pstr==0)
        throw NOT_A_LECROY_TRACE;
	
//	cout << buffer << endl;
//	cout.flush();
    
    unsigned startpos = pstr - buffer;
//	cout << "startpos: " << startpos << endl;
//	cout.flush();
//    (*WaveData<T>::pin).seekg(startpos,ios::beg);
    
	for(int i=startpos; i<127; i++)
	{
		((char*)&header)[i-startpos] = buffer[i];
	}
	
    (*WaveData<T>::pin).read((((char*)&header)+127-startpos), (346-127+startpos));
//	for(int i=0; i<346; i++)
//	{
//		cout << "header[" << i << "]: " << hex << (int)(unsigned char)((char*)&header)[i] << dec << endl;
//	}
//	cout << header.INSTRUMENT_NAME << "dt = " << header.HORIZ_INTERVAL << ", start = " << header.HORIZ_OFFSET << endl;

    // Byte ordering stuff:

    // Override WaveData's default processing of byte orders
    BinaryWaveData<T>::ASSUME_LOFIRST_BYTEORDER = false;
    // Subtlety in the process of reading an endian indicator using a data type
    // whose bytes might need to be switched. To be simple and robust, since
    // high-first is indicated by a zero, treat anything else as low-first.
    if(machine_endian() == MACHINE_LITTLE_ENDIAN && header.COMM_ORDER == COMM_ORDER_HIFIRST
        || machine_endian() == MACHINE_BIG_ENDIAN && header.COMM_ORDER != COMM_ORDER_HIFIRST)
    {
        headerSwitchEndian();
        BinaryWaveData<T>::REVERSE_BYTES_ANYWAY = true;
    }
    else
        BinaryWaveData<T>::REVERSE_BYTES_ANYWAY = false;

//	cout << "before seek to data: " << (*WaveData<T>::pin).tellg() << endl;
    // seek forward to beginning of binary data
//    (*WaveData<T>::pin).seekg(startpos+header.WAVE_DESCRIPTOR+header.USER_TEXT+header.RES_DESC1,ios::beg);

    // and set the wave's dt
    BinaryWaveData<T>::dt(header.HORIZ_INTERVAL);
}

template <class T>
void LecroyWaveData<T>::write()
{
    if(this->pout==0)
        throw NO_OUTPUT_STREAM;    

    if(machine_endian() == MACHINE_BIG_ENDIAN)
        headerSwitchEndian();

    this->pout->seekp(0);
    this->pout->write((const char*) &header, sizeof(header));
    BinaryWaveData<T>::write();
}

/** the truth is that I really don't want to type all those assignments */
#define HEADER_SET(val,foo) header.foo = val
template <class T>
void LecroyWaveData<T>::newHeader()
{
    HEADER_SET(COMM_TYPE_BYTE,COMM_TYPE);
    HEADER_SET(COMM_ORDER_LOFIRST,COMM_ORDER);
    HEADER_SET(0,WAVE_DESCRIPTOR);
    HEADER_SET(0,USER_TEXT);
    HEADER_SET(0,RES_DESC1);
    HEADER_SET(0,TRIGTIME_ARRAY);
    HEADER_SET(0,RIS_TIME_ARRAY);
    HEADER_SET(0,RES_ARRAY1);
    HEADER_SET((this->data1d).dim(),WAVE_ARRAY_1);
    HEADER_SET(0,WAVE_ARRAY_2);
    HEADER_SET(0,RES_ARRAY2);
    HEADER_SET(0,RES_ARRAY3);
    HEADER_SET(0,INSTRUMENT_NUMBER);
    HEADER_SET(0,RESERVED1);
    HEADER_SET(0,RESERVED2);
    HEADER_SET((this->data1d).dim(),WAVE_ARRAY_COUNT);
    HEADER_SET(0,PNTS_PER_SCREEN);
    HEADER_SET(0,FIRST_VALID_PNT);
    HEADER_SET((this->data1d).dim()-1,LAST_VALID_PNT);
    HEADER_SET(0,FIRST_POINT);
    HEADER_SET(0,SPARSING_FACTOR);
    HEADER_SET(0,SEGMENT_INDEX);
    HEADER_SET(0,SUBARRAY_COUNT);
    HEADER_SET(0,SWEEPS_PER_ACQ);
    HEADER_SET(0,POINTS_PER_PAIR);
    HEADER_SET(0,PAIR_OFFSET);
    HEADER_SET(1.0,VERTICAL_GAIN);
    HEADER_SET(0,VERTICAL_OFFSET);
    HEADER_SET(0,MAX_VALUE);
    HEADER_SET(0,MIN_VALUE);
    HEADER_SET(8,NOMINAL_BITS);	// note 8, not zero
    HEADER_SET(1,NOM_SUBARRAY_COUNT);
    HEADER_SET(this->dt(),HORIZ_INTERVAL);
    HEADER_SET(0,HORIZ_OFFSET);
    HEADER_SET(0,PIXEL_OFFSET);
    HEADER_SET(0,HORIZ_UNCERTAINTY);
    HEADER_SET(0,TRIGGER_TIME.seconds);
    HEADER_SET(0,TRIGGER_TIME.minutes);
    HEADER_SET(0,TRIGGER_TIME.hours);
    HEADER_SET(0,TRIGGER_TIME.days);
    HEADER_SET(0,TRIGGER_TIME.months);
    HEADER_SET(0,TRIGGER_TIME.year);
    HEADER_SET(0,TRIGGER_TIME.unused);
    HEADER_SET((this->data1d).dim()*(this->dt()),ACQ_DURATION);
    HEADER_SET(0,RECORD_TYPE);
    HEADER_SET(0,PROCESSING_DONE);
    HEADER_SET(0,RESERVED5);
    HEADER_SET(0,RIS_SWEEPS);
    HEADER_SET(0,TIMEBASE);
    HEADER_SET(0,VERT_COUPLING);
    HEADER_SET(0,PROBE_ATT);
    HEADER_SET(0,FIXED_VERT_GAIN);
    HEADER_SET(0,BANDWIDTH_LIMIT);
    HEADER_SET(0,VERTICAL_VERNIER);
    HEADER_SET(0,ACQ_VERT_OFFSET);
    HEADER_SET(0,WAVE_SOURCE);

    strcpy_pad_null(header.DESCRIPTOR_NAME,"WAVEDESC",16);
    strcpy_pad_null(header.TEMPLATE_NAME,"LECROY_2_3",16);
    strcpy_pad_null(header.INSTRUMENT_NAME,"C++ app, so ha!",16);
    strcpy_pad_null(header.TRACE_LABEL,"",16);
    strcpy_pad_null(header.VERTUNIT,"",48);
    strcpy_pad_null(header.HORUNIT,"",48);
}

/** the truth is that I really don't want to type all those assignments */
#define HEADER_REVERSE_IN_PLACE(type,foo) header.foo = reverse_bytes<type>(header.foo)

template <class T>
void LecroyWaveData<T>::headerSwitchEndian()
{
    HEADER_REVERSE_IN_PLACE(WORD,COMM_TYPE);
    HEADER_REVERSE_IN_PLACE(WORD,COMM_ORDER);
    HEADER_REVERSE_IN_PLACE(LONG,WAVE_DESCRIPTOR);
    HEADER_REVERSE_IN_PLACE(LONG,USER_TEXT);
    HEADER_REVERSE_IN_PLACE(LONG,RES_DESC1);
    HEADER_REVERSE_IN_PLACE(LONG,TRIGTIME_ARRAY);
    HEADER_REVERSE_IN_PLACE(LONG,RIS_TIME_ARRAY);
    HEADER_REVERSE_IN_PLACE(LONG,RES_ARRAY1);
    HEADER_REVERSE_IN_PLACE(LONG,WAVE_ARRAY_1);
    HEADER_REVERSE_IN_PLACE(LONG,WAVE_ARRAY_2);
    HEADER_REVERSE_IN_PLACE(LONG,RES_ARRAY2);
    HEADER_REVERSE_IN_PLACE(LONG,RES_ARRAY3);
    HEADER_REVERSE_IN_PLACE(LONG,INSTRUMENT_NUMBER);
    HEADER_REVERSE_IN_PLACE(WORD,RESERVED1);
    HEADER_REVERSE_IN_PLACE(WORD,RESERVED2);
    HEADER_REVERSE_IN_PLACE(LONG,WAVE_ARRAY_COUNT);
    HEADER_REVERSE_IN_PLACE(LONG,PNTS_PER_SCREEN);
    HEADER_REVERSE_IN_PLACE(LONG,FIRST_VALID_PNT);
    HEADER_REVERSE_IN_PLACE(LONG,LAST_VALID_PNT);
    HEADER_REVERSE_IN_PLACE(LONG,FIRST_POINT);
    HEADER_REVERSE_IN_PLACE(LONG,SPARSING_FACTOR);
    HEADER_REVERSE_IN_PLACE(LONG,SEGMENT_INDEX);
    HEADER_REVERSE_IN_PLACE(LONG,SUBARRAY_COUNT);
    HEADER_REVERSE_IN_PLACE(LONG,SWEEPS_PER_ACQ);
    HEADER_REVERSE_IN_PLACE(WORD,POINTS_PER_PAIR);
    HEADER_REVERSE_IN_PLACE(WORD,PAIR_OFFSET);
    HEADER_REVERSE_IN_PLACE(FLOAT,VERTICAL_GAIN);
    HEADER_REVERSE_IN_PLACE(FLOAT,VERTICAL_OFFSET);
    HEADER_REVERSE_IN_PLACE(FLOAT,MAX_VALUE);
    HEADER_REVERSE_IN_PLACE(FLOAT,MIN_VALUE);
    HEADER_REVERSE_IN_PLACE(WORD,NOMINAL_BITS);
    HEADER_REVERSE_IN_PLACE(WORD,NOM_SUBARRAY_COUNT);
    HEADER_REVERSE_IN_PLACE(FLOAT,HORIZ_INTERVAL);
    HEADER_REVERSE_IN_PLACE(DOUBLE,HORIZ_OFFSET);
    HEADER_REVERSE_IN_PLACE(DOUBLE,PIXEL_OFFSET);
    HEADER_REVERSE_IN_PLACE(FLOAT,HORIZ_UNCERTAINTY);
    HEADER_REVERSE_IN_PLACE(DOUBLE,TRIGGER_TIME.seconds);
    HEADER_REVERSE_IN_PLACE(BYTE,TRIGGER_TIME.minutes);
    HEADER_REVERSE_IN_PLACE(BYTE,TRIGGER_TIME.hours);
    HEADER_REVERSE_IN_PLACE(BYTE,TRIGGER_TIME.days);
    HEADER_REVERSE_IN_PLACE(BYTE,TRIGGER_TIME.months);
    HEADER_REVERSE_IN_PLACE(WORD,TRIGGER_TIME.year);
    HEADER_REVERSE_IN_PLACE(WORD,TRIGGER_TIME.unused);
    HEADER_REVERSE_IN_PLACE(FLOAT,ACQ_DURATION);
    HEADER_REVERSE_IN_PLACE(WORD,RECORD_TYPE);
    HEADER_REVERSE_IN_PLACE(WORD,PROCESSING_DONE);
    HEADER_REVERSE_IN_PLACE(WORD,RESERVED5);
    HEADER_REVERSE_IN_PLACE(WORD,RIS_SWEEPS);
    HEADER_REVERSE_IN_PLACE(WORD,TIMEBASE);
    HEADER_REVERSE_IN_PLACE(WORD,VERT_COUPLING);
    HEADER_REVERSE_IN_PLACE(FLOAT,PROBE_ATT);
    HEADER_REVERSE_IN_PLACE(WORD,FIXED_VERT_GAIN);
    HEADER_REVERSE_IN_PLACE(WORD,BANDWIDTH_LIMIT);
    HEADER_REVERSE_IN_PLACE(FLOAT,VERTICAL_VERNIER);
    HEADER_REVERSE_IN_PLACE(FLOAT,ACQ_VERT_OFFSET);
    HEADER_REVERSE_IN_PLACE(WORD,WAVE_SOURCE);
}

#endif
