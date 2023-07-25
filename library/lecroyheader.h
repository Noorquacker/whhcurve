/** \file
	Defines a useful subset of the Lecroy binary wave format.
*/
#ifndef _LECROY_FILE_FORMAT_H
#define _LECROY_FILE_FORMAT_H

namespace lecroyheader
{

typedef int8_t I8;
typedef int8_t BYTE;

typedef int16_t I16;
typedef int16_t WORD;

typedef int32_t I32;
typedef int32_t INT;

typedef int64_t I64;
typedef int64_t LONG;

typedef float F32;
typedef float FLOAT;

typedef double F64;
typedef double DOUBLE;

// packing necessary because 64-bit operating systems don't align on 32 bits
#pragma pack(1)
typedef struct
{
    DOUBLE seconds;
    BYTE minutes;
    BYTE hours;
    BYTE days;
    BYTE months;
    WORD year;
    WORD unused;
} TIMESTAMP;

const WORD COMM_TYPE_BYTE = 0;
const WORD COMM_TYPE_WORD = 1;

const WORD COMM_ORDER_HIFIRST = 0;
const WORD COMM_ORDER_LOFIRST = 1;

// packing necessary because 64-bit operating systems don't align on 32 bits
#pragma pack(1)
typedef struct
{
    char DESCRIPTOR_NAME[16];
    char TEMPLATE_NAME[16];
    WORD COMM_TYPE;
    WORD COMM_ORDER;
    
    I32 WAVE_DESCRIPTOR;   // byte 36
    I32 USER_TEXT;
    I32 RES_DESC1;

    I32 TRIGTIME_ARRAY;    // byte 48
    I32 RIS_TIME_ARRAY;
    I32 RES_ARRAY1;
    I32 WAVE_ARRAY_1;
    I32 WAVE_ARRAY_2;
    I32 RES_ARRAY2;
    I32 RES_ARRAY3;
    
    char INSTRUMENT_NAME[16];   // byte 76
    I32 INSTRUMENT_NUMBER;
    char TRACE_LABEL[16];
    WORD RESERVED1;
    WORD RESERVED2;
    
    I32 WAVE_ARRAY_COUNT;  // byte 116
    I32 PNTS_PER_SCREEN;
    I32 FIRST_VALID_PNT;
    I32 LAST_VALID_PNT;
    I32 FIRST_POINT;
    I32 SPARSING_FACTOR;
    I32 SEGMENT_INDEX;
    I32 SUBARRAY_COUNT;
    I32 SWEEPS_PER_ACQ;
    WORD POINTS_PER_PAIR;
    WORD PAIR_OFFSET;
    
    FLOAT VERTICAL_GAIN;    // byte 156
    FLOAT VERTICAL_OFFSET;  // float value = VERTICAL_GAIN*data-VERTICAL_OFFSET;

    FLOAT MAX_VALUE;
    FLOAT MIN_VALUE;
    WORD NOMINAL_BITS;
    WORD NOM_SUBARRAY_COUNT;

    FLOAT HORIZ_INTERVAL;   // byte 176
    DOUBLE HORIZ_OFFSET;

    DOUBLE PIXEL_OFFSET;
    char VERTUNIT[48];
    char HORUNIT[48];

    FLOAT HORIZ_UNCERTAINTY;    // byte 292
    TIMESTAMP TRIGGER_TIME;
    FLOAT ACQ_DURATION;
    
    WORD RECORD_TYPE;   // byte 316
    WORD PROCESSING_DONE;
    
    WORD RESERVED5;
    WORD RIS_SWEEPS;
    WORD TIMEBASE;
    WORD VERT_COUPLING;
    FLOAT PROBE_ATT;
    WORD FIXED_VERT_GAIN;
    
    WORD BANDWIDTH_LIMIT;   // byte 334
    FLOAT VERTICAL_VERNIER;
    FLOAT ACQ_VERT_OFFSET;
    WORD WAVE_SOURCE;
} WAVEDESC_BLOCK;
    
}

#endif

