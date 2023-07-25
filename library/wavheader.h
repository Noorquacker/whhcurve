/** \file
	Defines a useful subset of the WAV audio format.
*/
#ifndef _WAV_FILE_FORMAT_H
#define _WAV_FILE_FORMAT_H

namespace wavheader
{

typedef signed char I8;
typedef signed char BYTE;

typedef signed short int I16;
typedef signed short int WORD;

typedef signed long int I32;
typedef signed long int LONG;

typedef float F32;
typedef float FLOAT;

typedef double F64;
typedef double DOUBLE;


typedef struct
{
    // riff header
    char ChunkID[4];    // "RIFF"
    I32 ChunkSize;      // filesize - 8
    char Format[4];     // "WAVE"
    // fmt subchunk
    char Subchunk1ID[4] // "fmt "
    I32 Subchunk1Size;  // 16 for PCM
    I16 AudioFormat;    // 1 for PCM
    I16 NumChannels;    // 1=mono, 2=stereo
    I32 SampleRate;     // 8000, 44100, etc
    I32 ByteRate;       // SampleRate * NumChannels * BitsPerSample/8
    I16 BlockAlign;     // NumChannels * BitsPerSample/8
    I16 BitsPerSample;  // 8 or 16
    // data subchunk
    char Subchunk2ID[4] // "data"
    I32 Subchunk2Size;  // NumSamples * NumChannels * BitsPerSample/8
    // then data goes here for Subchunk2Size bytes
} WAV_HEADER;
    
}

#endif

