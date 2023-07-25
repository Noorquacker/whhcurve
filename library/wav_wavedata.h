/** read a wav file
*/
#ifndef _WAV_WAVE_DATA_H
#define _WAV_WAVE_DATA_H

#include "wavheader.h"
#include "binarywavedata.h"
#include "endian.h"

using namespace Wavheader;

extern const char* NOT_A_Wav_TRACE;

/** used to copy strings around, but keep unused memory tidy and null. */
void strcpy_pad_null(char* dest, const char* src, int destsize);

/** Subclass of BinaryWaveData to read a Wav file */
template <class T>
class WavWaveData : public WavWaveData<T>
{
    public:
        WavWaveData() : BinaryWaveData<T>::BinaryWaveData() { }
        WavWaveData(const char* fn) : BinaryWaveData<T>::BinaryWaveData(fn) { }
        WavWaveData(istream &xpin) : BinaryWaveData<T>::BinaryWaveData(xpin) { }
        WavWaveData(int pts) : BinaryWaveData<T>::BinaryWaveData(pts) { }
        
        void read() { readHeader(); BinaryWaveData<T>::read(); }
        void readHeader();
        void newHeader();
        void write();
        
        Wavheader::WAVEDESC_BLOCK header;
    
    protected:
        void headerSwitchEndian();
    
};

template <class T>
void WavWaveData<T>::readHeader()
{
    char buffer[128];
    
    (*WaveData<T>::pin).seekg(0,ios::beg);
    (*WaveData<T>::pin).read(buffer,127);
    buffer[127] = '\0';
    
    char* pstr = strstr(buffer,"RIFF");
    if(pstr==0)
        throw NOT_A_WAV_FILE;
    
    unsigned startpos = pstr - buffer;
    (*WaveData<T>::pin).seekg(startpos,ios::beg);
    
    (*WaveData<T>::pin).read((char*) &header, 346);

    // Byte ordering stuff:

    if(machine_endian()==MACHINE_BIG_ENDIAN)
    {
        headerSwitchEndian();
    }

    // seek forward to beginning of binary data
    (*WaveData<T>::pin).seekg(44);

    // and set the wave's dt
    BinaryWaveData<T>::dt(((T)1)/header.SampleRate);
}

template <class T>
void WavWaveData<T>::write()
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
void WavWaveData<T>::newHeader()
{
    HEADER_SET(LONG,ChunkSize);
    HEADER_SET(LONG,Subchunk1Size);
    HEADER_SET(WORD,AudioFormat);
    HEADER_SET(WORD,NumChannels);
    HEADER_SET(LONG,SampleRate);
    HEADER_SET(LONG,ByteRate);
    HEADER_SET(WORD,BlockAlign);
    HEADER_SET(WORD,BitsPerSample);
    HEADER_SET(LONG,Subchunk2Size);

    strcpy_pad_null(header.ChunkID,"RIFF",4);
    strcpy_pad_null(header.Format,"WAVE",4);
    strcpy_pad_null(header.Subchunk1ID,"fmt ",4);
    strcpy_pad_null(header.Subchunk2ID,"data",4);
}

/** the truth is that I really don't want to type all those assignments */
#define HEADER_REVERSE_IN_PLACE(type,foo) header.foo = reverse_bytes<type>(header.foo)

template <class T>
void WavWaveData<T>::headerSwitchEndian()
{
    HEADER_REVERSE_IN_PLACE(LONG,ChunkSize);
    HEADER_REVERSE_IN_PLACE(LONG,Subchunk1Size);
    HEADER_REVERSE_IN_PLACE(WORD,AudioFormat);
    HEADER_REVERSE_IN_PLACE(WORD,NumChannels);
    HEADER_REVERSE_IN_PLACE(LONG,SampleRate);
    HEADER_REVERSE_IN_PLACE(LONG,ByteRate);
    HEADER_REVERSE_IN_PLACE(WORD,BlockAlign);
    HEADER_REVERSE_IN_PLACE(WORD,BitsPerSample);
    HEADER_REVERSE_IN_PLACE(LONG,Subchunk2Size);
}

#endif
