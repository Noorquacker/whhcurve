#ifndef _TEXTWAVEDATA_H
#define _TEXTWAVEDATA_H

using namespace std;

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>

#include "wavedata.h"

#include "tnt/tnt.h"
using namespace TNT;

// TO-DO/BUGS LIST:
//  in open(fn), test for a bad file opening and chuck an error out
//  at the end of readText, iterate over that getline() until we make sure we 
//   get the whole thing out
//  needs a copy constructor
//  need to subclass this to read wave names
//  needs assignment/constructor from an array
//  needs ability to WRITE files
//  perhaps need to break this apart into a general superclass, binary, and text
//  bug: need bounds checking in getDataColumn


// Text load loads data1d into memory one line/column at a time using the
// extraction >> operator on a dereferenced pointer. Write your own
// custom extractor for class T, and the possibilities are endless.

template <class T>
class TextWaveData : public WaveData<T>
{
    public:
        TextWaveData() : WaveData<T>::WaveData() { create(); }
        TextWaveData(const char* fn) : WaveData<T>::WaveData(fn) { create(); }
        TextWaveData(istream &xpin) : WaveData<T>::WaveData(xpin) { create(); }
        TextWaveData(int cols, int rows)
            : WaveData<T>::WaveData(cols,rows) { create(); }
        
        void read(int cols = -1, int maxlength = -1);
        void write();
        
    protected:
        void create() { }
    
};

template <class T>
void TextWaveData<T>::read(int cols, int maxlength)
{
    if(this->pin==0)
        throw NO_INPUT_STREAM;
	
	istringstream sin;
	string *filein = new string;
	string linein;
	while(!((*(this->pin)).eof()) && !((*(this->pin)).bad()))
	{
		getline(*(this->pin), linein);
		filein->append(linein);
		filein->append("\n");
	}
	sin.str(*filein);
	delete filein;
    
    // Notice that readText starts out a lot like readBinary. The game plan
    // is to read the entire file at once, figure out the column structure,
    // then istringstream it into the data array. As with readBinary,
    // the read begins at the current file position, which may be set before
    // begining the read
    int start = sin.tellg();
    int length;

    char c;
	// get length of file:
    sin.seekg(0, ios::end);
    length = sin.tellg();    // this line seperated from the next b/c ISO C++
    length -= start;            // complains about an operator-() type conversion
    sin.seekg(start, ios::beg);

    // Ugly, but easier than fishing around through regex. Dissect the first line.
    c = sin.get();  // character iterator
    bool in_whitespace = false;
    int column_count = 0;
    
    // loop until we hit the end of the line or the end of the string.
    for(; c != '\r' && c != '\n' && c != '\0' && !sin.eof(); c = sin.get())
    {
        // check for whitespace
        switch(c)
        {
            case '\t':
            case ' ':
                if(!in_whitespace)
                {
                    column_count++;
                    in_whitespace = true;
                }
                break;

            default:
                in_whitespace = false;
                break;
        }
    }
    if(!in_whitespace)  // correction applies if newline follows data directly
        column_count++;
        
    // go back to start
    sin.seekg(start, ios::beg);

    // We now have column_count columns, and the application wants cols
    // of them (or "all" if cols == -1). We also have row_count rows,
    // and the app wants maxlength of them (or "all"/-1)
    
    // Rows and columns get treated differently. We have to make sure to throw
    // out any unused columns, but we'll stop the iteration before we even
    // get to any extra rows.
    
    if(cols == -1 || cols > column_count) cols = column_count;

    bool hasHeadings = false;
    c = sin.get();
    switch(c)
	{
		case '1':
		case '2':
		case '3':
		case '4':
		case '5':
		case '6':
		case '7':
		case '8':
		case '9':
		case '0':
		case '.':
		case '+':
		case '-':
			break;
		default:
			hasHeadings = true;
	}
    // return to start
    sin.seekg(start, ios::beg);

    // so far not DOING anything with the heading row, just getting it out of the way
    string headingRow;
    if(hasHeadings)
        getline(sin, headingRow);
	
//	cerr << "okay, starting data portion. hasHeadings == " << hasHeadings << ", start == " << start << ", cols == " << cols << endl;

    // some iterators we're about to use
    int i,j=0;
    T x;
    Array1D< vector<T> > tempcolumns(cols);

    for(i=0; (i < maxlength || maxlength==-1) && !(sin.eof()) && (sin.good()); i++)
    {
        // read the columns we want
        for(j=0; j < cols && (!sin.eof()) && (sin.good()); j++)
        {
            sin >> x;
//			cout << "reading " << x << ", eof() == " << sin.eof() << ", good() == " << sin.good() << endl;
			if(sin.good())
			{
				tempcolumns[j].push_back(x);
			}
			else
				j--;
		}
		
        // if there's more, read the rest of the line (which we don't want)
        if(column_count - cols)
        {
            for(c = 'a'; c != '\r' && c != '\n' && c != '\0' && !(sin.eof()); c = sin.get());
        }
    }
    if(j < cols)
        i--;    // didn't get all the data from the row, so put last one back

    // data2d points to an array of pointers that point to
    // the head of each column.
    (this->data2d) = Array2D<T>(cols,i);
    (this->data1d) = Array1D<T>((this->data2d).dim2(),(this->data2d)[0]);

//cerr << "data1d dim n=" << (this->data1d).dim() << endl;
//cerr << "data2d dim n1=" << (this->data2d).dim1() << ", n2=" << (this->data2d).dim2() << endl;
    // copy data in
    for(i=0; i<(this->data2d).dim1(); i++)
        for(j=0; j<(this->data2d).dim2(); j++)
            (this->data2d)[i][j] = tempcolumns[i][j];
}

template <class T>
void TextWaveData<T>::write()
{
    if(this->pout==0)
        throw NO_OUTPUT_STREAM;    
    
    for(int j=0; j < (this->data2d).dim2(); j++)
    {
        for(int i=0; i < (this->data2d).dim1(); i++)
        {
            if(i>0)
                (*(this->pout)) << "\t";    
            (*(this->pout)) << (this->data2d)[i][j];
        }
        (*(this->pout)) << endl;
    }
}

#endif // _WAVEDATA_H
