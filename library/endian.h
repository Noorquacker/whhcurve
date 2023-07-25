#ifndef MACHINE_ENDIAN_H
#define MACHINE_ENDIAN_H

const unsigned MACHINE_BIG_ENDIAN = 1;
const unsigned MACHINE_LITTLE_ENDIAN = 2;

extern const char* ERROR_UNKNOWN_ENDIAN;

unsigned machine_endian();

/**	Effectively change the endian of the target.
	Note that this function may be grossly misused, such as applying it to
	any sort of data structure larger than a primitive.
*/
template <class T>
T reverse_bytes(T data)
{
    // automatic size detection for whatever type T we have
    unsigned n = sizeof(data);
    unsigned char* p = (unsigned char*) &data;
    unsigned char temp;
    for(unsigned i=0; i < (n/2); i++)
    {
        temp = *(p+i);
        *(p+i) = *(p+n-i-1);
        *(p+n-i-1) = temp;
    }
    return data;
}

#endif
