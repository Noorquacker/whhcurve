using namespace std;

#include "endian.h"

const char* ERROR_UNKNOWN_ENDIAN = "Processor/compiler has strange byte ordering or type sizes. Do a reality check first.";


unsigned machine_endian()
{
    unsigned short test = 258;
    unsigned char *byte_one, *byte_two;
    byte_one = (unsigned char*) &test;
    byte_two = byte_one+1;
    
    if((*byte_one) == 1 && (*byte_two) == 2)
        return MACHINE_BIG_ENDIAN;
    if((*byte_one) == 2 && (*byte_two) == 1)
        return MACHINE_LITTLE_ENDIAN;
    else
        throw ERROR_UNKNOWN_ENDIAN;

    return 0;
}

