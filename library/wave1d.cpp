#include "wave1d.h"

const char* DIMENSION_DISAGREEMENT = "Dimensions of waves must agree in order to perform arithmetic on waves.";
const char* CANNOT_CONSTRUCT_FROM_WAVE1D = "You shouldn't be constructing a TWave1D or FWave1D from a generic Wave1D. Either use the proper T/F type initially, or if you think you know what you're doing, begin with an Array1D.";
const char* FOURIER_MUST_BE_POWER_OF_TWO = "Waves in Fourier space must have their dimension a power of two. You might pad the rest with zeros.";
const char* INCOMPATIBLE_TWAVE_FWAVE = "TWaves and FWaves are incompatible types. If you wish to convert between them, use the assignment operator (=).";

/**
When performing FFT's, we should pad the data in some fashion up to some
power of two. nextPowerOfTwo() finds the largest power of two => than n.
*/
int nextPowerOfTwo(int n)
{
    if(n==0) return 0;
    int q = 1;
    for(; q<n; q*=2);
    return q;
}
