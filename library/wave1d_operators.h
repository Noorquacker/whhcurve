#ifndef _WAVE1D_OPERATORS
#define _WAVE1D_OPERATORS

#include "twave1d.h"
#include "fwave1d.h"

// define a bunch of wave operators for type double. each set of operators 
// performs an operation according to the types of its operands

/*       const               Wave1D               TWave1D                 FWave1D
const    built-in            member-by-member     member-by-member        +- mag, * / complex
Wave1D   member-by-member    each (check dim)     each (check dim)        half-width+1 (mag)
TWave1D  member-by-member    each (check dim)     each (check dim/time)   NO OP (exception)
FWave1D  +- mag, * / each     half-width+1 (mag)   NO OP (exception)       complex
*/

//      c   W   T   F
// c    cc  cW  cT  cF
// W    Wc  WW  WT  WF
// T    Tc  TW  TT  TF
// F    Fc  FW  FT  FF

// operations needed
// + - * /
// += -= *= /=

#define WAVE1D Wave1D<T,TimeType>
#define TWAVE1D TWave1D<T,TimeType>
#define FWAVE1D FWave1D<T,TimeType>


// conventions for defines: first variable a, scond variable b
// primary variable is the one used as template for return
#define DEFINE_OP_COPY(RET_TYPE,OP,TYPE1,TYPE2,PRI_VARIABLE) \
template<class T,class TimeType> \
RET_TYPE operator OP (const TYPE1 & a, const TYPE2 & b) { \
RET_TYPE q(PRI_VARIABLE.dim()); q.copyTime(PRI_VARIABLE);

#define DEFINE_OP_INPLACE(OP,TYPE1,TYPE2) \
template<class T,class TimeType> \
TYPE1 & operator OP (TYPE1 & a, const TYPE2 & b) { \

#define DEFINE_OP_EXCEPTION(OP,TYPE1,TYPE2,EXCEPTION) \
template<class T,class TimeType> \
void operator OP (const TYPE1 & a, const TYPE2 & b) { throw EXCEPTION; }

#define CHECK_BOUNDS if(!a.isCompatible(b)) throw DIMENSION_DISAGREEMENT;

#define LOOP_EVERY_COPY(N_INC) for(int i=0; i<q.dim(); i+=N_INC) {
#define LOOP_EVERY_INPLACE(N_INC) for(int i=0; i<a.dim(); i+=N_INC) {

#define END_OP_COPY } return q; }
#define END_OP_INPLACE } return a; }

// Writing 100+ functions does not sound like fun. There's a lot of functionality
// overlap here. How do we best eliminate duplicate code (and make it easier to
// debug?
//
// Common loops:
//  Member-by-member (copy)
//  Member-by-member (in-place)
//  Member with constant (copy)
//  Member with constant (in-place)
//  Two indices mapped onto one (copy)
//  Two indices mapped onto one (in-place)
//  Two indices mapped onto constant (copy)
//  Two indices mapped onto constant (in-place)
//  One index mapped onto two (copy)
//  One index mapped onto two (in-place)


// We're going to start at the upper-left of the table, and work our way down
// and to the right

// const const (built-in, or not our concern)

// const Wave1D     (member-by-member)
// Wave1D const     (member-by-member)
DEFINE_OP_COPY(WAVE1D,+,T,WAVE1D,b) LOOP_EVERY_COPY(1) q[i]=a+b[i]; END_OP_COPY
DEFINE_OP_COPY(WAVE1D,-,T,WAVE1D,b) LOOP_EVERY_COPY(1) q[i]=a-b[i]; END_OP_COPY
DEFINE_OP_COPY(WAVE1D,*,T,WAVE1D,b) LOOP_EVERY_COPY(1) q[i]=a*b[i]; END_OP_COPY
DEFINE_OP_COPY(WAVE1D,/,T,WAVE1D,b) LOOP_EVERY_COPY(1) q[i]=a/b[i]; END_OP_COPY

DEFINE_OP_COPY(WAVE1D,+,WAVE1D,T,a) LOOP_EVERY_COPY(1) q[i]=a[i]+b; END_OP_COPY
DEFINE_OP_COPY(WAVE1D,-,WAVE1D,T,a) LOOP_EVERY_COPY(1) q[i]=a[i]-b; END_OP_COPY
DEFINE_OP_COPY(WAVE1D,*,WAVE1D,T,a) LOOP_EVERY_COPY(1) q[i]=a[i]*b; END_OP_COPY
DEFINE_OP_COPY(WAVE1D,/,WAVE1D,T,a) LOOP_EVERY_COPY(1) q[i]=a[i]/b; END_OP_COPY

DEFINE_OP_INPLACE(+=,WAVE1D,T) LOOP_EVERY_INPLACE(1) a[i]+=b; END_OP_INPLACE
DEFINE_OP_INPLACE(-=,WAVE1D,T) LOOP_EVERY_INPLACE(1) a[i]-=b; END_OP_INPLACE
DEFINE_OP_INPLACE(*=,WAVE1D,T) LOOP_EVERY_INPLACE(1) a[i]*=b; END_OP_INPLACE
DEFINE_OP_INPLACE(/=,WAVE1D,T) LOOP_EVERY_INPLACE(1) a[i]/=b; END_OP_INPLACE


// Wave1D Wave1D    each (check dim) PALOOZA
DEFINE_OP_COPY(WAVE1D,+,WAVE1D,WAVE1D,a) CHECK_BOUNDS LOOP_EVERY_COPY(1) q[i]=a[i]+b[i]; END_OP_COPY
DEFINE_OP_COPY(WAVE1D,-,WAVE1D,WAVE1D,a) CHECK_BOUNDS LOOP_EVERY_COPY(1) q[i]=a[i]-b[i]; END_OP_COPY
DEFINE_OP_COPY(WAVE1D,*,WAVE1D,WAVE1D,a) CHECK_BOUNDS LOOP_EVERY_COPY(1) q[i]=a[i]*b[i]; END_OP_COPY
DEFINE_OP_COPY(WAVE1D,/,WAVE1D,WAVE1D,a) CHECK_BOUNDS LOOP_EVERY_COPY(1) q[i]=a[i]/b[i]; END_OP_COPY

DEFINE_OP_INPLACE(+=,WAVE1D,WAVE1D) CHECK_BOUNDS LOOP_EVERY_INPLACE(1) a[i]+=b[i]; END_OP_INPLACE
DEFINE_OP_INPLACE(-=,WAVE1D,WAVE1D) CHECK_BOUNDS LOOP_EVERY_INPLACE(1) a[i]-=b[i]; END_OP_INPLACE
DEFINE_OP_INPLACE(*=,WAVE1D,WAVE1D) CHECK_BOUNDS LOOP_EVERY_INPLACE(1) a[i]*=b[i]; END_OP_INPLACE
DEFINE_OP_INPLACE(/=,WAVE1D,WAVE1D) CHECK_BOUNDS LOOP_EVERY_INPLACE(1) a[i]/=b[i]; END_OP_INPLACE

// const TWave1D    (member-by-member)
// TWave1D const    (member-by-member)
DEFINE_OP_COPY(TWAVE1D,+,T,TWAVE1D,b) LOOP_EVERY_COPY(1) q[i]=a+b[i]; END_OP_COPY
DEFINE_OP_COPY(TWAVE1D,-,T,TWAVE1D,b) LOOP_EVERY_COPY(1) q[i]=a-b[i]; END_OP_COPY
DEFINE_OP_COPY(TWAVE1D,*,T,TWAVE1D,b) LOOP_EVERY_COPY(1) q[i]=a*b[i]; END_OP_COPY
DEFINE_OP_COPY(TWAVE1D,/,T,TWAVE1D,b) LOOP_EVERY_COPY(1) q[i]=a/b[i]; END_OP_COPY

DEFINE_OP_COPY(TWAVE1D,+,TWAVE1D,T,a) LOOP_EVERY_COPY(1) q[i]=a[i]+b; END_OP_COPY
DEFINE_OP_COPY(TWAVE1D,-,TWAVE1D,T,a) LOOP_EVERY_COPY(1) q[i]=a[i]-b; END_OP_COPY
DEFINE_OP_COPY(TWAVE1D,*,TWAVE1D,T,a) LOOP_EVERY_COPY(1) q[i]=a[i]*b; END_OP_COPY
DEFINE_OP_COPY(TWAVE1D,/,TWAVE1D,T,a) LOOP_EVERY_COPY(1) q[i]=a[i]/b; END_OP_COPY

DEFINE_OP_INPLACE(+=,TWAVE1D,T) LOOP_EVERY_INPLACE(1) a[i]+=b; END_OP_INPLACE
DEFINE_OP_INPLACE(-=,TWAVE1D,T) LOOP_EVERY_INPLACE(1) a[i]-=b; END_OP_INPLACE
DEFINE_OP_INPLACE(*=,TWAVE1D,T) LOOP_EVERY_INPLACE(1) a[i]*=b; END_OP_INPLACE
DEFINE_OP_INPLACE(/=,TWAVE1D,T) LOOP_EVERY_INPLACE(1) a[i]/=b; END_OP_INPLACE

// Wave1D TWave1D   (member-by-member)
// TWave1D Wave1D   (member-by-member)
DEFINE_OP_COPY(TWAVE1D,+,WAVE1D,TWAVE1D,a) CHECK_BOUNDS LOOP_EVERY_COPY(1) q[i]=a[i]+b[i]; END_OP_COPY
DEFINE_OP_COPY(TWAVE1D,-,WAVE1D,TWAVE1D,a) CHECK_BOUNDS LOOP_EVERY_COPY(1) q[i]=a[i]-b[i]; END_OP_COPY
DEFINE_OP_COPY(TWAVE1D,*,WAVE1D,TWAVE1D,a) CHECK_BOUNDS LOOP_EVERY_COPY(1) q[i]=a[i]*b[i]; END_OP_COPY
DEFINE_OP_COPY(TWAVE1D,/,WAVE1D,TWAVE1D,a) CHECK_BOUNDS LOOP_EVERY_COPY(1) q[i]=a[i]/b[i]; END_OP_COPY

DEFINE_OP_COPY(TWAVE1D,+,TWAVE1D,WAVE1D,a) CHECK_BOUNDS LOOP_EVERY_COPY(1) q[i]=a[i]+b[i]; END_OP_COPY
DEFINE_OP_COPY(TWAVE1D,-,TWAVE1D,WAVE1D,a) CHECK_BOUNDS LOOP_EVERY_COPY(1) q[i]=a[i]-b[i]; END_OP_COPY
DEFINE_OP_COPY(TWAVE1D,*,TWAVE1D,WAVE1D,a) CHECK_BOUNDS LOOP_EVERY_COPY(1) q[i]=a[i]*b[i]; END_OP_COPY
DEFINE_OP_COPY(TWAVE1D,/,TWAVE1D,WAVE1D,a) CHECK_BOUNDS LOOP_EVERY_COPY(1) q[i]=a[i]/b[i]; END_OP_COPY

DEFINE_OP_INPLACE(+=,TWAVE1D,WAVE1D) CHECK_BOUNDS LOOP_EVERY_INPLACE(1) a[i]+=b[i]; END_OP_INPLACE
DEFINE_OP_INPLACE(-=,TWAVE1D,WAVE1D) CHECK_BOUNDS LOOP_EVERY_INPLACE(1) a[i]-=b[i]; END_OP_INPLACE
DEFINE_OP_INPLACE(*=,TWAVE1D,WAVE1D) CHECK_BOUNDS LOOP_EVERY_INPLACE(1) a[i]*=b[i]; END_OP_INPLACE
DEFINE_OP_INPLACE(/=,TWAVE1D,WAVE1D) CHECK_BOUNDS LOOP_EVERY_INPLACE(1) a[i]/=b[i]; END_OP_INPLACE

// const FWave1D    mag only (complex div)
// FWave1D const    mag only
// lots of subtleties- when the magnitude is zero, default phase is real
DEFINE_OP_COPY(FWAVE1D,+,T,FWAVE1D,b) LOOP_EVERY_COPY(2)
    T mag=sqrt(b[i]*b[i]+b[i+1]*b[i+1]);
    q[i]=(mag!=0 ? b[i]*(mag+a)/mag : a);
    q[i+1]=(mag!=0 ? b[i+1]*(mag+a)/mag : 0);
END_OP_COPY
DEFINE_OP_COPY(FWAVE1D,-,T,FWAVE1D,b) LOOP_EVERY_COPY(2)
    T mag=sqrt(b[i]*b[i]+b[i+1]*b[i+1]);
    q[i]=(mag!=0 ? b[i]*(mag-a)/mag : a);
    q[i+1]=(mag!=0 ? b[i+1]*(mag-a)/mag : 0);
END_OP_COPY
DEFINE_OP_COPY(FWAVE1D,*,T,FWAVE1D,b) LOOP_EVERY_COPY(2)
    q[i]=b[i]*a;
    q[i+1]=b[i+1]*a;
END_OP_COPY
DEFINE_OP_COPY(FWAVE1D,/,T,FWAVE1D,b) LOOP_EVERY_COPY(2)
    // really annoying complex division
    T denom = b[i]*b[i]-b[i+1]*b[i+1];
    q[i] = a*b[i]/denom;
    q[i+1] = -a*b[i+1]/denom;
END_OP_COPY

DEFINE_OP_COPY(FWAVE1D,+,FWAVE1D,T,a) LOOP_EVERY_COPY(2)
    T mag=sqrt(a[i]*a[i]+a[i+1]*a[i+1]);
    q[i]=(mag!=0 ? a[i]*(mag+b)/mag : b);
    q[i+1]=(mag!=0 ? a[i+1]*(mag+b)/mag : 0);
END_OP_COPY
DEFINE_OP_COPY(FWAVE1D,-,FWAVE1D,T,a) LOOP_EVERY_COPY(2)
    T mag=sqrt(a[i]*a[i]+a[i+1]*a[i+1]);
    q[i]=(mag!=0 ? a[i]*(mag-b)/mag : b);
    q[i+1]=(mag!=0 ? a[i+1]*(mag-b)/mag : 0);
END_OP_COPY
DEFINE_OP_COPY(FWAVE1D,*,FWAVE1D,T,a) LOOP_EVERY_COPY(2)
    q[i]=a[i]*b;
    q[i+1]=a[i+1]*b;
END_OP_COPY
DEFINE_OP_COPY(FWAVE1D,/,FWAVE1D,T,a) LOOP_EVERY_COPY(2)
    q[i]=a[i]/b;
    q[i+1]=a[i+1]/b;
END_OP_COPY

DEFINE_OP_INPLACE(+=,FWAVE1D,T) LOOP_EVERY_INPLACE(2)
    T mag=sqrt(a[i]*a[i]+a[i+1]*a[i+1]);
    a[i]=(mag!=0 ? a[i]*(mag+b)/mag : b);
    a[i+1]=(mag!=0 ? a[i+1]*(mag+b)/mag : 0);
END_OP_INPLACE
DEFINE_OP_INPLACE(-=,FWAVE1D,T) LOOP_EVERY_INPLACE(2)
    T mag=sqrt(a[i]*a[i]+a[i+1]*a[i+1]);
    a[i]=(mag!=0 ? a[i]*(mag-b)/mag : b);
    a[i+1]=(mag!=0 ? a[i+1]*(mag-b)/mag : 0);
END_OP_INPLACE
DEFINE_OP_INPLACE(*=,FWAVE1D,T) LOOP_EVERY_INPLACE(2)
    a[i]=a[i]*b;
    a[i+1]=a[i+1]*b;
END_OP_INPLACE
DEFINE_OP_INPLACE(/=,FWAVE1D,T) LOOP_EVERY_INPLACE(2)
    a[i]=a[i]/b;
    a[i+1]=a[i+1]/b;
END_OP_INPLACE

// TWave1D TWave1D  each (check dim, time)
DEFINE_OP_COPY(TWAVE1D,+,TWAVE1D,TWAVE1D,a) CHECK_BOUNDS LOOP_EVERY_COPY(1) q[i]=a[i]+b[i]; END_OP_COPY
DEFINE_OP_COPY(TWAVE1D,-,TWAVE1D,TWAVE1D,a) CHECK_BOUNDS LOOP_EVERY_COPY(1) q[i]=a[i]-b[i]; END_OP_COPY
DEFINE_OP_COPY(TWAVE1D,*,TWAVE1D,TWAVE1D,a) CHECK_BOUNDS LOOP_EVERY_COPY(1) q[i]=a[i]*b[i]; END_OP_COPY
DEFINE_OP_COPY(TWAVE1D,/,TWAVE1D,TWAVE1D,a) CHECK_BOUNDS LOOP_EVERY_COPY(1) q[i]=a[i]/b[i]; END_OP_COPY

DEFINE_OP_INPLACE(+=,TWAVE1D,TWAVE1D) CHECK_BOUNDS LOOP_EVERY_INPLACE(1) a[i]+=b[i]; END_OP_INPLACE
DEFINE_OP_INPLACE(-=,TWAVE1D,TWAVE1D) CHECK_BOUNDS LOOP_EVERY_INPLACE(1) a[i]-=b[i]; END_OP_INPLACE
DEFINE_OP_INPLACE(*=,TWAVE1D,TWAVE1D) CHECK_BOUNDS LOOP_EVERY_INPLACE(1) a[i]*=b[i]; END_OP_INPLACE
DEFINE_OP_INPLACE(/=,TWAVE1D,TWAVE1D) CHECK_BOUNDS LOOP_EVERY_INPLACE(1) a[i]/=b[i]; END_OP_INPLACE

// Wave1D FWave1D   half-width+1 (mag, complex)
// FWave1D Wave1D   half-width+1 (mag)

// TWave1D FWave1D  disallowed (exception)
// FWave1D TWave1D  disallowed (exception)
DEFINE_OP_EXCEPTION(+,TWAVE1D,TWAVE1D,INCOMPATIBLE_TWAVE_FWAVE)
DEFINE_OP_EXCEPTION(-,TWAVE1D,TWAVE1D,INCOMPATIBLE_TWAVE_FWAVE)
DEFINE_OP_EXCEPTION(*,TWAVE1D,TWAVE1D,INCOMPATIBLE_TWAVE_FWAVE)
DEFINE_OP_EXCEPTION(/,TWAVE1D,TWAVE1D,INCOMPATIBLE_TWAVE_FWAVE)

DEFINE_OP_EXCEPTION(+,FWAVE1D,TWAVE1D,INCOMPATIBLE_TWAVE_FWAVE)
DEFINE_OP_EXCEPTION(-,FWAVE1D,TWAVE1D,INCOMPATIBLE_TWAVE_FWAVE)
DEFINE_OP_EXCEPTION(*,FWAVE1D,TWAVE1D,INCOMPATIBLE_TWAVE_FWAVE)
DEFINE_OP_EXCEPTION(/,FWAVE1D,TWAVE1D,INCOMPATIBLE_TWAVE_FWAVE)

DEFINE_OP_EXCEPTION(+=,TWAVE1D,TWAVE1D,INCOMPATIBLE_TWAVE_FWAVE)
DEFINE_OP_EXCEPTION(-=,TWAVE1D,TWAVE1D,INCOMPATIBLE_TWAVE_FWAVE)
DEFINE_OP_EXCEPTION(*=,TWAVE1D,TWAVE1D,INCOMPATIBLE_TWAVE_FWAVE)
DEFINE_OP_EXCEPTION(/=,TWAVE1D,TWAVE1D,INCOMPATIBLE_TWAVE_FWAVE)

DEFINE_OP_EXCEPTION(+=,FWAVE1D,TWAVE1D,INCOMPATIBLE_TWAVE_FWAVE)
DEFINE_OP_EXCEPTION(-=,FWAVE1D,TWAVE1D,INCOMPATIBLE_TWAVE_FWAVE)
DEFINE_OP_EXCEPTION(*=,FWAVE1D,TWAVE1D,INCOMPATIBLE_TWAVE_FWAVE)
DEFINE_OP_EXCEPTION(/=,FWAVE1D,TWAVE1D,INCOMPATIBLE_TWAVE_FWAVE)

// FWave1D FWave1D  complex

#endif
